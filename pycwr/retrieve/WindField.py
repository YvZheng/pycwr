"""Single-radar wind retrieval helpers for pycwr.

This module provides two related horizontal-wind retrieval paths:

* VAD: ring-by-ring harmonic fitting on a single PPI sweep
* VVP: local moving-window least-squares retrieval on one PPI sweep

The implementations target practical single-radar workflows where a user has
already read a volume as a :class:`pycwr.core.NRadar.PRD` object and wants a
scientifically traceable wind product without leaving the pycwr API.

References
----------
* Browning, K. A., and R. Wexler, 1968:
  "The Determination of Kinematic Properties of a Wind Field Using Doppler
  Radar", Journal of Applied Meteorology.
  https://doi.org/10.1175/1520-0450(1968)007<0105:TDOKPO>2.0.CO;2
* Waldteufel, P., and H. Corbin, 1979:
  "On the Analysis of Single-Doppler Radar Data", Journal of Applied
  Meteorology, 18, 532-542.
  https://doi.org/10.1175/1520-0450(1979)018<0532:OTAOSD>2.0.CO;2
"""

import concurrent.futures
import json
import os

import numpy as np
import xarray as xr
from scipy.spatial import cKDTree

from ..configure.default_config import DEFAULT_METADATA


WIND_REFERENCE_NOTES = {
    "VAD": [
        "Browning and Wexler (1968), Journal of Applied Meteorology, doi:10.1175/1520-0450(1968)007<0105:TDOKPO>2.0.CO;2",
    ],
    "VVP": [
        "Waldteufel and Corbin (1979), Journal of Applied Meteorology, doi:10.1175/1520-0450(1979)018<0532:OTAOSD>2.0.CO;2",
    ],
    "WIND_VOLUME": [
        "Waldteufel and Corbin (1979), Journal of Applied Meteorology, doi:10.1175/1520-0450(1979)018<0532:OTAOSD>2.0.CO;2",
        "Holleman (2003), Doppler Radar Wind Profiles, KNMI publication.",
        "Holleman (2005), Quality Control and Verification of Weather Radar Wind Profiles, KNMI publication.",
    ],
}


def _as_float_array(values, fillvalue=None):
    """Return a float64 array where fill values are replaced by NaN."""
    array = np.asarray(values, dtype=np.float64)
    if fillvalue is None:
        return array
    return np.where(array == float(fillvalue), np.nan, array)


def _azimuth_coverage_deg(azimuth):
    """Return the azimuth coverage in degrees after removing the largest gap."""
    azimuth = np.asarray(azimuth, dtype=np.float64)
    azimuth = np.mod(azimuth[np.isfinite(azimuth)], 360.0)
    if azimuth.size < 2:
        return 0.0
    azimuth = np.unique(np.sort(azimuth))
    if azimuth.size < 2:
        return 0.0
    gaps = np.diff(np.concatenate((azimuth, [azimuth[0] + 360.0])))
    return float(max(0.0, 360.0 - np.max(gaps)))


def _azimuth_sector_count(azimuth, sector_width_deg=45.0):
    """Return the number of occupied azimuth sectors."""
    azimuth = np.asarray(azimuth, dtype=np.float64)
    azimuth = np.mod(azimuth[np.isfinite(azimuth)], 360.0)
    if azimuth.size == 0:
        return 0
    sector_width_deg = float(sector_width_deg)
    if sector_width_deg <= 0.0:
        raise ValueError("sector_width_deg must be positive")
    sector_index = np.floor(azimuth / sector_width_deg).astype(np.int32)
    return int(np.unique(sector_index).size)


def _effective_min_sector_count(requested_count, azimuth_window, sector_width_deg=45.0):
    """Cap a requested sector threshold by the local window capacity."""
    requested_count = int(requested_count)
    if requested_count <= 0:
        return 0
    capacity = _azimuth_sector_count(azimuth_window, sector_width_deg=sector_width_deg)
    return int(min(requested_count, max(capacity, 1)))


def _design_matrix(azimuth_deg, elevation_deg):
    """Build the single-Doppler horizontal-wind design matrix."""
    azimuth_rad = np.deg2rad(np.asarray(azimuth_deg, dtype=np.float64))
    elevation_rad = np.deg2rad(np.asarray(elevation_deg, dtype=np.float64))
    if elevation_rad.ndim == 0:
        elevation_rad = np.full(azimuth_rad.shape, float(elevation_rad), dtype=np.float64)
    cos_elevation = np.cos(elevation_rad)
    return np.column_stack(
        (
            np.sin(azimuth_rad) * cos_elevation,
            np.cos(azimuth_rad) * cos_elevation,
            np.ones(azimuth_rad.shape, dtype=np.float64),
        )
    )


def _wind_direction_from_uv(u, v):
    """Return meteorological wind direction (degrees from north, clockwise)."""
    direction = (270.0 - np.rad2deg(np.arctan2(v, u))) % 360.0
    return direction


def _solve_weighted_least_squares(design, response, weights=None, max_condition_number=1.0e6):
    """Solve a weighted least-squares system and return coefficients plus diagnostics."""
    response = np.asarray(response, dtype=np.float64).reshape(-1)
    design = np.asarray(design, dtype=np.float64)
    if weights is None:
        weights = np.ones(response.shape, dtype=np.float64)
    else:
        weights = np.asarray(weights, dtype=np.float64).reshape(-1)

    valid = (
        np.all(np.isfinite(design), axis=1)
        & np.isfinite(response)
        & np.isfinite(weights)
        & (weights > 0.0)
    )
    if not np.any(valid):
        return None

    design_valid = design[valid]
    response_valid = response[valid]
    weight_valid = weights[valid]
    def _solve_once(current_weights):
        sqrt_weight = np.sqrt(current_weights)
        design_weighted = design_valid * sqrt_weight[:, np.newaxis]
        response_weighted = response_valid * sqrt_weight
        try:
            singular_values = np.linalg.svd(design_weighted, compute_uv=False)
        except np.linalg.LinAlgError:
            singular_values = np.array([np.nan], dtype=np.float64)
        smallest_singular_value = float(np.min(singular_values)) if singular_values.size else np.nan
        normal = design_weighted.T.dot(design_weighted)
        condition_number = float(np.linalg.cond(normal))
        if (not np.isfinite(condition_number)) or condition_number > float(max_condition_number):
            return None, condition_number, smallest_singular_value
        try:
            coefficients = np.linalg.solve(normal, design_weighted.T.dot(response_weighted))
        except np.linalg.LinAlgError:
            return None, condition_number, smallest_singular_value
        return coefficients, condition_number, smallest_singular_value

    coefficients, condition_number, smallest_singular_value = _solve_once(weight_valid)
    if coefficients is None:
        return {
            "coefficients": None,
            "condition_number": condition_number,
            "smallest_singular_value": smallest_singular_value,
            "rmse": np.nan,
            "fitted": None,
            "valid_count": int(valid.sum()),
        }

    robust_weights = weight_valid.copy()
    for _ in range(5):
        fitted = design_valid.dot(coefficients)
        residual = response_valid - fitted
        scale = 1.4826 * np.median(np.abs(residual - np.median(residual)))
        if (not np.isfinite(scale)) or scale < 1.0e-6:
            break
        normalized = np.abs(residual) / (1.5 * scale)
        huber = np.where(normalized <= 1.0, 1.0, 1.0 / normalized)
        updated_weights = weight_valid * huber
        new_coefficients, new_condition_number, new_smallest_singular_value = _solve_once(updated_weights)
        if new_coefficients is None:
            break
        if np.allclose(new_coefficients, coefficients, rtol=1.0e-5, atol=1.0e-5):
            coefficients = new_coefficients
            condition_number = new_condition_number
            smallest_singular_value = new_smallest_singular_value
            robust_weights = updated_weights
            break
        coefficients = new_coefficients
        condition_number = new_condition_number
        smallest_singular_value = new_smallest_singular_value
        robust_weights = updated_weights
    fitted = design_valid.dot(coefficients)
    rmse = float(np.sqrt(np.mean((response_valid - fitted) ** 2))) if response_valid.size else np.nan
    return {
        "coefficients": coefficients,
        "condition_number": condition_number,
        "smallest_singular_value": smallest_singular_value,
        "rmse": rmse,
        "fitted": fitted,
        "valid_count": int(valid.sum()),
        "weights": robust_weights,
    }


def _empty_fit_result(valid_count, valid_fraction, azimuth_coverage_deg):
    """Return a standard empty fit record."""
    return {
        "u": np.nan,
        "v": np.nan,
        "wind_speed": np.nan,
        "wind_direction": np.nan,
        "radial_offset": np.nan,
        "fit_rmse": np.nan,
        "condition_number": np.nan,
        "smallest_singular_value": np.nan,
        "valid_count": int(valid_count),
        "valid_fraction": float(valid_fraction),
        "azimuth_coverage_deg": float(azimuth_coverage_deg),
        "azimuth_sector_count": 0,
    }


def fit_vad_ring(
    azimuth,
    elevation,
    radial_velocity,
    fillvalue=-999.0,
    min_valid_fraction=0.5,
    min_valid_count=16,
    min_azimuth_coverage_deg=240.0,
    min_azimuth_sector_count=4,
    max_condition_number=1.0e6,
    weights=None,
):
    """
    Fit a single range ring using a first-harmonic VAD model.

    Parameters
    ----------
    azimuth : array-like
        Azimuth angles in degrees.
    elevation : float or array-like
        Sweep elevation angle in degrees. A scalar is typical for one PPI.
    radial_velocity : array-like
        Radial velocity samples for one ring.
    fillvalue : float
        Missing-value marker.

    Returns
    -------
    dict
        Retrieval record containing horizontal wind, fit quality, and sample
        coverage diagnostics.
    """
    azimuth = np.asarray(azimuth, dtype=np.float64).reshape(-1)
    radial_velocity = _as_float_array(radial_velocity, fillvalue=fillvalue).reshape(-1)
    if azimuth.shape != radial_velocity.shape:
        raise ValueError("azimuth and radial_velocity must have the same shape")
    if np.asarray(elevation).ndim == 0:
        elevation = np.full(azimuth.shape, float(elevation), dtype=np.float64)
    else:
        elevation = np.asarray(elevation, dtype=np.float64).reshape(-1)
    if elevation.shape != azimuth.shape:
        raise ValueError("elevation must be scalar or match the azimuth shape")

    valid = np.isfinite(radial_velocity)
    valid_count = int(valid.sum())
    valid_fraction = float(valid_count / radial_velocity.size) if radial_velocity.size else 0.0
    azimuth_coverage = _azimuth_coverage_deg(azimuth[valid])
    azimuth_sector_count = _azimuth_sector_count(azimuth[valid])
    if (
        valid_count < int(min_valid_count)
        or valid_fraction < float(min_valid_fraction)
        or azimuth_coverage < float(min_azimuth_coverage_deg)
        or azimuth_sector_count < int(min_azimuth_sector_count)
    ):
        result = _empty_fit_result(valid_count, valid_fraction, azimuth_coverage)
        result["azimuth_sector_count"] = int(azimuth_sector_count)
        return result

    fit = _solve_weighted_least_squares(
        _design_matrix(azimuth, elevation),
        radial_velocity,
        weights=weights,
        max_condition_number=max_condition_number,
    )
    if fit is None or fit["coefficients"] is None:
        result = _empty_fit_result(valid_count, valid_fraction, azimuth_coverage)
        if fit is not None:
            result["condition_number"] = fit["condition_number"]
            result["smallest_singular_value"] = fit.get("smallest_singular_value", np.nan)
        result["azimuth_sector_count"] = int(azimuth_sector_count)
        return result

    u = float(fit["coefficients"][0])
    v = float(fit["coefficients"][1])
    return {
        "u": u,
        "v": v,
        "wind_speed": float(np.hypot(u, v)),
        "wind_direction": float(_wind_direction_from_uv(u, v)),
        "radial_offset": float(fit["coefficients"][2]),
        "fit_rmse": float(fit["rmse"]),
        "condition_number": float(fit["condition_number"]),
        "smallest_singular_value": float(fit.get("smallest_singular_value", np.nan)),
        "valid_count": valid_count,
        "valid_fraction": valid_fraction,
        "azimuth_coverage_deg": azimuth_coverage,
        "azimuth_sector_count": int(azimuth_sector_count),
    }


def _resolve_velocity_field_name(prd, sweep, field_name=None):
    """Choose a velocity field for one sweep, preferring corrected velocity."""
    available = set(prd.available_fields(sweep=int(sweep), range_mode="aligned"))
    if field_name is not None:
        if field_name not in available:
            raise KeyError("field %s is not available on sweep %s" % (field_name, sweep))
        return str(field_name)
    for candidate in ("Vc", "V"):
        if candidate in available:
            return candidate
    raise KeyError("No usable velocity field is available on sweep %s" % sweep)


def select_velocity_field(prd, sweep, field_name=None):
    """Public helper that resolves the velocity field used for retrieval."""
    return _resolve_velocity_field_name(prd, sweep=sweep, field_name=field_name)


def _init_vad_output(shape):
    """Allocate padded output arrays for VAD retrieval."""
    arrays = {}
    for name in (
        "u",
        "v",
        "wind_speed",
        "wind_direction",
        "radial_offset",
        "fit_rmse",
        "condition_number",
        "range",
        "height",
        "valid_fraction",
        "azimuth_coverage_deg",
    ):
        arrays[name] = np.full(shape, np.nan, dtype=np.float64)
    arrays["valid_count"] = np.zeros(shape, dtype=np.int32)
    return arrays


def _weighted_quantile(values, quantile, weights):
    """Compute a weighted quantile for one-dimensional finite data."""
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    weights = np.asarray(weights, dtype=np.float64).reshape(-1)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(valid):
        return np.nan
    values = values[valid]
    weights = weights[valid]
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    if cumulative[-1] <= 0.0:
        return np.nan
    target = float(quantile) * cumulative[-1]
    return float(np.interp(target, cumulative, values))


def _robust_weighted_location(values, weights):
    """Return a robust weighted location estimate for one-dimensional data."""
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    weights = np.asarray(weights, dtype=np.float64).reshape(-1)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(valid):
        return np.nan
    values = values[valid]
    weights = weights[valid]
    center = _weighted_quantile(values, 0.5, weights)
    spread = _weighted_quantile(np.abs(values - center), 0.5, weights)
    if (not np.isfinite(spread)) or spread < 1.0e-6:
        return center
    normalized = np.abs(values - center) / (3.0 * 1.4826 * spread)
    robust = np.where(normalized <= 1.0, (1.0 - normalized ** 2) ** 2, 0.0)
    total_weights = weights * robust
    if not np.any(total_weights > 0.0):
        return center
    return float(np.sum(values * total_weights) / np.sum(total_weights))


def _vertical_smooth(values, weights, window=3):
    """Apply a short weighted vertical smoother to a 1-D profile."""
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    weights = np.asarray(weights, dtype=np.float64).reshape(-1)
    if window <= 1:
        return values.copy()
    half_window = int(window) // 2
    smoothed = values.copy()
    for index in range(values.size):
        start = max(0, index - half_window)
        end = min(values.size, index + half_window + 1)
        smoothed[index] = _robust_weighted_location(values[start:end], weights[start:end])
    return smoothed


def _interpolate_short_gaps(values, max_gap_bins=2):
    """Linearly fill short interior gaps of a one-dimensional profile."""
    values = np.asarray(values, dtype=np.float64).copy()
    finite = np.isfinite(values)
    if finite.sum() < 2:
        return values
    index = np.arange(values.size, dtype=np.int32)
    gap_start = None
    for idx in range(values.size):
        if not finite[idx] and gap_start is None:
            gap_start = idx
        if finite[idx] and gap_start is not None:
            gap_length = idx - gap_start
            if gap_start > 0 and gap_length <= int(max_gap_bins):
                left = gap_start - 1
                right = idx
                values[gap_start:right] = np.interp(index[gap_start:right], [left, right], [values[left], values[right]])
            gap_start = None
    return values


def _gate_heights_from_field(field):
    """Return one representative height per range gate."""
    if "z" not in field.coords:
        return np.full(field.range.size, np.nan, dtype=np.float64)
    z_values = np.asarray(field["z"].values, dtype=np.float64)
    if z_values.ndim == 1:
        return z_values
    return np.nanmean(z_values, axis=0)


def _apply_range_limit(field, max_range_km=None):
    """Optionally clip a sweep field to a requested maximum range."""
    if max_range_km is None:
        return field
    max_range_m = float(max_range_km) * 1000.0
    return field.sel(range=field.range.values <= max_range_m)


def _as_1d_float(values):
    """Return a one-dimensional float64 array."""
    return np.asarray(values, dtype=np.float64).reshape(-1)


def _normalize_level_heights(level_heights):
    """Validate and normalize requested target heights."""
    level_heights = _as_1d_float(level_heights)
    if level_heights.size == 0:
        raise ValueError("level_heights must contain at least one value")
    if not np.all(np.isfinite(level_heights)):
        raise ValueError("level_heights must be finite")
    if level_heights.size > 1 and np.any(np.diff(level_heights) <= 0.0):
        raise ValueError("level_heights must be strictly increasing")
    return level_heights


def _infer_grid_spacing(axis):
    """Infer the representative spacing of a one-dimensional grid axis."""
    axis = _as_1d_float(axis)
    if axis.size < 2:
        return np.nan
    diffs = np.abs(np.diff(axis))
    finite = diffs[np.isfinite(diffs) & (diffs > 0.0)]
    if finite.size == 0:
        return np.nan
    return float(np.nanmedian(finite))


def _default_horizontal_radius_m(x_axis, y_axis):
    """Choose a conservative horizontal search radius from grid spacing."""
    x_spacing = _infer_grid_spacing(x_axis)
    y_spacing = _infer_grid_spacing(y_axis)
    spacing = np.nanmax([x_spacing, y_spacing])
    if not np.isfinite(spacing) or spacing <= 0.0:
        spacing = 2000.0
    return float(max(1500.0, 2.5 * spacing))


def _default_vertical_tolerance_m(level_heights):
    """Choose a conservative vertical tolerance from target level spacing."""
    level_heights = _normalize_level_heights(level_heights)
    if level_heights.size < 2:
        return 500.0
    spacing = np.nanmedian(np.diff(level_heights))
    if not np.isfinite(spacing) or spacing <= 0.0:
        return 500.0
    return float(max(250.0, 0.5 * spacing))


def _default_max_vertical_gap_m(level_heights):
    """Choose a conservative maximum vertical interpolation gap."""
    level_heights = _normalize_level_heights(level_heights)
    if level_heights.size < 2:
        return 1000.0
    spacing = np.nanmedian(np.diff(level_heights))
    if not np.isfinite(spacing) or spacing <= 0.0:
        return 1000.0
    return float(max(1000.0, 2.5 * spacing))


def _resolve_worker_count(workers, task_count):
    """Resolve the effective worker count for a parallel stage."""
    if workers is None:
        return 1
    workers = int(workers)
    if workers <= 0:
        raise ValueError("workers must be a positive integer")
    available = max(1, int(os.cpu_count() or 1))
    if task_count < 2:
        return 1
    return max(1, min(workers, available, int(task_count)))


def _resolve_requested_sweeps(prd, sweeps):
    """Normalize the requested sweep selection mode for wind-volume retrieval."""
    if sweeps is None:
        return "auto", [int(sweep) for sweep in np.asarray(prd.scan_info.sweep.values, dtype=np.int32)]
    if isinstance(sweeps, str):
        lowered = sweeps.strip().lower()
        if lowered == "all":
            return "all", [int(sweep) for sweep in np.asarray(prd.scan_info.sweep.values, dtype=np.int32)]
        if lowered == "auto":
            return "auto", [int(sweep) for sweep in np.asarray(prd.scan_info.sweep.values, dtype=np.int32)]
        raise ValueError("sweeps string mode must be 'auto' or 'all'")
    if np.isscalar(sweeps):
        return "explicit", [int(sweeps)]
    return "explicit", [int(sweep) for sweep in sweeps]


def _summary_jsonable(summary):
    """Convert one sweep summary dict into JSON-friendly values."""
    clean = {}
    for key, value in summary.items():
        if isinstance(value, np.generic):
            value = value.item()
        if isinstance(value, float) and not np.isfinite(value):
            value = None
        clean[key] = value
    return clean


def _extract_result_array(retrieval, key):
    """Extract an ndarray from either xarray or dict-like retrieval results."""
    value = retrieval[key]
    if hasattr(value, "values"):
        value = value.values
    return np.asarray(value, dtype=np.float64)


def _extract_result_int_array(retrieval, key):
    """Extract an int ndarray from either xarray or dict-like retrieval results."""
    value = retrieval[key]
    if hasattr(value, "values"):
        value = value.values
    return np.asarray(value)


def _flatten_valid_vvp_samples(retrieval):
    """Return flattened valid VVP samples with coordinates and diagnostics."""
    x = _extract_result_array(retrieval, "x").reshape(-1)
    y = _extract_result_array(retrieval, "y").reshape(-1)
    z = _extract_result_array(retrieval, "z").reshape(-1)
    u = _extract_result_array(retrieval, "u").reshape(-1)
    v = _extract_result_array(retrieval, "v").reshape(-1)
    fit_rmse = _extract_result_array(retrieval, "fit_rmse").reshape(-1)
    valid_count = _extract_result_array(retrieval, "valid_count").reshape(-1)
    valid_fraction = _extract_result_array(retrieval, "valid_fraction").reshape(-1)
    condition_number = _extract_result_array(retrieval, "condition_number").reshape(-1)
    azimuth_sector_count = _extract_result_array(retrieval, "azimuth_sector_count").reshape(-1)
    mask = (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(z)
        & np.isfinite(u)
        & np.isfinite(v)
        & np.isfinite(fit_rmse)
        & np.isfinite(valid_count)
        & np.isfinite(valid_fraction)
    )
    return {
        "x": x[mask],
        "y": y[mask],
        "z": z[mask],
        "u": u[mask],
        "v": v[mask],
        "fit_rmse": fit_rmse[mask],
        "valid_count": valid_count[mask],
        "valid_fraction": valid_fraction[mask],
        "condition_number": condition_number[mask],
        "azimuth_sector_count": azimuth_sector_count[mask],
    }


def _aggregate_neighbors(samples, neighbors, target_xy, power=2.0, expanded_radius=False):
    """Aggregate local VVP neighbors onto one target horizontal point."""
    if len(neighbors) == 0:
        return None
    indices = np.asarray(neighbors, dtype=np.int32)
    dx = samples["x"][indices] - float(target_xy[0])
    dy = samples["y"][indices] - float(target_xy[1])
    distance = np.hypot(dx, dy)
    quality = (
        np.maximum(samples["valid_count"][indices], 1.0)
        * np.clip(samples["valid_fraction"][indices], 0.0, 1.0)
        / np.maximum(samples["fit_rmse"][indices], 0.5)
    )
    distance = np.maximum(distance, 1.0)
    weights = quality / np.power(distance, float(power))
    weights = np.where(np.isfinite(weights) & (weights > 0.0), weights, 0.0)
    if not np.any(weights > 0.0):
        return None
    total = float(np.sum(weights))
    return {
        "u": float(np.sum(samples["u"][indices] * weights) / total),
        "v": float(np.sum(samples["v"][indices] * weights) / total),
        "z": float(np.sum(samples["z"][indices] * weights) / total),
        "fit_rmse": float(np.sum(samples["fit_rmse"][indices] * weights) / total),
        "valid_count": int(np.round(np.sum(samples["valid_count"][indices] * weights) / total)),
        "neighbor_count": int(indices.size),
        "condition_number": float(np.sum(samples["condition_number"][indices] * weights) / total),
        "azimuth_sector_count": int(np.round(np.sum(samples["azimuth_sector_count"][indices] * weights) / total)),
        "expanded_radius": bool(expanded_radius),
    }


def _grid_vvp_retrieval_xy(
    retrieval,
    target_x,
    target_y,
    horizontal_radius_m,
    horizontal_min_neighbors=3,
    max_horizontal_radius_m=None,
    distance_power=2.0,
):
    """Grid one sweep VVP analysis onto the target horizontal grid."""
    target_x = np.asarray(target_x, dtype=np.float64)
    target_y = np.asarray(target_y, dtype=np.float64)
    if target_x.shape != target_y.shape:
        raise ValueError("target_x and target_y must have the same shape")
    shape = target_x.shape
    outputs = {
        "u": np.full(shape, np.nan, dtype=np.float64),
        "v": np.full(shape, np.nan, dtype=np.float64),
        "z": np.full(shape, np.nan, dtype=np.float64),
        "fit_rmse": np.full(shape, np.nan, dtype=np.float64),
        "valid_count": np.zeros(shape, dtype=np.int32),
        "neighbor_count": np.zeros(shape, dtype=np.int32),
        "condition_number": np.full(shape, np.nan, dtype=np.float64),
        "azimuth_sector_count": np.zeros(shape, dtype=np.int32),
        "expanded_radius": np.zeros(shape, dtype=np.uint8),
    }
    samples = _flatten_valid_vvp_samples(retrieval)
    if samples["x"].size == 0:
        return outputs
    points = np.column_stack((samples["x"], samples["y"]))
    tree = cKDTree(points)
    target_points = np.column_stack((target_x.reshape(-1), target_y.reshape(-1)))
    initial_radius_m = float(horizontal_radius_m)
    if max_horizontal_radius_m is None:
        max_horizontal_radius_m = max(initial_radius_m, 2.0 * initial_radius_m)
    max_horizontal_radius_m = float(max_horizontal_radius_m)
    neighborhoods = tree.query_ball_point(target_points, r=initial_radius_m)
    for flat_index, neighbors in enumerate(neighborhoods):
        expanded_radius = False
        if len(neighbors) < int(horizontal_min_neighbors) and max_horizontal_radius_m > initial_radius_m:
            neighbors = tree.query_ball_point(target_points[flat_index], r=max_horizontal_radius_m)
            expanded_radius = True
        if len(neighbors) < int(horizontal_min_neighbors):
            continue
        aggregated = _aggregate_neighbors(
            samples,
            neighbors,
            target_points[flat_index],
            power=distance_power,
            expanded_radius=expanded_radius,
        )
        if aggregated is None:
            continue
        x_index, y_index = np.unravel_index(flat_index, shape)
        for name in ("u", "v", "z", "fit_rmse", "condition_number"):
            outputs[name][x_index, y_index] = aggregated[name]
        outputs["valid_count"][x_index, y_index] = aggregated["valid_count"]
        outputs["neighbor_count"][x_index, y_index] = aggregated["neighbor_count"]
        outputs["azimuth_sector_count"][x_index, y_index] = aggregated["azimuth_sector_count"]
        outputs["expanded_radius"][x_index, y_index] = 1 if aggregated["expanded_radius"] else 0
    return outputs


def _vertical_profile_from_sweeps(
    heights,
    u_values,
    v_values,
    rmse_values,
    count_values,
    neighbor_counts,
    expanded_radius_flags,
    level_heights,
    vertical_tolerance_m,
    max_vertical_gap_m,
):
    """Reconstruct one target-grid wind profile from sweep-level samples."""
    level_heights = _normalize_level_heights(level_heights)
    output = {
        "u": np.full(level_heights.shape, np.nan, dtype=np.float64),
        "v": np.full(level_heights.shape, np.nan, dtype=np.float64),
        "fit_rmse": np.full(level_heights.shape, np.nan, dtype=np.float64),
        "valid_count": np.zeros(level_heights.shape, dtype=np.int32),
        "source_sweep_count": np.zeros(level_heights.shape, dtype=np.int32),
        "neighbor_count": np.zeros(level_heights.shape, dtype=np.int32),
        "effective_vertical_support": np.zeros(level_heights.shape, dtype=np.int32),
        "sweep_spread_m": np.full(level_heights.shape, np.nan, dtype=np.float64),
        "quality_score": np.zeros(level_heights.shape, dtype=np.float64),
        "quality_flag": np.zeros(level_heights.shape, dtype=np.uint8),
    }
    valid = (
        np.isfinite(heights)
        & np.isfinite(u_values)
        & np.isfinite(v_values)
        & np.isfinite(rmse_values)
        & np.isfinite(count_values)
        & np.isfinite(neighbor_counts)
        & np.isfinite(expanded_radius_flags)
        & (count_values > 0)
    )
    if not np.any(valid):
        return output
    heights = np.asarray(heights[valid], dtype=np.float64)
    u_values = np.asarray(u_values[valid], dtype=np.float64)
    v_values = np.asarray(v_values[valid], dtype=np.float64)
    rmse_values = np.asarray(rmse_values[valid], dtype=np.float64)
    count_values = np.asarray(count_values[valid], dtype=np.float64)
    neighbor_counts = np.asarray(neighbor_counts[valid], dtype=np.float64)
    expanded_radius_flags = np.asarray(expanded_radius_flags[valid], dtype=np.float64)
    order = np.argsort(heights)
    heights = heights[order]
    u_values = u_values[order]
    v_values = v_values[order]
    rmse_values = rmse_values[order]
    count_values = count_values[order]
    neighbor_counts = neighbor_counts[order]
    expanded_radius_flags = expanded_radius_flags[order]
    for level_index, level_height in enumerate(level_heights):
        score = 0.0
        flag = 0
        near = np.abs(heights - level_height) <= float(vertical_tolerance_m)
        if np.any(near):
            distance = np.maximum(np.abs(heights[near] - level_height), 1.0)
            weights = count_values[near] / (np.maximum(rmse_values[near], 0.5) * distance)
            weights = np.where(np.isfinite(weights) & (weights > 0.0), weights, 0.0)
            if np.any(weights > 0.0):
                spread = float(np.max(heights[near]) - np.min(heights[near])) if np.sum(near) > 1 else 0.0
                expanded = bool(np.any(expanded_radius_flags[near] > 0.5))
                output["u"][level_index] = float(np.sum(u_values[near] * weights) / np.sum(weights))
                output["v"][level_index] = float(np.sum(v_values[near] * weights) / np.sum(weights))
                output["fit_rmse"][level_index] = float(np.sum(rmse_values[near] * weights) / np.sum(weights))
                output["valid_count"][level_index] = int(np.round(np.sum(count_values[near] * weights) / np.sum(weights)))
                output["source_sweep_count"][level_index] = int(np.sum(near))
                output["neighbor_count"][level_index] = int(np.round(np.sum(neighbor_counts[near] * weights) / np.sum(weights)))
                output["effective_vertical_support"][level_index] = int(np.sum(near))
                output["sweep_spread_m"][level_index] = spread
                score, flag = _compute_voxel_quality(
                    source_sweep_count=output["source_sweep_count"][level_index],
                    valid_count=output["valid_count"][level_index],
                    fit_rmse=output["fit_rmse"][level_index],
                    neighbor_count=output["neighbor_count"][level_index],
                    sweep_spread_m=spread,
                    max_vertical_gap_m=max_vertical_gap_m,
                    expanded_radius=expanded,
                )
                output["quality_score"][level_index] = score
                output["quality_flag"][level_index] = flag
                continue
        below = np.where(heights < level_height)[0]
        above = np.where(heights > level_height)[0]
        if below.size == 0 or above.size == 0:
            continue
        lower = below[-1]
        upper = above[0]
        span = heights[upper] - heights[lower]
        if (not np.isfinite(span)) or span <= 0.0 or span > float(max_vertical_gap_m):
            continue
        fraction = float((level_height - heights[lower]) / span)
        expanded = bool(expanded_radius_flags[lower] > 0.5 or expanded_radius_flags[upper] > 0.5)
        output["u"][level_index] = float((1.0 - fraction) * u_values[lower] + fraction * u_values[upper])
        output["v"][level_index] = float((1.0 - fraction) * v_values[lower] + fraction * v_values[upper])
        output["fit_rmse"][level_index] = float((1.0 - fraction) * rmse_values[lower] + fraction * rmse_values[upper])
        output["valid_count"][level_index] = int(
            np.round((1.0 - fraction) * count_values[lower] + fraction * count_values[upper])
        )
        output["source_sweep_count"][level_index] = 2
        output["neighbor_count"][level_index] = int(
            np.round((1.0 - fraction) * neighbor_counts[lower] + fraction * neighbor_counts[upper])
        )
        output["effective_vertical_support"][level_index] = 2
        output["sweep_spread_m"][level_index] = float(span)
        score, flag = _compute_voxel_quality(
            source_sweep_count=2,
            valid_count=output["valid_count"][level_index],
            fit_rmse=output["fit_rmse"][level_index],
            neighbor_count=output["neighbor_count"][level_index],
            sweep_spread_m=float(span),
            max_vertical_gap_m=max_vertical_gap_m,
            expanded_radius=expanded,
        )
        output["quality_score"][level_index] = score
        output["quality_flag"][level_index] = flag
    return output


def _compute_voxel_quality(
    source_sweep_count,
    valid_count,
    fit_rmse,
    neighbor_count,
    sweep_spread_m,
    max_vertical_gap_m,
    expanded_radius=False,
):
    """Return a simple quality score/flag pair for one reconstructed voxel."""
    if source_sweep_count <= 0 or valid_count <= 0 or not np.isfinite(fit_rmse):
        return 0.0, np.uint8(0)
    support_component = min(1.0, float(source_sweep_count) / 3.0)
    sample_component = min(1.0, float(valid_count) / 60.0)
    neighbor_component = min(1.0, float(max(neighbor_count, 0)) / 6.0)
    rmse_component = float(np.exp(-max(float(fit_rmse), 0.0) / 6.0))
    if np.isfinite(sweep_spread_m) and max_vertical_gap_m > 0.0:
        gap_component = max(0.0, 1.0 - float(sweep_spread_m) / float(max_vertical_gap_m))
    else:
        gap_component = 1.0
    score = 100.0 * (
        0.30 * support_component
        + 0.20 * sample_component
        + 0.20 * neighbor_component
        + 0.20 * rmse_component
        + 0.10 * gap_component
    )
    if expanded_radius:
        score *= 0.9
    score = float(np.clip(score, 0.0, 100.0))
    if score <= 0.0:
        flag = 0
    elif score < 40.0:
        flag = 1
    elif score < 70.0:
        flag = 2
    else:
        flag = 3
    return score, np.uint8(flag)


def _build_wind_volume_dataset_xy(
    x_axis,
    y_axis,
    level_heights,
    u_volume,
    v_volume,
    fit_rmse_volume,
    valid_count_volume,
    source_sweep_count_volume,
    neighbor_count_volume,
    effective_vertical_support_volume,
    sweep_spread_m_volume,
    quality_score_volume,
    quality_flag_volume,
    attrs,
):
    """Build a standard xy wind-volume dataset."""
    wind_speed = np.hypot(u_volume, v_volume)
    wind_direction = _wind_direction_from_uv(u_volume, v_volume)
    dataset = xr.Dataset(
        data_vars={
            "u": (("z_wind", "x_wind", "y_wind"), u_volume),
            "v": (("z_wind", "x_wind", "y_wind"), v_volume),
            "wind_speed": (("z_wind", "x_wind", "y_wind"), wind_speed),
            "wind_direction": (("z_wind", "x_wind", "y_wind"), wind_direction),
            "fit_rmse": (("z_wind", "x_wind", "y_wind"), fit_rmse_volume),
            "valid_count": (("z_wind", "x_wind", "y_wind"), valid_count_volume),
            "source_sweep_count": (("z_wind", "x_wind", "y_wind"), source_sweep_count_volume),
            "neighbor_count": (("z_wind", "x_wind", "y_wind"), neighbor_count_volume),
            "effective_vertical_support": (
                ("z_wind", "x_wind", "y_wind"),
                effective_vertical_support_volume,
            ),
            "sweep_spread_m": (("z_wind", "x_wind", "y_wind"), sweep_spread_m_volume),
            "quality_score": (("z_wind", "x_wind", "y_wind"), quality_score_volume),
            "quality_flag": (("z_wind", "x_wind", "y_wind"), quality_flag_volume),
        },
        coords={
            "z_wind": level_heights,
            "x_wind": x_axis,
            "y_wind": y_axis,
        },
        attrs=dict(attrs),
    )
    dataset["z_wind"].attrs = {
        "units": "meters",
        "standard_name": "height",
        "long_name": "target wind-volume height above mean sea level",
        "positive": "up",
    }
    dataset["x_wind"].attrs = {
        "units": "meters",
        "standard_name": "wind_volume_x_axis",
        "long_name": "east distance from radar",
        "axis": "xy_coordinate",
    }
    dataset["y_wind"].attrs = {
        "units": "meters",
        "standard_name": "wind_volume_y_axis",
        "long_name": "north distance from radar",
        "axis": "xy_coordinate",
    }
    dataset["u"].attrs = {
        "units": "meters_per_second",
        "standard_name": "eastward_wind_component",
        "long_name": "gridded eastward wind component from single-radar VVP",
    }
    dataset["v"].attrs = {
        "units": "meters_per_second",
        "standard_name": "northward_wind_component",
        "long_name": "gridded northward wind component from single-radar VVP",
    }
    dataset["wind_speed"].attrs = {
        "units": "meters_per_second",
        "long_name": "gridded horizontal wind speed",
    }
    dataset["wind_direction"].attrs = {
        "units": "degrees",
        "long_name": "gridded meteorological wind direction",
    }
    dataset["fit_rmse"].attrs = {
        "units": "meters_per_second",
        "long_name": "gridded VVP fit root-mean-square residual",
    }
    dataset["valid_count"].attrs = {
        "units": "count",
        "long_name": "effective valid sample count contributing to the target voxel",
    }
    dataset["source_sweep_count"].attrs = {
        "units": "count",
        "long_name": "number of sweeps contributing to the target voxel",
    }
    dataset["neighbor_count"].attrs = {
        "units": "count",
        "long_name": "effective horizontal neighbors contributing to the target voxel",
    }
    dataset["effective_vertical_support"].attrs = {
        "units": "count",
        "long_name": "number of vertically supporting sweep samples used at the target voxel",
    }
    dataset["sweep_spread_m"].attrs = {
        "units": "meters",
        "long_name": "vertical spread of sweep samples contributing to the target voxel",
    }
    dataset["quality_score"].attrs = {
        "units": "1",
        "long_name": "heuristic voxel quality score from 0 to 100",
    }
    dataset["quality_flag"].attrs = {
        "units": "1",
        "long_name": "voxel quality flag",
        "flag_values": np.array([0, 1, 2, 3], dtype=np.uint8),
        "flag_meanings": "invalid low_confidence usable high_confidence",
    }
    return dataset


def _build_wind_volume_attrs(
    range_mode,
    requested_sweeps,
    selected_sweeps,
    rejected_sweeps,
    rejected_sweep_reasons,
    selection_mode,
    velocity_field_by_sweep,
    az_num,
    bin_num,
    azimuth_step,
    range_step,
    max_range_km,
    horizontal_radius_m,
    max_horizontal_radius_m,
    horizontal_min_neighbors,
    vertical_tolerance_m,
    max_vertical_gap_m,
    workers_requested,
    workers_used,
    parallel_mode,
):
    """Build shared dataset attrs for gridded wind-volume retrievals."""
    return {
        "method": "single_radar_horizontal_wind_volume",
        "source_method": "VVP",
        "range_mode": range_mode,
        "selection_mode": selection_mode,
        "requested_sweeps": json.dumps([int(sweep) for sweep in requested_sweeps]),
        "selected_sweeps": json.dumps([int(sweep) for sweep in selected_sweeps]),
        "rejected_sweeps": json.dumps([int(sweep) for sweep in rejected_sweeps]),
        "rejected_sweep_reasons": json.dumps({str(k): v for k, v in rejected_sweep_reasons.items()}, sort_keys=True),
        "velocity_field_used_by_sweep": json.dumps(
            {str(sweep): field for sweep, field in velocity_field_by_sweep.items()},
            sort_keys=True,
        ),
        "window_azimuth": int(az_num),
        "window_range": int(bin_num),
        "azimuth_step": int(azimuth_step),
        "range_step": int(range_step),
        "max_range_km": None if max_range_km is None else float(max_range_km),
        "horizontal_radius_m": float(horizontal_radius_m),
        "max_horizontal_radius_m": float(max_horizontal_radius_m),
        "horizontal_min_neighbors": int(horizontal_min_neighbors),
        "vertical_tolerance_m": float(vertical_tolerance_m),
        "max_vertical_gap_m": float(max_vertical_gap_m),
        "horizontal_interpolation": "idw",
        "vertical_interpolation": "linear_with_local_tolerance_aggregation",
        "workers_requested": int(workers_requested),
        "workers_used": int(workers_used),
        "parallel_mode": parallel_mode,
        "assumptions": (
            "Horizontal wind is retrieved locally with single-sweep VVP and "
            "then reconstructed onto fixed height levels without vertical-wind "
            "retrieval or extrapolation beyond the sampled height envelope."
        ),
        "references": json.dumps(WIND_REFERENCE_NOTES["WIND_VOLUME"]),
    }


def _extract_subset_dict(full_results, azimuth, elevation, range_values, x, y, z, azimuth_step, range_step):
    """Subset full VVP arrays and coordinates according to output steps."""
    azimuth_index = np.arange(0, azimuth.size, int(azimuth_step), dtype=np.int32)
    range_index = np.arange(0, range_values.size, int(range_step), dtype=np.int32)
    subset = {
        "azimuth": azimuth[azimuth_index],
        "elevation": elevation[azimuth_index],
        "range": range_values[range_index],
        "x": x[np.ix_(azimuth_index, range_index)],
        "y": y[np.ix_(azimuth_index, range_index)],
        "z": z[np.ix_(azimuth_index, range_index)],
    }
    for name, values in full_results.items():
        subset[name] = values[np.ix_(azimuth_index, range_index)]
    return subset


def _run_vvp_worker(task):
    """Run one VVP retrieval task for wind-volume generation."""
    full_results = _run_local_vvp(
        task["azimuth"],
        task["elevation"],
        task["values"],
        az_num=task["az_num"],
        bin_num=task["bin_num"],
        fillvalue=task["fillvalue"],
        min_valid_fraction=task["min_valid_fraction"],
        min_valid_count=task["min_valid_count"],
        min_azimuth_coverage_deg=task["min_azimuth_coverage_deg"],
        min_azimuth_sector_count=task["min_azimuth_sector_count"],
        max_condition_number=task["max_condition_number"],
    )
    subset = _extract_subset_dict(
        full_results,
        task["azimuth"],
        task["elevation"],
        task["range"],
        task["x"],
        task["y"],
        task["z"],
        task["azimuth_step"],
        task["range_step"],
    )
    subset["sweep"] = int(task["sweep"])
    subset["selected_field"] = task["selected_field"]
    subset["range_mode"] = task["range_mode"]
    return subset


def _summarize_vvp_result(retrieval, sweep):
    """Summarize one VVP retrieval for auto sweep selection."""
    u = _extract_result_array(retrieval, "u")
    fit_rmse = _extract_result_array(retrieval, "fit_rmse")
    valid_count = _extract_result_array(retrieval, "valid_count")
    valid_fraction = _extract_result_array(retrieval, "valid_fraction")
    condition_number = _extract_result_array(retrieval, "condition_number")
    azimuth_sector_count = _extract_result_array(retrieval, "azimuth_sector_count")
    x = _extract_result_array(retrieval, "x")
    y = _extract_result_array(retrieval, "y")
    elevation = _extract_result_array(retrieval, "elevation")
    finite = np.isfinite(u)
    finite_count = int(np.sum(finite))
    finite_fraction = float(finite_count / u.size) if u.size else 0.0
    sampled_range = np.hypot(x[finite], y[finite]) if finite_count > 0 else np.array([], dtype=np.float64)
    return {
        "sweep": int(sweep),
        "field_name": retrieval["selected_field"] if isinstance(retrieval, dict) else retrieval.attrs.get("velocity_field_used", ""),
        "finite_count": finite_count,
        "finite_fraction": finite_fraction,
        "median_fit_rmse": float(np.nanmedian(fit_rmse[finite])) if finite_count > 0 else np.inf,
        "median_valid_count": float(np.nanmedian(valid_count[finite])) if finite_count > 0 else 0.0,
        "median_valid_fraction": float(np.nanmedian(valid_fraction[finite])) if finite_count > 0 else 0.0,
        "median_condition_number": float(np.nanmedian(condition_number[finite])) if finite_count > 0 else np.inf,
        "median_sector_count": float(np.nanmedian(azimuth_sector_count[finite])) if finite_count > 0 else 0.0,
        "sampled_range_m": float(np.nanmax(sampled_range)) if sampled_range.size else 0.0,
        "fixed_angle_deg": float(np.nanmean(elevation)) if elevation.size else np.nan,
    }


def _score_sweep_summary(summary, reference_range_m):
    """Return a heuristic auto-selection score for one sweep summary."""
    if summary["finite_count"] <= 0:
        return 0.0
    range_component = min(1.0, float(summary["sampled_range_m"]) / max(float(reference_range_m), 1.0))
    fraction_component = min(1.0, float(summary["finite_fraction"]) / 0.15)
    rmse_component = float(np.exp(-max(float(summary["median_fit_rmse"]), 0.0) / 8.0))
    sector_component = min(1.0, float(summary["median_sector_count"]) / 6.0)
    elevation_component = 1.0 / (1.0 + max(float(summary["fixed_angle_deg"]), 0.0) / 8.0)
    return 100.0 * (
        0.25 * range_component
        + 0.25 * fraction_component
        + 0.20 * rmse_component
        + 0.15 * sector_component
        + 0.15 * elevation_component
    )


def _auto_select_sweeps(
    retrievals_by_sweep,
    *,
    max_selected_sweeps,
    sweep_min_valid_fraction,
    sweep_min_valid_count,
    sweep_min_azimuth_sector_count,
    sweep_max_median_fit_rmse,
    sweep_max_condition_number,
):
    """Select the best sweeps for wind-volume reconstruction."""
    summaries = {int(sweep): _summarize_vvp_result(retrieval, sweep) for sweep, retrieval in retrievals_by_sweep.items()}
    reference_range_m = max((summary["sampled_range_m"] for summary in summaries.values()), default=1.0)
    rejected_reasons = {}
    candidates = []
    for sweep, summary in summaries.items():
        reasons = []
        if summary["finite_count"] < int(sweep_min_valid_count):
            reasons.append("finite_count")
        if summary["finite_fraction"] < float(sweep_min_valid_fraction):
            reasons.append("finite_fraction")
        if summary["median_sector_count"] < float(sweep_min_azimuth_sector_count):
            reasons.append("azimuth_sector_count")
        if summary["median_fit_rmse"] > float(sweep_max_median_fit_rmse):
            reasons.append("fit_rmse")
        if summary["median_condition_number"] > float(sweep_max_condition_number):
            reasons.append("condition_number")
        summary["selection_score"] = _score_sweep_summary(summary, reference_range_m=reference_range_m)
        if reasons:
            rejected_reasons[int(sweep)] = reasons
            continue
        candidates.append(summary)
    candidates.sort(key=lambda item: (-item["selection_score"], item["fixed_angle_deg"], -item["sampled_range_m"]))
    selected = [int(item["sweep"]) for item in candidates[: int(max_selected_sweeps)]]
    if not selected:
        fallback = [
            item
            for item in summaries.values()
            if item["finite_count"] > 0
        ]
        fallback.sort(key=lambda item: (-item["selection_score"], item["fixed_angle_deg"], -item["sampled_range_m"]))
        if fallback:
            selected = [int(fallback[0]["sweep"])]
            rejected_reasons.pop(selected[0], None)
    rejected = [int(sweep) for sweep in summaries if int(sweep) not in selected]
    for sweep in rejected:
        rejected_reasons.setdefault(int(sweep), ["deprioritized"])
    return {
        "selected": selected,
        "rejected": rejected,
        "rejected_reasons": rejected_reasons,
        "summaries": summaries,
    }


def _reconstruct_wind_volume_chunk(
    start_index,
    stop_index,
    sweep_z,
    sweep_u,
    sweep_v,
    sweep_fit_rmse,
    sweep_valid_count,
    sweep_neighbor_count,
    sweep_expanded_radius,
    level_heights,
    vertical_tolerance_m,
    max_vertical_gap_m,
):
    """Reconstruct one x-slab of the target wind volume."""
    nz = int(level_heights.size)
    chunk_nx = int(stop_index - start_index)
    ny = int(sweep_z.shape[2])
    output = {
        "u": np.full((nz, chunk_nx, ny), np.nan, dtype=np.float64),
        "v": np.full((nz, chunk_nx, ny), np.nan, dtype=np.float64),
        "fit_rmse": np.full((nz, chunk_nx, ny), np.nan, dtype=np.float64),
        "valid_count": np.zeros((nz, chunk_nx, ny), dtype=np.int32),
        "source_sweep_count": np.zeros((nz, chunk_nx, ny), dtype=np.int32),
        "neighbor_count": np.zeros((nz, chunk_nx, ny), dtype=np.int32),
        "effective_vertical_support": np.zeros((nz, chunk_nx, ny), dtype=np.int32),
        "sweep_spread_m": np.full((nz, chunk_nx, ny), np.nan, dtype=np.float64),
        "quality_score": np.zeros((nz, chunk_nx, ny), dtype=np.float64),
        "quality_flag": np.zeros((nz, chunk_nx, ny), dtype=np.uint8),
    }
    for local_ix, ix in enumerate(range(int(start_index), int(stop_index))):
        for iy in range(ny):
            profile = _vertical_profile_from_sweeps(
                sweep_z[:, ix, iy],
                sweep_u[:, ix, iy],
                sweep_v[:, ix, iy],
                sweep_fit_rmse[:, ix, iy],
                sweep_valid_count[:, ix, iy],
                sweep_neighbor_count[:, ix, iy],
                sweep_expanded_radius[:, ix, iy],
                level_heights=level_heights,
                vertical_tolerance_m=float(vertical_tolerance_m),
                max_vertical_gap_m=float(max_vertical_gap_m),
            )
            for name in output:
                output[name][:, local_ix, iy] = profile[name]
    return int(start_index), int(stop_index), output


def _retrieve_wind_volume_to_target_grid(
    prd,
    target_x,
    target_y,
    level_heights,
    sweeps=None,
    field_name=None,
    range_mode="aligned",
    az_num=91,
    bin_num=9,
    azimuth_step=3,
    range_step=3,
    max_range_km=None,
    fillvalue=-999.0,
    min_valid_fraction=0.7,
    min_valid_count=20,
    min_azimuth_coverage_deg=0.0,
    min_azimuth_sector_count=2,
    max_condition_number=1.0e9,
    horizontal_radius_m=None,
    max_horizontal_radius_m=None,
    horizontal_min_neighbors=3,
    vertical_tolerance_m=None,
    max_vertical_gap_m=None,
    workers=1,
    max_selected_sweeps=6,
    sweep_min_valid_fraction=0.01,
    sweep_min_valid_count=8,
    sweep_min_azimuth_sector_count=2,
    sweep_max_median_fit_rmse=10.0,
    sweep_max_condition_number=None,
):
    """Compute the wind-volume arrays on an arbitrary 2-D target grid."""
    target_x = np.asarray(target_x, dtype=np.float64)
    target_y = np.asarray(target_y, dtype=np.float64)
    if target_x.shape != target_y.shape:
        raise ValueError("target_x and target_y must have the same shape")
    if target_x.ndim != 2:
        raise ValueError("target_x and target_y must be two-dimensional target grids")
    level_heights = _normalize_level_heights(level_heights)
    if horizontal_radius_m is None:
        horizontal_radius_m = _default_horizontal_radius_m(target_x[:, 0], target_y[0, :])
    if float(horizontal_radius_m) <= 0.0:
        raise ValueError("horizontal_radius_m must be positive")
    if max_horizontal_radius_m is None:
        max_horizontal_radius_m = max(float(horizontal_radius_m), 2.0 * float(horizontal_radius_m))
    if float(max_horizontal_radius_m) < float(horizontal_radius_m):
        raise ValueError("max_horizontal_radius_m must be greater than or equal to horizontal_radius_m")
    if int(horizontal_min_neighbors) <= 0:
        raise ValueError("horizontal_min_neighbors must be positive")
    if vertical_tolerance_m is None:
        vertical_tolerance_m = _default_vertical_tolerance_m(level_heights)
    if float(vertical_tolerance_m) <= 0.0:
        raise ValueError("vertical_tolerance_m must be positive")
    if max_vertical_gap_m is None:
        max_vertical_gap_m = _default_max_vertical_gap_m(level_heights)
    if float(max_vertical_gap_m) <= 0.0:
        raise ValueError("max_vertical_gap_m must be positive")
    selection_mode, requested_sweeps = _resolve_requested_sweeps(prd, sweeps)
    if not requested_sweeps:
        raise ValueError("at least one sweep is required for wind-volume retrieval")
    if sweep_max_condition_number is None:
        sweep_max_condition_number = max_condition_number
    velocity_field_by_sweep = {}
    tasks = []
    for sweep in requested_sweeps:
        selected_field = _resolve_velocity_field_name(prd, sweep=sweep, field_name=field_name)
        velocity_field_by_sweep[int(sweep)] = selected_field
        field = prd.get_sweep_field(int(sweep), selected_field, range_mode=range_mode, sort_by_azimuth=True)
        field = _apply_range_limit(field, max_range_km=max_range_km)
        tasks.append(
            {
                "sweep": int(sweep),
                "selected_field": selected_field,
                "range_mode": range_mode,
                "azimuth": np.asarray(field.azimuth.values, dtype=np.float64),
                "elevation": np.asarray(field.elevation.values, dtype=np.float64),
                "range": np.asarray(field.range.values, dtype=np.float64),
                "x": np.asarray(field["x"].values, dtype=np.float64),
                "y": np.asarray(field["y"].values, dtype=np.float64),
                "z": np.asarray(field["z"].values, dtype=np.float64),
                "values": np.asarray(field.values, dtype=np.float64),
                "az_num": int(az_num),
                "bin_num": int(bin_num),
                "fillvalue": float(fillvalue),
                "min_valid_fraction": float(min_valid_fraction),
                "min_valid_count": int(min_valid_count),
                "min_azimuth_coverage_deg": float(min_azimuth_coverage_deg),
                "min_azimuth_sector_count": int(min_azimuth_sector_count),
                "max_condition_number": float(max_condition_number),
                "azimuth_step": int(azimuth_step),
                "range_step": int(range_step),
            }
        )
    sweep_workers = _resolve_worker_count(workers, len(tasks))
    if sweep_workers > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=sweep_workers) as executor:
            retrieval_results = list(executor.map(_run_vvp_worker, tasks))
        parallel_mode = "sweeps"
    else:
        retrieval_results = [_run_vvp_worker(task) for task in tasks]
        parallel_mode = "serial"
    retrievals_by_sweep = {int(result["sweep"]): result for result in retrieval_results}
    if selection_mode == "explicit":
        selected_sweeps = [int(sweep) for sweep in requested_sweeps]
        rejected_sweeps = []
        rejected_reasons = {}
    elif selection_mode == "all":
        selected_sweeps = [int(sweep) for sweep in requested_sweeps]
        rejected_sweeps = []
        rejected_reasons = {}
    else:
        selection = _auto_select_sweeps(
            retrievals_by_sweep,
            max_selected_sweeps=int(max_selected_sweeps),
            sweep_min_valid_fraction=float(sweep_min_valid_fraction),
            sweep_min_valid_count=int(sweep_min_valid_count),
            sweep_min_azimuth_sector_count=int(sweep_min_azimuth_sector_count),
            sweep_max_median_fit_rmse=float(sweep_max_median_fit_rmse),
            sweep_max_condition_number=float(sweep_max_condition_number),
        )
        selected_sweeps = selection["selected"]
        rejected_sweeps = selection["rejected"]
        rejected_reasons = selection["rejected_reasons"]
    if not selected_sweeps:
        raise ValueError("wind-volume retrieval could not find any usable sweeps")
    gridded_sweeps = []
    for sweep in selected_sweeps:
        retrieval = retrievals_by_sweep[int(sweep)]
        gridded_sweeps.append(
            _grid_vvp_retrieval_xy(
                retrieval,
                target_x=target_x,
                target_y=target_y,
                horizontal_radius_m=float(horizontal_radius_m),
                max_horizontal_radius_m=float(max_horizontal_radius_m),
                horizontal_min_neighbors=int(horizontal_min_neighbors),
            )
        )

    nsweeps = len(gridded_sweeps)
    nx, ny = target_x.shape
    nz = level_heights.size
    sweep_u = np.stack([item["u"] for item in gridded_sweeps], axis=0)
    sweep_v = np.stack([item["v"] for item in gridded_sweeps], axis=0)
    sweep_z = np.stack([item["z"] for item in gridded_sweeps], axis=0)
    sweep_fit_rmse = np.stack([item["fit_rmse"] for item in gridded_sweeps], axis=0)
    sweep_valid_count = np.stack([item["valid_count"] for item in gridded_sweeps], axis=0).astype(np.float64)
    sweep_neighbor_count = np.stack([item["neighbor_count"] for item in gridded_sweeps], axis=0).astype(np.float64)
    sweep_expanded_radius = np.stack([item["expanded_radius"] for item in gridded_sweeps], axis=0).astype(np.float64)

    u_volume = np.full((nz, nx, ny), np.nan, dtype=np.float64)
    v_volume = np.full((nz, nx, ny), np.nan, dtype=np.float64)
    fit_rmse_volume = np.full((nz, nx, ny), np.nan, dtype=np.float64)
    valid_count_volume = np.zeros((nz, nx, ny), dtype=np.int32)
    source_sweep_count_volume = np.zeros((nz, nx, ny), dtype=np.int32)
    neighbor_count_volume = np.zeros((nz, nx, ny), dtype=np.int32)
    effective_vertical_support_volume = np.zeros((nz, nx, ny), dtype=np.int32)
    sweep_spread_m_volume = np.full((nz, nx, ny), np.nan, dtype=np.float64)
    quality_score_volume = np.zeros((nz, nx, ny), dtype=np.float64)
    quality_flag_volume = np.zeros((nz, nx, ny), dtype=np.uint8)
    grid_task_count = nx if (nx * ny) >= 64 else 1
    grid_workers = _resolve_worker_count(workers, grid_task_count)
    if grid_workers > 1:
        chunk_edges = np.linspace(0, nx, grid_workers + 1, dtype=np.int32)
        with concurrent.futures.ProcessPoolExecutor(max_workers=grid_workers) as executor:
            futures = [
                executor.submit(
                    _reconstruct_wind_volume_chunk,
                    int(chunk_edges[index]),
                    int(chunk_edges[index + 1]),
                    sweep_z,
                    sweep_u,
                    sweep_v,
                    sweep_fit_rmse,
                    sweep_valid_count,
                    sweep_neighbor_count,
                    sweep_expanded_radius,
                    level_heights,
                    float(vertical_tolerance_m),
                    float(max_vertical_gap_m),
                )
                for index in range(grid_workers)
                if int(chunk_edges[index + 1]) > int(chunk_edges[index])
            ]
            for future in concurrent.futures.as_completed(futures):
                start_index, stop_index, chunk = future.result()
                target_slice = slice(start_index, stop_index)
                u_volume[:, target_slice, :] = chunk["u"]
                v_volume[:, target_slice, :] = chunk["v"]
                fit_rmse_volume[:, target_slice, :] = chunk["fit_rmse"]
                valid_count_volume[:, target_slice, :] = chunk["valid_count"]
                source_sweep_count_volume[:, target_slice, :] = chunk["source_sweep_count"]
                neighbor_count_volume[:, target_slice, :] = chunk["neighbor_count"]
                effective_vertical_support_volume[:, target_slice, :] = chunk["effective_vertical_support"]
                sweep_spread_m_volume[:, target_slice, :] = chunk["sweep_spread_m"]
                quality_score_volume[:, target_slice, :] = chunk["quality_score"]
                quality_flag_volume[:, target_slice, :] = chunk["quality_flag"]
        parallel_mode = "sweeps+grid" if parallel_mode == "sweeps" else "grid"
    else:
        _, _, chunk = _reconstruct_wind_volume_chunk(
            0,
            nx,
            sweep_z,
            sweep_u,
            sweep_v,
            sweep_fit_rmse,
            sweep_valid_count,
            sweep_neighbor_count,
            sweep_expanded_radius,
            level_heights,
            float(vertical_tolerance_m),
            float(max_vertical_gap_m),
        )
        u_volume[:, :, :] = chunk["u"]
        v_volume[:, :, :] = chunk["v"]
        fit_rmse_volume[:, :, :] = chunk["fit_rmse"]
        valid_count_volume[:, :, :] = chunk["valid_count"]
        source_sweep_count_volume[:, :, :] = chunk["source_sweep_count"]
        neighbor_count_volume[:, :, :] = chunk["neighbor_count"]
        effective_vertical_support_volume[:, :, :] = chunk["effective_vertical_support"]
        sweep_spread_m_volume[:, :, :] = chunk["sweep_spread_m"]
        quality_score_volume[:, :, :] = chunk["quality_score"]
        quality_flag_volume[:, :, :] = chunk["quality_flag"]

    attrs = _build_wind_volume_attrs(
        range_mode=range_mode,
        requested_sweeps=requested_sweeps,
        selected_sweeps=selected_sweeps,
        rejected_sweeps=rejected_sweeps,
        rejected_sweep_reasons=rejected_reasons,
        selection_mode=selection_mode,
        velocity_field_by_sweep=velocity_field_by_sweep,
        az_num=az_num,
        bin_num=bin_num,
        azimuth_step=azimuth_step,
        range_step=range_step,
        max_range_km=max_range_km,
        horizontal_radius_m=float(horizontal_radius_m),
        max_horizontal_radius_m=float(max_horizontal_radius_m),
        horizontal_min_neighbors=int(horizontal_min_neighbors),
        vertical_tolerance_m=float(vertical_tolerance_m),
        max_vertical_gap_m=float(max_vertical_gap_m),
        workers_requested=int(workers),
        workers_used=max(sweep_workers, grid_workers),
        parallel_mode=parallel_mode,
    )
    attrs["sweep_count"] = int(nsweeps)
    return {
        "u": u_volume,
        "v": v_volume,
        "fit_rmse": fit_rmse_volume,
        "valid_count": valid_count_volume,
        "source_sweep_count": source_sweep_count_volume,
        "neighbor_count": neighbor_count_volume,
        "effective_vertical_support": effective_vertical_support_volume,
        "sweep_spread_m": sweep_spread_m_volume,
        "quality_score": quality_score_volume,
        "quality_flag": quality_flag_volume,
        "attrs": attrs,
    }


def retrieve_vad(
    prd,
    sweeps=None,
    field_name=None,
    range_mode="aligned",
    gate_step=1,
    max_range_km=None,
    fillvalue=-999.0,
    min_valid_fraction=0.5,
    min_valid_count=16,
    min_azimuth_coverage_deg=240.0,
    max_condition_number=1.0e6,
):
    """
    Retrieve horizontal wind ring-by-ring from one or more PPI sweeps.

    Parameters
    ----------
    prd : pycwr.core.NRadar.PRD
        Radar volume.
    sweeps : int or sequence of int, optional
        Sweep indices to retrieve. ``None`` means all sweeps.
    field_name : str, optional
        Velocity field to use. Defaults to ``Vc`` when present, else ``V``.
    range_mode : str
        Range mode for the velocity field. Velocity retrieval defaults to the
        aligned grid because Doppler moments share the Doppler range.
    gate_step : int
        Subsampling step along range gates.

    Returns
    -------
    xarray.Dataset
        Retrieval dataset with padded ``(sweep, gate)`` dimensions.
    """
    if gate_step <= 0:
        raise ValueError("gate_step must be a positive integer")
    if sweeps is None:
        sweep_indices = [int(sweep) for sweep in np.asarray(prd.scan_info.sweep.values, dtype=np.int32)]
    elif np.isscalar(sweeps):
        sweep_indices = [int(sweeps)]
    else:
        sweep_indices = [int(sweep) for sweep in sweeps]
    if not sweep_indices:
        raise ValueError("at least one sweep is required for VAD retrieval")

    fields = []
    velocity_field_by_sweep = {}
    max_gate_count = 0
    for sweep in sweep_indices:
        selected_field = _resolve_velocity_field_name(prd, sweep=sweep, field_name=field_name)
        velocity_field_by_sweep[sweep] = selected_field
        field = prd.get_sweep_field(sweep, selected_field, range_mode=range_mode, sort_by_azimuth=True)
        field = _apply_range_limit(field, max_range_km=max_range_km)
        field = field.isel(range=slice(0, None, int(gate_step)))
        fields.append((sweep, selected_field, field))
        max_gate_count = max(max_gate_count, int(field.sizes["range"]))

    output = _init_vad_output((len(fields), max_gate_count))
    fixed_angles = np.full(len(fields), np.nan, dtype=np.float64)
    rays_per_sweep = np.zeros(len(fields), dtype=np.int32)

    for row_index, (sweep, selected_field, field) in enumerate(fields):
        azimuth = np.asarray(field.azimuth.values, dtype=np.float64)
        elevation = np.asarray(field.elevation.values, dtype=np.float64)
        values = _as_float_array(field.values, fillvalue=fillvalue)
        ranges = np.asarray(field.range.values, dtype=np.float64)
        heights = _gate_heights_from_field(field)
        fixed_angles[row_index] = float(np.nanmean(elevation))
        rays_per_sweep[row_index] = int(field.sizes["azimuth"])
        gate_count = int(field.sizes["range"])
        output["range"][row_index, :gate_count] = ranges
        output["height"][row_index, :gate_count] = heights
        for gate_index in range(gate_count):
            result = fit_vad_ring(
                azimuth,
                elevation,
                values[:, gate_index],
                fillvalue=np.nan,
                min_valid_fraction=min_valid_fraction,
                min_valid_count=min_valid_count,
                min_azimuth_coverage_deg=min_azimuth_coverage_deg,
                max_condition_number=max_condition_number,
            )
            for name in (
                "u",
                "v",
                "wind_speed",
                "wind_direction",
                "radial_offset",
                "fit_rmse",
                "condition_number",
                "valid_fraction",
                "azimuth_coverage_deg",
            ):
                output[name][row_index, gate_index] = result[name]
            output["valid_count"][row_index, gate_index] = int(result["valid_count"])

    dataset = xr.Dataset(
        data_vars={
            "u": (("sweep", "gate"), output["u"]),
            "v": (("sweep", "gate"), output["v"]),
            "wind_speed": (("sweep", "gate"), output["wind_speed"]),
            "wind_direction": (("sweep", "gate"), output["wind_direction"]),
            "radial_offset": (("sweep", "gate"), output["radial_offset"]),
            "fit_rmse": (("sweep", "gate"), output["fit_rmse"]),
            "condition_number": (("sweep", "gate"), output["condition_number"]),
            "valid_count": (("sweep", "gate"), output["valid_count"]),
            "valid_fraction": (("sweep", "gate"), output["valid_fraction"]),
            "azimuth_coverage_deg": (("sweep", "gate"), output["azimuth_coverage_deg"]),
            "range": (("sweep", "gate"), output["range"]),
            "height": (("sweep", "gate"), output["height"]),
            "fixed_angle": (("sweep",), fixed_angles),
            "rays_per_sweep": (("sweep",), rays_per_sweep),
        },
        coords={
            "sweep": np.asarray(sweep_indices, dtype=np.int32),
            "gate": np.arange(max_gate_count, dtype=np.int32),
        },
        attrs={
            "method": "VAD",
            "algorithm": "first_harmonic_weighted_least_squares",
            "range_mode": range_mode,
            "fillvalue": float(fillvalue),
            "gate_step": int(gate_step),
            "max_range_km": None if max_range_km is None else float(max_range_km),
            "velocity_field_used_by_sweep": json.dumps(
                {str(sweep): name for sweep, name in velocity_field_by_sweep.items()},
                sort_keys=True,
            ),
            "assumptions": (
                "Uniform horizontal wind within each sampled VAD ring; "
                "intercept term is retained as a radial offset rather than "
                "interpreted as vertical air motion."
            ),
            "references": json.dumps(WIND_REFERENCE_NOTES["VAD"]),
        },
    )
    dataset["u"].attrs = {
        "units": "meters_per_second",
        "standard_name": "eastward_wind_component",
        "long_name": "retrieved eastward wind component from single-radar VAD",
    }
    dataset["v"].attrs = {
        "units": "meters_per_second",
        "standard_name": "northward_wind_component",
        "long_name": "retrieved northward wind component from single-radar VAD",
    }
    dataset["wind_speed"].attrs = {
        "units": "meters_per_second",
        "long_name": "horizontal wind speed",
    }
    dataset["wind_direction"].attrs = {
        "units": "degrees",
        "long_name": "meteorological wind direction",
    }
    dataset["radial_offset"].attrs = {
        "units": "meters_per_second",
        "long_name": "fitted radial offset term",
    }
    dataset["fit_rmse"].attrs = {
        "units": "meters_per_second",
        "long_name": "root-mean-square fit residual",
    }
    dataset["condition_number"].attrs = {
        "units": "unitless",
        "long_name": "normal-equation condition number",
    }
    dataset["valid_count"].attrs = {"long_name": "number of valid radial-velocity samples in the fitted ring"}
    dataset["valid_fraction"].attrs = {"units": "1", "long_name": "fraction of valid samples in the fitted ring"}
    dataset["azimuth_coverage_deg"].attrs = {"units": "degrees", "long_name": "azimuthal coverage of valid samples"}
    dataset["range"].attrs = dict(DEFAULT_METADATA["range"])
    dataset["height"].attrs = {
        "units": "meters",
        "standard_name": "height",
        "long_name": "mean gate height above radar mean sea level",
    }
    dataset["fixed_angle"].attrs = dict(DEFAULT_METADATA["fixed_angle"])
    return dataset


def retrieve_vwp(
    prd,
    sweeps=None,
    field_name=None,
    range_mode="aligned",
    height_bins=None,
    height_step=250.0,
    max_height_m=None,
    gate_step=1,
    max_range_km=None,
    fillvalue=-999.0,
    min_valid_fraction=0.4,
    min_valid_count=12,
    min_azimuth_coverage_deg=180.0,
    max_condition_number=1.0e6,
    vertical_smoothing_bins=3,
    interpolate_gaps=True,
    max_gap_bins=2,
):
    """
    Build a robust vertical wind profile product from VAD retrievals.
    """
    generated_height_bins = height_bins is None
    vad = retrieve_vad(
        prd,
        sweeps=sweeps,
        field_name=field_name,
        range_mode=range_mode,
        gate_step=gate_step,
        max_range_km=max_range_km,
        fillvalue=fillvalue,
        min_valid_fraction=min_valid_fraction,
        min_valid_count=min_valid_count,
        min_azimuth_coverage_deg=min_azimuth_coverage_deg,
        max_condition_number=max_condition_number,
    )
    height = np.asarray(vad["height"].values, dtype=np.float64)
    if height_bins is None:
        finite_height = height[np.isfinite(height)]
        if finite_height.size == 0:
            raise ValueError("VWP requires at least one finite gate height")
        step = float(height_step)
        if step <= 0.0:
            raise ValueError("height_step must be positive")
        height_max = float(np.nanmax(finite_height)) if max_height_m is None else float(max_height_m)
        height_bins = np.arange(0.0, height_max + step, step, dtype=np.float64)
        if height_bins.size < 2:
            height_bins = np.array([0.0, step], dtype=np.float64)
    else:
        height_bins = np.asarray(height_bins, dtype=np.float64).reshape(-1)
        if height_bins.size < 2:
            raise ValueError("height_bins must contain at least two edges")
    height_center = 0.5 * (height_bins[:-1] + height_bins[1:])

    u_values = np.asarray(vad["u"].values, dtype=np.float64).reshape(-1)
    v_values = np.asarray(vad["v"].values, dtype=np.float64).reshape(-1)
    rmse = np.asarray(vad["fit_rmse"].values, dtype=np.float64).reshape(-1)
    valid_fraction = np.asarray(vad["valid_fraction"].values, dtype=np.float64).reshape(-1)
    coverage = np.asarray(vad["azimuth_coverage_deg"].values, dtype=np.float64).reshape(-1)
    height_values = height.reshape(-1)

    base_weight = valid_fraction * np.clip(coverage / 360.0, 0.0, 1.0) / np.maximum(rmse, 0.5)
    base_weight = np.where(np.isfinite(base_weight), base_weight, 0.0)

    profile_u = np.full(height_center.shape, np.nan, dtype=np.float64)
    profile_v = np.full(height_center.shape, np.nan, dtype=np.float64)
    profile_count = np.zeros(height_center.shape, dtype=np.int32)
    profile_weight = np.zeros(height_center.shape, dtype=np.float64)
    profile_rmse = np.full(height_center.shape, np.nan, dtype=np.float64)

    for index in range(height_center.size):
        lower = height_bins[index]
        upper = height_bins[index + 1]
        mask = (
            np.isfinite(height_values)
            & (height_values >= lower)
            & (height_values < upper)
            & np.isfinite(u_values)
            & np.isfinite(v_values)
            & (base_weight > 0.0)
        )
        profile_count[index] = int(mask.sum())
        if profile_count[index] == 0:
            continue
        weights = base_weight[mask]
        profile_weight[index] = float(np.sum(weights))
        profile_u[index] = _robust_weighted_location(u_values[mask], weights)
        profile_v[index] = _robust_weighted_location(v_values[mask], weights)
        if np.any(np.isfinite(rmse[mask])):
            profile_rmse[index] = float(np.average(rmse[mask], weights=np.maximum(weights, 1.0e-12)))

    if vertical_smoothing_bins and int(vertical_smoothing_bins) > 1:
        profile_u = _vertical_smooth(profile_u, np.maximum(profile_weight, 1.0e-12), window=int(vertical_smoothing_bins))
        profile_v = _vertical_smooth(profile_v, np.maximum(profile_weight, 1.0e-12), window=int(vertical_smoothing_bins))
    if interpolate_gaps:
        profile_u = _interpolate_short_gaps(profile_u, max_gap_bins=max_gap_bins)
        profile_v = _interpolate_short_gaps(profile_v, max_gap_bins=max_gap_bins)

    wind_speed = np.hypot(profile_u, profile_v)
    wind_direction = _wind_direction_from_uv(profile_u, profile_v)
    dataset = xr.Dataset(
        data_vars={
            "u": (("height",), profile_u),
            "v": (("height",), profile_v),
            "wind_speed": (("height",), wind_speed),
            "wind_direction": (("height",), wind_direction),
            "sample_count": (("height",), profile_count),
            "sample_weight": (("height",), profile_weight),
            "fit_rmse": (("height",), profile_rmse),
            "height_lower": (("height",), height_bins[:-1]),
            "height_upper": (("height",), height_bins[1:]),
        },
        coords={"height": height_center},
        attrs={
            "method": "VWP",
            "source_method": "VAD",
            "range_mode": range_mode,
            "height_step": float(height_step) if generated_height_bins else np.nan,
            "gate_step": int(gate_step),
            "interpolate_gaps": "true" if interpolate_gaps else "false",
            "max_gap_bins": int(max_gap_bins),
            "vertical_smoothing_bins": int(vertical_smoothing_bins),
            "references": json.dumps(WIND_REFERENCE_NOTES["VAD"]),
            "velocity_field_used_by_sweep": vad.attrs["velocity_field_used_by_sweep"],
            "assumptions": (
                "Profile is aggregated from VAD-retrieved horizontal winds; "
                "short vertical gaps may be interpolated when adjacent bins are valid."
            ),
        },
    )
    dataset["height"].attrs = {
        "units": "meters",
        "standard_name": "height",
        "long_name": "wind-profile bin center height",
        "positive": "up",
    }
    dataset["u"].attrs = {
        "units": "meters_per_second",
        "standard_name": "eastward_wind_component",
        "long_name": "vertical wind profile eastward component",
    }
    dataset["v"].attrs = {
        "units": "meters_per_second",
        "standard_name": "northward_wind_component",
        "long_name": "vertical wind profile northward component",
    }
    dataset["wind_speed"].attrs = {"units": "meters_per_second", "long_name": "vertical wind profile speed"}
    dataset["wind_direction"].attrs = {"units": "degrees", "long_name": "vertical wind profile direction"}
    dataset["sample_count"].attrs = {"units": "count", "long_name": "number of VAD samples contributing to the bin"}
    dataset["sample_weight"].attrs = {"units": "unitless", "long_name": "total robust aggregation weight in the bin"}
    dataset["fit_rmse"].attrs = {"units": "meters_per_second", "long_name": "weighted mean VAD fit RMSE in the bin"}
    dataset["height_lower"].attrs = {"units": "meters", "long_name": "lower bound of the profile bin"}
    dataset["height_upper"].attrs = {"units": "meters", "long_name": "upper bound of the profile bin"}
    return dataset


def _vvp_window_weights(az_num, bin_num):
    """Return separable Gaussian-like window weights for local VVP retrieval."""
    az_offsets = np.arange(az_num, dtype=np.float64) - (az_num - 1) / 2.0
    bin_offsets = np.arange(bin_num, dtype=np.float64) - (bin_num - 1) / 2.0
    sigma_az = max(1.0, (az_num - 1) / 3.0)
    sigma_bin = max(1.0, (bin_num - 1) / 3.0)
    az_weight = np.exp(-0.5 * (az_offsets / sigma_az) ** 2)
    bin_weight = np.exp(-0.5 * (bin_offsets / sigma_bin) ** 2)
    return np.outer(az_weight, bin_weight)


def _run_local_vvp(
    azimuth,
    elevation,
    radial_velocity,
    az_num,
    bin_num,
    fillvalue=-999.0,
    min_valid_fraction=0.7,
    min_valid_count=20,
    min_azimuth_coverage_deg=0.0,
    min_azimuth_sector_count=2,
    max_condition_number=1.0e9,
):
    """Run the local moving-window VVP retrieval on one sweep array."""
    if az_num % 2 != 1:
        raise ValueError("az_num must be odd")
    if bin_num % 2 != 1:
        raise ValueError("bin_num must be odd")
    radial_velocity = _as_float_array(radial_velocity, fillvalue=fillvalue)
    azimuth = np.asarray(azimuth, dtype=np.float64).reshape(-1)
    if azimuth.size != radial_velocity.shape[0]:
        raise ValueError("azimuth length must match the first radial_velocity dimension")
    if np.asarray(elevation).ndim == 0:
        elevation = np.full(azimuth.shape, float(elevation), dtype=np.float64)
    else:
        elevation = np.asarray(elevation, dtype=np.float64).reshape(-1)
    if elevation.shape != azimuth.shape:
        raise ValueError("elevation must be scalar or match the azimuth shape")

    naz, nrange = radial_velocity.shape
    half_az = az_num // 2
    half_bin = bin_num // 2
    padded_velocity = np.pad(radial_velocity, ((half_az, half_az), (0, 0)), mode="wrap")
    padded_azimuth = np.pad(azimuth, (half_az, half_az), mode="wrap")
    padded_elevation = np.pad(elevation, (half_az, half_az), mode="wrap")
    window_weights = _vvp_window_weights(az_num, bin_num)

    outputs = {
        "u": np.full((naz, nrange), np.nan, dtype=np.float64),
        "v": np.full((naz, nrange), np.nan, dtype=np.float64),
        "wind_speed": np.full((naz, nrange), np.nan, dtype=np.float64),
        "wind_direction": np.full((naz, nrange), np.nan, dtype=np.float64),
        "radial_offset": np.full((naz, nrange), np.nan, dtype=np.float64),
        "fit_rmse": np.full((naz, nrange), np.nan, dtype=np.float64),
        "condition_number": np.full((naz, nrange), np.nan, dtype=np.float64),
        "valid_fraction": np.full((naz, nrange), np.nan, dtype=np.float64),
        "azimuth_coverage_deg": np.full((naz, nrange), np.nan, dtype=np.float64),
        "azimuth_sector_count": np.zeros((naz, nrange), dtype=np.int32),
        "valid_count": np.zeros((naz, nrange), dtype=np.int32),
    }

    for az_index in range(naz):
        az_slice = slice(az_index, az_index + az_num)
        az_window = padded_azimuth[az_slice]
        el_window = padded_elevation[az_slice]
        required_sector_count = _effective_min_sector_count(
            min_azimuth_sector_count,
            az_window,
        )
        for range_index in range(half_bin, nrange - half_bin):
            window = padded_velocity[az_slice, range_index - half_bin : range_index + half_bin + 1]
            valid = np.isfinite(window)
            valid_count = int(valid.sum())
            valid_fraction = float(valid_count / window.size)
            azimuth_coverage = _azimuth_coverage_deg(np.repeat(az_window, bin_num)[valid.reshape(-1)])
            azimuth_sector_count = _azimuth_sector_count(np.repeat(az_window, bin_num)[valid.reshape(-1)])
            outputs["valid_count"][az_index, range_index] = valid_count
            outputs["valid_fraction"][az_index, range_index] = valid_fraction
            outputs["azimuth_coverage_deg"][az_index, range_index] = azimuth_coverage
            outputs["azimuth_sector_count"][az_index, range_index] = azimuth_sector_count
            if (
                valid_count < int(min_valid_count)
                or valid_fraction < float(min_valid_fraction)
                or azimuth_coverage < float(min_azimuth_coverage_deg)
                or azimuth_sector_count < int(required_sector_count)
            ):
                continue

            design = _design_matrix(np.repeat(az_window, bin_num), np.repeat(el_window, bin_num))
            fit = _solve_weighted_least_squares(
                design,
                window.reshape(-1),
                weights=window_weights.reshape(-1),
                max_condition_number=max_condition_number,
            )
            if fit is None or fit["coefficients"] is None:
                if fit is not None:
                    outputs["condition_number"][az_index, range_index] = fit["condition_number"]
                continue
            u = float(fit["coefficients"][0])
            v = float(fit["coefficients"][1])
            outputs["u"][az_index, range_index] = u
            outputs["v"][az_index, range_index] = v
            outputs["wind_speed"][az_index, range_index] = float(np.hypot(u, v))
            outputs["wind_direction"][az_index, range_index] = float(_wind_direction_from_uv(u, v))
            outputs["radial_offset"][az_index, range_index] = float(fit["coefficients"][2])
            outputs["fit_rmse"][az_index, range_index] = float(fit["rmse"])
            outputs["condition_number"][az_index, range_index] = float(fit["condition_number"])
    return outputs


def retrieve_vvp(
    prd,
    sweep,
    field_name=None,
    range_mode="aligned",
    az_num=91,
    bin_num=9,
    azimuth_step=1,
    range_step=1,
    max_range_km=None,
    fillvalue=-999.0,
    min_valid_fraction=0.7,
    min_valid_count=20,
    min_azimuth_coverage_deg=0.0,
    min_azimuth_sector_count=2,
    max_condition_number=1.0e9,
):
    """
    Retrieve a local horizontal wind field from one PPI sweep.

    Returns
    -------
    xarray.Dataset
        Local VVP analysis on the requested sweep.

    Notes
    -----
    Single-sweep VVP becomes ill-conditioned when the azimuth window is too
    narrow. The default ``az_num=91`` is intentionally wider than a cosmetic
    smoothing window so that the horizontal-wind fit remains stable on real
    radar PPIs.
    """
    if azimuth_step <= 0 or range_step <= 0:
        raise ValueError("azimuth_step and range_step must be positive integers")
    sweep = int(sweep)
    selected_field = _resolve_velocity_field_name(prd, sweep=sweep, field_name=field_name)
    field = prd.get_sweep_field(sweep, selected_field, range_mode=range_mode, sort_by_azimuth=True)
    field = _apply_range_limit(field, max_range_km=max_range_km)
    full_results = _run_local_vvp(
        np.asarray(field.azimuth.values, dtype=np.float64),
        np.asarray(field.elevation.values, dtype=np.float64),
        np.asarray(field.values, dtype=np.float64),
        az_num=az_num,
        bin_num=bin_num,
        fillvalue=fillvalue,
        min_valid_fraction=min_valid_fraction,
        min_valid_count=min_valid_count,
        min_azimuth_coverage_deg=min_azimuth_coverage_deg,
        min_azimuth_sector_count=min_azimuth_sector_count,
        max_condition_number=max_condition_number,
    )
    azimuth_index = np.arange(0, field.sizes["azimuth"], int(azimuth_step), dtype=np.int32)
    range_index = np.arange(0, field.sizes["range"], int(range_step), dtype=np.int32)
    subset = field.isel(azimuth=azimuth_index, range=range_index)

    dataset = xr.Dataset(coords=subset.coords)
    for name, values in full_results.items():
        dataset[name] = (subset.dims, values[np.ix_(azimuth_index, range_index)])
    dataset.attrs.update(
        {
            "method": "VVP",
            "algorithm": "local_weighted_least_squares",
            "sweep": sweep,
            "range_mode": range_mode,
            "velocity_field_used": selected_field,
            "window_azimuth": int(az_num),
            "window_range": int(bin_num),
            "azimuth_step": int(azimuth_step),
            "range_step": int(range_step),
            "max_range_km": None if max_range_km is None else float(max_range_km),
            "min_azimuth_sector_count": int(min_azimuth_sector_count),
            "assumptions": (
                "Uniform horizontal wind within each local moving window; "
                "the fitted intercept is retained as a radial offset term."
            ),
            "references": json.dumps(WIND_REFERENCE_NOTES["VVP"]),
        }
    )
    dataset["u"].attrs = {
        "units": "meters_per_second",
        "standard_name": "eastward_wind_component",
        "long_name": "retrieved eastward wind component from local single-radar VVP",
    }
    dataset["v"].attrs = {
        "units": "meters_per_second",
        "standard_name": "northward_wind_component",
        "long_name": "retrieved northward wind component from local single-radar VVP",
    }
    dataset["wind_speed"].attrs = {
        "units": "meters_per_second",
        "long_name": "horizontal wind speed",
    }
    dataset["wind_direction"].attrs = {
        "units": "degrees",
        "long_name": "meteorological wind direction",
    }
    dataset["radial_offset"].attrs = {
        "units": "meters_per_second",
        "long_name": "fitted radial offset term",
    }
    dataset["fit_rmse"].attrs = {
        "units": "meters_per_second",
        "long_name": "root-mean-square fit residual",
    }
    dataset["condition_number"].attrs = {
        "units": "unitless",
        "long_name": "normal-equation condition number",
    }
    dataset["valid_fraction"].attrs = {"units": "1", "long_name": "fraction of valid samples in the local window"}
    dataset["azimuth_coverage_deg"].attrs = {"units": "degrees", "long_name": "azimuthal coverage of valid samples"}
    dataset["azimuth_sector_count"].attrs = {"units": "count", "long_name": "number of occupied azimuth sectors in the local window"}
    dataset["valid_count"].attrs = {"long_name": "number of valid samples in the local retrieval window"}
    return dataset


def retrieve_wind_volume_xy(
    prd,
    XRange,
    YRange,
    level_heights,
    sweeps=None,
    field_name=None,
    range_mode="aligned",
    az_num=91,
    bin_num=9,
    azimuth_step=3,
    range_step=3,
    max_range_km=None,
    fillvalue=-999.0,
    min_valid_fraction=0.7,
    min_valid_count=20,
    min_azimuth_coverage_deg=0.0,
    min_azimuth_sector_count=2,
    max_condition_number=1.0e9,
    horizontal_radius_m=None,
    max_horizontal_radius_m=None,
    horizontal_min_neighbors=3,
    vertical_tolerance_m=None,
    max_vertical_gap_m=None,
    workers=1,
    max_selected_sweeps=6,
    sweep_min_valid_fraction=0.01,
    sweep_min_valid_count=8,
    sweep_min_azimuth_sector_count=2,
    sweep_max_median_fit_rmse=10.0,
    sweep_max_condition_number=None,
):
    """
    Retrieve a gridded 3-D horizontal wind volume on a Cartesian grid.

    The returned dataset stores one wind profile at each horizontal target
    grid point. This is a single-radar ``u/v`` wind-volume product derived
    from local VVP retrievals and fixed-height reconstruction.
    """
    x_axis = _as_1d_float(XRange)
    y_axis = _as_1d_float(YRange)
    if x_axis.size == 0 or y_axis.size == 0:
        raise ValueError("XRange and YRange must contain at least one value")
    level_heights = _normalize_level_heights(level_heights)
    target_x, target_y = np.meshgrid(x_axis, y_axis, indexing="ij")
    outputs = _retrieve_wind_volume_to_target_grid(
        prd,
        target_x=target_x,
        target_y=target_y,
        level_heights=level_heights,
        sweeps=sweeps,
        field_name=field_name,
        range_mode=range_mode,
        az_num=az_num,
        bin_num=bin_num,
        azimuth_step=azimuth_step,
        range_step=range_step,
        max_range_km=max_range_km,
        fillvalue=fillvalue,
        min_valid_fraction=min_valid_fraction,
        min_valid_count=min_valid_count,
        min_azimuth_coverage_deg=min_azimuth_coverage_deg,
        min_azimuth_sector_count=min_azimuth_sector_count,
        max_condition_number=max_condition_number,
        horizontal_radius_m=horizontal_radius_m,
        max_horizontal_radius_m=max_horizontal_radius_m,
        horizontal_min_neighbors=horizontal_min_neighbors,
        vertical_tolerance_m=vertical_tolerance_m,
        max_vertical_gap_m=max_vertical_gap_m,
        workers=workers,
        max_selected_sweeps=max_selected_sweeps,
        sweep_min_valid_fraction=sweep_min_valid_fraction,
        sweep_min_valid_count=sweep_min_valid_count,
        sweep_min_azimuth_sector_count=sweep_min_azimuth_sector_count,
        sweep_max_median_fit_rmse=sweep_max_median_fit_rmse,
        sweep_max_condition_number=sweep_max_condition_number,
    )
    attrs = dict(outputs["attrs"])
    attrs["grid_definition"] = "xy"
    return _build_wind_volume_dataset_xy(
        x_axis=x_axis,
        y_axis=y_axis,
        level_heights=level_heights,
        u_volume=outputs["u"],
        v_volume=outputs["v"],
        fit_rmse_volume=outputs["fit_rmse"],
        valid_count_volume=outputs["valid_count"],
        source_sweep_count_volume=outputs["source_sweep_count"],
        neighbor_count_volume=outputs["neighbor_count"],
        effective_vertical_support_volume=outputs["effective_vertical_support"],
        sweep_spread_m_volume=outputs["sweep_spread_m"],
        quality_score_volume=outputs["quality_score"],
        quality_flag_volume=outputs["quality_flag"],
        attrs=attrs,
    )


def retrieve_wind_volume_lonlat(
    prd,
    XLon,
    YLat,
    level_heights,
    sweeps=None,
    field_name=None,
    range_mode="aligned",
    **kwargs
):
    """Retrieve a gridded 3-D horizontal wind volume on a lon/lat grid."""
    x_lon = _as_1d_float(XLon)
    y_lat = _as_1d_float(YLat)
    level_heights = _normalize_level_heights(level_heights)
    if x_lon.size == 0 or y_lat.size == 0:
        raise ValueError("XLon and YLat must contain at least one value")
    grid_lon, grid_lat = np.meshgrid(x_lon, y_lat, indexing="ij")
    proj = prd._get_local_projection()
    grid_x, grid_y = proj(grid_lon, grid_lat, inverse=False)
    outputs = _retrieve_wind_volume_to_target_grid(
        prd,
        target_x=np.asarray(grid_x, dtype=np.float64),
        target_y=np.asarray(grid_y, dtype=np.float64),
        level_heights=level_heights,
        sweeps=sweeps,
        field_name=field_name,
        range_mode=range_mode,
        **kwargs
    )
    attrs = dict(outputs["attrs"])
    attrs["grid_definition"] = "lonlat"
    dataset = xr.Dataset(
        data_vars={
            "u": (("z_wind_geo", "lon_wind", "lat_wind"), np.asarray(outputs["u"], dtype=np.float64)),
            "v": (("z_wind_geo", "lon_wind", "lat_wind"), np.asarray(outputs["v"], dtype=np.float64)),
            "wind_speed": (
                ("z_wind_geo", "lon_wind", "lat_wind"),
                np.hypot(np.asarray(outputs["u"], dtype=np.float64), np.asarray(outputs["v"], dtype=np.float64)),
            ),
            "wind_direction": (
                ("z_wind_geo", "lon_wind", "lat_wind"),
                _wind_direction_from_uv(np.asarray(outputs["u"], dtype=np.float64), np.asarray(outputs["v"], dtype=np.float64)),
            ),
            "fit_rmse": (
                ("z_wind_geo", "lon_wind", "lat_wind"),
                np.asarray(outputs["fit_rmse"], dtype=np.float64),
            ),
            "valid_count": (
                ("z_wind_geo", "lon_wind", "lat_wind"),
                np.asarray(outputs["valid_count"], dtype=np.int32),
            ),
            "source_sweep_count": (
                ("z_wind_geo", "lon_wind", "lat_wind"),
                np.asarray(outputs["source_sweep_count"], dtype=np.int32),
            ),
            "neighbor_count": (
                ("z_wind_geo", "lon_wind", "lat_wind"),
                np.asarray(outputs["neighbor_count"], dtype=np.int32),
            ),
            "effective_vertical_support": (
                ("z_wind_geo", "lon_wind", "lat_wind"),
                np.asarray(outputs["effective_vertical_support"], dtype=np.int32),
            ),
            "sweep_spread_m": (
                ("z_wind_geo", "lon_wind", "lat_wind"),
                np.asarray(outputs["sweep_spread_m"], dtype=np.float64),
            ),
            "quality_score": (
                ("z_wind_geo", "lon_wind", "lat_wind"),
                np.asarray(outputs["quality_score"], dtype=np.float64),
            ),
            "quality_flag": (
                ("z_wind_geo", "lon_wind", "lat_wind"),
                np.asarray(outputs["quality_flag"], dtype=np.uint8),
            ),
        },
        coords={
            "z_wind_geo": level_heights,
            "lon_wind": x_lon,
            "lat_wind": y_lat,
        },
        attrs=attrs,
    )
    dataset["z_wind_geo"].attrs = {
        "units": "meters",
        "standard_name": "height",
        "long_name": "target wind-volume height above mean sea level",
        "positive": "up",
    }
    dataset["lon_wind"].attrs = {
        "units": "degrees",
        "standard_name": "longitude",
        "long_name": "wind-volume longitude",
        "axis": "lonlat_coordinate",
    }
    dataset["lat_wind"].attrs = {
        "units": "degrees",
        "standard_name": "latitude",
        "long_name": "wind-volume latitude",
        "axis": "lonlat_coordinate",
    }
    dataset["u"].attrs = {
        "units": "meters_per_second",
        "standard_name": "eastward_wind_component",
        "long_name": "gridded eastward wind component from single-radar VVP",
    }
    dataset["v"].attrs = {
        "units": "meters_per_second",
        "standard_name": "northward_wind_component",
        "long_name": "gridded northward wind component from single-radar VVP",
    }
    dataset["wind_speed"].attrs = {
        "units": "meters_per_second",
        "long_name": "gridded horizontal wind speed",
    }
    dataset["wind_direction"].attrs = {
        "units": "degrees",
        "long_name": "gridded meteorological wind direction",
    }
    dataset["fit_rmse"].attrs = {
        "units": "meters_per_second",
        "long_name": "gridded VVP fit root-mean-square residual",
    }
    dataset["valid_count"].attrs = {
        "units": "count",
        "long_name": "effective valid sample count contributing to the target voxel",
    }
    dataset["source_sweep_count"].attrs = {
        "units": "count",
        "long_name": "number of sweeps contributing to the target voxel",
    }
    dataset["neighbor_count"].attrs = {
        "units": "count",
        "long_name": "effective horizontal neighbors contributing to the target voxel",
    }
    dataset["effective_vertical_support"].attrs = {
        "units": "count",
        "long_name": "number of vertically supporting sweep samples used at the target voxel",
    }
    dataset["sweep_spread_m"].attrs = {
        "units": "meters",
        "long_name": "vertical spread of sweep samples contributing to the target voxel",
    }
    dataset["quality_score"].attrs = {
        "units": "1",
        "long_name": "heuristic voxel quality score from 0 to 100",
    }
    dataset["quality_flag"].attrs = {
        "units": "1",
        "long_name": "voxel quality flag",
        "flag_values": np.array([0, 1, 2, 3], dtype=np.uint8),
        "flag_meanings": "invalid low_confidence usable high_confidence",
    }
    return dataset


def vad(azimuth, elevation, PPI_Vr, fillvalue=-999.0):
    """
    Backward-compatible single-ring VAD helper.

    Returns
    -------
    tuple[float, float]
        Retrieved ``(u, v)`` wind components.
    """
    result = fit_vad_ring(azimuth, elevation, PPI_Vr, fillvalue=fillvalue)
    if not np.isfinite(result["u"]) or not np.isfinite(result["v"]):
        raise ValueError("VAD retrieval failed because the ring does not satisfy the fit criteria")
    return result["u"], result["v"]


def VAD(azimuth, elevation, PPI_Vr, fillvalue=-999.0):
    """Backward-compatible alias for ``vad``."""
    return vad(azimuth, elevation, PPI_Vr, fillvalue=fillvalue)


def vvp(azimuth, elevation, PPI_Vr, az_num, bin_num, fillvalue=-999.0):
    """
    Backward-compatible local VVP helper.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Retrieved eastward and northward wind components.
    """
    outputs = _run_local_vvp(
        azimuth,
        elevation,
        PPI_Vr,
        az_num=az_num,
        bin_num=bin_num,
        fillvalue=fillvalue,
    )
    wind_u = np.where(np.isfinite(outputs["u"]), outputs["u"], float(fillvalue)).astype(np.float32)
    wind_v = np.where(np.isfinite(outputs["v"]), outputs["v"], float(fillvalue)).astype(np.float32)
    return wind_u, wind_v


def VVP(azimuth, elevation, PPI_Vr, az_num, bin_num, fillvalue=-999.0):
    """Backward-compatible alias for ``vvp``."""
    return vvp(azimuth, elevation, PPI_Vr, az_num, bin_num, fillvalue=fillvalue)
