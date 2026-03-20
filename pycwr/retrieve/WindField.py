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

import json

import numpy as np
import xarray as xr

from ..configure.default_config import DEFAULT_METADATA


WIND_REFERENCE_NOTES = {
    "VAD": [
        "Browning and Wexler (1968), Journal of Applied Meteorology, doi:10.1175/1520-0450(1968)007<0105:TDOKPO>2.0.CO;2",
    ],
    "VVP": [
        "Waldteufel and Corbin (1979), Journal of Applied Meteorology, doi:10.1175/1520-0450(1979)018<0532:OTAOSD>2.0.CO;2",
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
        normal = design_weighted.T.dot(design_weighted)
        condition_number = float(np.linalg.cond(normal))
        if (not np.isfinite(condition_number)) or condition_number > float(max_condition_number):
            return None, condition_number
        try:
            coefficients = np.linalg.solve(normal, design_weighted.T.dot(response_weighted))
        except np.linalg.LinAlgError:
            return None, condition_number
        return coefficients, condition_number

    coefficients, condition_number = _solve_once(weight_valid)
    if coefficients is None:
        return {
            "coefficients": None,
            "condition_number": condition_number,
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
        new_coefficients, new_condition_number = _solve_once(updated_weights)
        if new_coefficients is None:
            break
        if np.allclose(new_coefficients, coefficients, rtol=1.0e-5, atol=1.0e-5):
            coefficients = new_coefficients
            condition_number = new_condition_number
            robust_weights = updated_weights
            break
        coefficients = new_coefficients
        condition_number = new_condition_number
        robust_weights = updated_weights
    fitted = design_valid.dot(coefficients)
    rmse = float(np.sqrt(np.mean((response_valid - fitted) ** 2))) if response_valid.size else np.nan
    return {
        "coefficients": coefficients,
        "condition_number": condition_number,
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
        "valid_count": int(valid_count),
        "valid_fraction": float(valid_fraction),
        "azimuth_coverage_deg": float(azimuth_coverage_deg),
    }


def fit_vad_ring(
    azimuth,
    elevation,
    radial_velocity,
    fillvalue=-999.0,
    min_valid_fraction=0.5,
    min_valid_count=16,
    min_azimuth_coverage_deg=240.0,
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
    if (
        valid_count < int(min_valid_count)
        or valid_fraction < float(min_valid_fraction)
        or azimuth_coverage < float(min_azimuth_coverage_deg)
    ):
        return _empty_fit_result(valid_count, valid_fraction, azimuth_coverage)

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
        "valid_count": valid_count,
        "valid_fraction": valid_fraction,
        "azimuth_coverage_deg": azimuth_coverage,
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
        "valid_count": np.zeros((naz, nrange), dtype=np.int32),
    }

    for az_index in range(naz):
        az_slice = slice(az_index, az_index + az_num)
        az_window = padded_azimuth[az_slice]
        el_window = padded_elevation[az_slice]
        for range_index in range(half_bin, nrange - half_bin):
            window = padded_velocity[az_slice, range_index - half_bin : range_index + half_bin + 1]
            valid = np.isfinite(window)
            valid_count = int(valid.sum())
            valid_fraction = float(valid_count / window.size)
            azimuth_coverage = _azimuth_coverage_deg(np.repeat(az_window, bin_num)[valid.reshape(-1)])
            outputs["valid_count"][az_index, range_index] = valid_count
            outputs["valid_fraction"][az_index, range_index] = valid_fraction
            outputs["azimuth_coverage_deg"][az_index, range_index] = azimuth_coverage
            if (
                valid_count < int(min_valid_count)
                or valid_fraction < float(min_valid_fraction)
                or azimuth_coverage < float(min_azimuth_coverage_deg)
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
    dataset["valid_count"].attrs = {"long_name": "number of valid samples in the local retrieval window"}
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
