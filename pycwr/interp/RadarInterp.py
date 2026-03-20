from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timedelta, timezone
import json
import multiprocessing as mp
import os
from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import BoundaryNorm, ListedColormap
try:
    import tomllib
except ImportError:  # pragma: no cover
    tomllib = None
try:
    from zoneinfo import ZoneInfo
except ImportError:  # pragma: no cover
    ZoneInfo = None

import numpy as np
import pyproj
import xarray as xr
from scipy import spatial

from ..configure.config import cfg
from ..configure.default_config import CINRAD_COLORMAP, CINRAD_field_normvar
from ..core.RadarProduct import (
    PRODUCT_REFERENCE_NOTES,
    derive_cr as _derive_cr_impl,
    derive_et as _derive_et_impl,
    derive_vil as _derive_vil_impl,
)
from ..core.transforms import cartesian_xyz_to_antenna
from ..io import read_auto
from ..io.util import get_radar_info

try:
    from ..core.RadarGridC import get_CAPPI_3d
except ImportError:
    from ..core.RadarGrid import get_CAPPI_3d


_FILENAME_TIME_RE = re.compile(r"(\d{14})")
_REFERENCE_CANVAS_PX = (920, 790)
_REFERENCE_CR_BOUNDS = np.arange(-5.0, 70.0, 5.0, dtype=np.float64)
_REFERENCE_CR_TICKS = _REFERENCE_CR_BOUNDS.copy()
_REFERENCE_CR_COLORS = [
    (0 / 255.0, 172 / 255.0, 164 / 255.0),
    (192 / 255.0, 192 / 255.0, 254 / 255.0),
    (122 / 255.0, 114 / 255.0, 238 / 255.0),
    (30 / 255.0, 38 / 255.0, 208 / 255.0),
    (166 / 255.0, 252 / 255.0, 168 / 255.0),
    (0 / 255.0, 234 / 255.0, 0 / 255.0),
    (16 / 255.0, 146 / 255.0, 26 / 255.0),
    (252 / 255.0, 244 / 255.0, 100 / 255.0),
    (200 / 255.0, 200 / 255.0, 2 / 255.0),
    (140 / 255.0, 140 / 255.0, 0 / 255.0),
    (254 / 255.0, 172 / 255.0, 172 / 255.0),
    (254 / 255.0, 100 / 255.0, 92 / 255.0),
    (238 / 255.0, 2 / 255.0, 48 / 255.0),
    (212 / 255.0, 142 / 255.0, 254 / 255.0),
]
_REFERENCE_FONT_CANDIDATES = (
    "Noto Sans CJK JP",
    "Noto Sans CJK SC",
    "SimHei",
    "Microsoft YaHei",
    "WenQuanYi Zen Hei",
    "Arial Unicode MS",
    "DejaVu Sans",
)
_PRODUCT_REFERENCE_NOTES = PRODUCT_REFERENCE_NOTES

def get_weight(dist, r, method="barnes"):
    """
    Return interpolation weights for the requested method.
    :param dist: distance from source points to the target point.
    :param r: influence radius.
    :param method: interpolation method.
    :return: weight for each source point.
    """
    if method == "barnes":
        weight = np.exp(-4*dist**2/r**2)
    elif method == "cressman":
        weight = (r ** 2 - dist ** 2) / (dist ** 2 + r ** 2)
    else:
        raise Exception("Unidentified method!, must be cressman, barnes")
    return weight

def _get_interp_around_point(point_old, point_new, around_r):
    """
    Find all source points within ``around_r`` of each target point.
    :param point_old: source coordinates, ``np.c_[x, y, z]`` or ``np.c_[x, y]``.
    :param point_new: target coordinates, ``np.c_[x', y', z']`` or ``np.c_[x', y']``.
    :param around_r: search radius in meters.
    """
    kdtree = spatial.cKDTree(point_old)
    index_nearest = kdtree.query_ball_point(point_new, around_r)
    dist = [np.sqrt(np.sum(np.square(point_old[i,...] - itarget), axis=1)) for i, itarget \
            in zip(index_nearest, point_new)]
    return index_nearest, dist

def radar_interp2d(points, values, xi, around_r,  influence_radius=None, method="barnes", fill_value=np.nan):
    """
    Interpolate unstructured D-dimensional data.
    Parameters
    ----------
    points : ndarray of floats, shape (n, D)
        Data point coordinates. Can either be an array of
        shape (n, D), or a tuple of `ndim` arrays.
    values : ndarray of float or complex, shape (n,)
        Data values.
    xi : 2-D ndarray of float or tuple of 1-D array, shape (M, D)
        Points at which to interpolate data.
    around_r: interpolate from source points within ``around_r``.
    influence_radius: influence radius used by the weighting function.
    method : {'barnes', 'cressman'}
        Method of interpolation. One of 'barnes', 'cressman'
    fill_value : float, optional
        Value used to fill in for requested points outside of the
        convex hull of the input points.  If not provided, then the
        default is ``nan``. This option has no effect for the
        'nearest' method.
    Returns
    -------
    ndarray
        Array of interpolated values.
    """
    if influence_radius is None:
        influence_radius = around_r
    grid_shape = xi[0].shape
    target = np.column_stack([xi_grid.ravel() for xi_grid in xi])
    index, distance = _get_interp_around_point(points, target, around_r)
    nrows, _ = target.shape
    grid_vals = np.empty(nrows)
    for i in range(nrows):
        if index[i]:
            weight = get_weight(distance[i], influence_radius, method=method)
            grid_vals[i] = np.dot(values[index[i]], weight)/np.sum(weight)
        else:
            grid_vals[i] = fill_value
    return grid_vals.reshape(grid_shape)

def _get_interp_around_point_var(point_old, point_new, bandwidth=1):
    """
    Find neighbors with a range-dependent radius of influence.
    :param point_old:
    :param point_new:
    :param bandwidth: beam width in degrees.
    """
    min_roi = cfg.interp.mroi
    kdtree = spatial.cKDTree(point_old)
    index_nearest = []
    nrows = point_new.shape[0]
    roi = np.empty(nrows)
    for i, itarget in enumerate(point_new):
        roi[i] = min_roi + ((itarget[0] / 1000.) ** 2 + (itarget[1] / 1000.) ** 2) ** 0.5 * bandwidth * cfg.interp.coeff
        index_nearest.append(kdtree.query_ball_point(itarget, roi[i]))
    distance = [np.sqrt(np.sum(np.square(point_old[i, :] - j), axis=1)) for i, j \
            in zip(index_nearest, point_new)]
    return index_nearest, distance, roi

def radar_interp2d_var(points, values, xi, bandwidth=1, method="barnes", fill_value=np.nan):
    """
    Interpolate with a radius of influence that increases with range.
    :param points: source points.
    :param values:
    :param xi:
    :param bandwidth: beam width.
    :param method:
    :param fill_value:
    :return:
    """

    grid_shape = xi[0].shape
    target = np.column_stack([xi_grid.ravel() for xi_grid in xi])
    index, distance, roi = _get_interp_around_point_var(points, target, bandwidth=bandwidth)
    nrows, _ = target.shape
    grid_vals = np.empty(nrows)
    for i in range(nrows):
        if index[i]:
            weight = get_weight(distance[i], roi[i], method=method)
            grid_vals[i] = np.dot(values[index[i]], weight) / np.sum(weight)
        else:
            grid_vals[i] = fill_value
    return grid_vals.reshape(grid_shape)


def _coerce_field_names(field_names, network_cfg=None):
    """Normalize field names to a non-empty list of strings."""
    if field_names is None:
        field_names = _resolve_network_value("field_names", None, network_cfg=network_cfg)
    if field_names is None:
        return ["dBZ"]
    if isinstance(field_names, str):
        field_names = [field_names]
    field_names = [str(name) for name in field_names if str(name)]
    if not field_names:
        raise ValueError("field_names must contain at least one field.")
    return field_names


def _resolve_network_value(name, value, network_cfg=None):
    """Read a network parameter from config when not provided explicitly."""
    if value is not None:
        return value
    if network_cfg is None:
        network_cfg = getattr(cfg, "network", None)
    if network_cfg is None:
        return None
    if isinstance(network_cfg, dict):
        return network_cfg.get(name)
    return getattr(network_cfg, name, None)


def _read_network_config_file(config_path):
    """Load 3D radar network options from a JSON/TOML/YAML file."""
    config_path = Path(config_path)
    suffix = config_path.suffix.lower()
    if suffix == ".json":
        with config_path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
    elif suffix == ".toml":
        if tomllib is None:
            raise ValueError("TOML configuration requires Python 3.11+ or tomllib.")
        with config_path.open("rb") as handle:
            data = tomllib.load(handle)
    elif suffix in (".yml", ".yaml"):
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover
            raise ValueError("YAML configuration requires PyYAML to be installed.") from exc
        with config_path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
    else:
        raise ValueError("Unsupported config file suffix %s; use .json, .toml, .yml, or .yaml." % suffix)
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError("Network config file must contain a mapping/object.")
    if "network" in data:
        data = data["network"]
    if not isinstance(data, dict):
        raise ValueError("The 'network' section must contain a mapping/object.")
    return data


def load_network_config(config_path, inplace=False):
    """Load a 3D radar network configuration file."""
    network_config = _read_network_config_file(config_path)
    if inplace:
        cfg.network.update(network_config)
    return network_config


def _coerce_field_range_mode(field_names, field_range_mode, network_cfg=None):
    """Resolve per-field range extraction modes."""
    field_names = _coerce_field_names(field_names, network_cfg=network_cfg)
    if field_range_mode is None:
        field_range_mode = _resolve_network_value("field_range_mode", None, network_cfg=network_cfg)
    field_modes = {field_name: "aligned" for field_name in field_names}
    if isinstance(field_range_mode, str):
        for field_name in field_names:
            field_modes[field_name] = field_range_mode
    elif field_range_mode:
        for field_name in field_names:
            if field_name in field_range_mode:
                field_modes[field_name] = str(field_range_mode[field_name])
    if "dBZ" in field_modes and field_modes["dBZ"] == "aligned":
        field_modes["dBZ"] = "native"
    for field_name, range_mode in field_modes.items():
        if range_mode not in ("aligned", "native"):
            raise ValueError("range_mode for %s must be 'aligned' or 'native'." % field_name)
    return field_modes


def _resolve_qc_field_name(field_name, use_qc):
    """Map requested fields to QC-corrected field names when QC is enabled."""
    if not use_qc:
        return field_name
    corrected = {
        "dBZ": "Zc",
        "ZDR": "ZDRc",
        "KDP": "KDPc",
    }
    return corrected.get(field_name, field_name)


def _normalize_blind_method(blind_method):
    blind_method = "mask" if blind_method is None else str(blind_method).lower()
    valid_methods = {"mask", "nearest_gate", "lowest_sweep", "hybrid", "nearest_valid"}
    if blind_method not in valid_methods:
        raise ValueError(
            "blind_method must be one of 'mask', 'nearest_gate', 'lowest_sweep', "
            "'hybrid', or 'nearest_valid'."
        )
    return blind_method


def _radar_supports_dualpol_qc(prd):
    """Return whether a PRD volume carries the moments required by dual-pol QC."""
    for sweep in range(prd.nsweeps):
        dataset = prd.fields[sweep]
        if "dBZ" not in dataset:
            return False
        if ("KDP" not in dataset) and ("PhiDP" not in dataset):
            return False
    return True


def _get_prd_station_info(prd):
    """Extract radar site latitude, longitude, and altitude from PRD metadata."""
    scan_info = getattr(prd, "scan_info", None)
    if scan_info is None:
        raise ValueError("Radar object does not provide scan_info.")
    altitude = float(scan_info["altitude"].values) if "altitude" in scan_info else 0.0
    return (
        float(scan_info["latitude"].values),
        float(scan_info["longitude"].values),
        altitude,
    )


def parse_radar_time_from_filename(path):
    """Extract the scan time from a radar filename."""
    match = _FILENAME_TIME_RE.search(Path(path).name)
    if match is None:
        raise ValueError("Unable to parse scan time from filename: %s" % path)
    return datetime.strptime(match.group(1), "%Y%m%d%H%M%S")


def discover_radar_files(radar_dirs, pattern="*.bin*"):
    """Collect candidate radar files under each configured radar directory."""
    if not radar_dirs:
        raise ValueError("radar_dirs must contain at least one directory.")
    radar_files = {}
    for radar_dir in radar_dirs:
        directory = Path(radar_dir)
        if not directory.is_dir():
            continue
        matches = sorted(path for path in directory.glob(pattern) if path.is_file())
        if matches:
            radar_files[directory.name] = matches
    return radar_files


def select_radar_files(radar_dirs, target_time, tolerance_minutes=10, pattern="*.bin*"):
    """Select the nearest-in-time file from each radar directory."""
    if target_time is None:
        raise ValueError("target_time is required.")
    if isinstance(target_time, str):
        target_time = datetime.fromisoformat(target_time)
    tolerance_seconds = float(tolerance_minutes) * 60.0
    selected = []
    for radar_id, files in discover_radar_files(radar_dirs, pattern=pattern).items():
        best_path = None
        best_delta = None
        best_time = None
        for path in files:
            scan_time = parse_radar_time_from_filename(path)
            delta_seconds = abs((scan_time - target_time).total_seconds())
            if (best_delta is None) or (delta_seconds < best_delta):
                best_delta = delta_seconds
                best_path = path
                best_time = scan_time
        if best_path is not None and best_delta <= tolerance_seconds:
            selected.append(
                {
                    "radar_id": radar_id,
                    "path": str(best_path),
                    "scan_time": best_time,
                    "delta_seconds": best_delta,
                }
            )
    if not selected:
        raise ValueError("No radar files matched the target time within tolerance.")
    return selected


def build_latlon_grid(lon_min, lon_max, lat_min, lat_max, lon_res_deg, lat_res_deg):
    """Build a regular longitude/latitude target grid."""
    lon_res_deg = float(lon_res_deg)
    lat_res_deg = float(lat_res_deg)
    if lon_res_deg <= 0 or lat_res_deg <= 0:
        raise ValueError("lon_res_deg and lat_res_deg must be positive.")
    if lon_max < lon_min or lat_max < lat_min:
        raise ValueError("Grid bounds must satisfy min <= max.")
    lon = np.arange(float(lon_min), float(lon_max) + lon_res_deg * 0.5, lon_res_deg, dtype=np.float64)
    lat = np.arange(float(lat_min), float(lat_max) + lat_res_deg * 0.5, lat_res_deg, dtype=np.float64)
    grid_lon, grid_lat = np.meshgrid(lon, lat, indexing="xy")
    return lon, lat, grid_lon, grid_lat


def _coerce_output_products(output_products, network_cfg=None):
    output_products = _resolve_network_value("output_products", output_products, network_cfg=network_cfg)
    if output_products is None:
        return ["CR", "VIL", "ET"]
    if isinstance(output_products, str):
        output_products = [output_products]
    valid = []
    for product_name in output_products:
        product_name = str(product_name).upper()
        if product_name not in {"CR", "VIL", "ET"}:
            raise ValueError("Unsupported output product %s." % product_name)
        if product_name not in valid:
            valid.append(product_name)
    return valid


def _coerce_plot_height_levels(level_heights, plot_height_levels, network_cfg=None):
    plot_height_levels = _resolve_network_value("plot_height_levels", plot_height_levels, network_cfg=network_cfg)
    if level_heights is None:
        return np.asarray([], dtype=np.float64)
    if plot_height_levels is None or len(plot_height_levels) == 0:
        return np.asarray(level_heights, dtype=np.float64)
    return np.asarray(plot_height_levels, dtype=np.float64)


def _coerce_product_level_heights(level_heights, product_level_heights, network_cfg=None):
    product_level_heights = _resolve_network_value(
        "product_level_heights",
        product_level_heights,
        network_cfg=network_cfg,
    )
    if level_heights is None:
        return np.asarray([], dtype=np.float64)
    level_heights = np.asarray(level_heights, dtype=np.float64)
    if level_heights.size == 0:
        return level_heights
    if product_level_heights is None or len(product_level_heights) == 0:
        if level_heights.size == 1:
            return level_heights.copy()
        level_min = float(np.min(level_heights))
        level_max = float(np.max(level_heights))
        return np.arange(level_min, level_max + 500.0 * 0.5, 500.0, dtype=np.float64)
    return np.asarray(product_level_heights, dtype=np.float64)


def _coerce_quality_weight_terms(quality_weight_terms, network_cfg=None):
    quality_weight_terms = _resolve_network_value(
        "quality_weight_terms",
        quality_weight_terms,
        network_cfg=network_cfg,
    )
    if quality_weight_terms is None:
        quality_weight_terms = ["range", "beam", "vertical", "qc", "blockage"]
    if isinstance(quality_weight_terms, str):
        quality_weight_terms = [quality_weight_terms]
    valid_terms = []
    for term in quality_weight_terms:
        normalized = str(term).strip().lower()
        if normalized not in {"range", "beam", "vertical", "qc", "blockage"}:
            raise ValueError("Unsupported quality weight term %s." % term)
        if normalized not in valid_terms:
            valid_terms.append(normalized)
    return valid_terms


def _validate_network_field_names(field_names):
    unsupported = {"V", "VR", "VELOCITY"}
    for field_name in field_names:
        if str(field_name).upper() in unsupported:
            raise ValueError(
                "3D network mosaic does not currently support velocity fields; "
                "remove %s from field_names." % field_name
            )


def _resolve_plot_output_dir(plot_output_dir, output_path, output_dir):
    if plot_output_dir is not None:
        return str(plot_output_dir)
    if output_path is not None:
        return str(Path(output_path).parent)
    if output_dir is not None:
        return str(output_dir)
    return None


def _dBZ_to_linear(dbz):
    return np.power(10.0, dbz / 10.0)


def _derive_cr(volume, fillvalue=-999.0):
    return _derive_cr_impl(volume, fillvalue=fillvalue)


def _derive_vil(volume, level_heights, fillvalue=-999.0, min_dbz=18.0, max_dbz_cap=56.0):
    """
    Derive VIL from a reflectivity volume.

    The liquid-water relation follows Greene and Clark (1972) and the common
    WDTD/operational VIL practice of capping reflectivity at 56 dBZ. The lower
    reflectivity cutoff remains a configurable pycwr assumption and should be
    treated as a workflow parameter rather than a universal standard.
    """
    return _derive_vil_impl(
        volume,
        level_heights,
        fillvalue=fillvalue,
        min_dbz=min_dbz,
        max_dbz_cap=max_dbz_cap,
    )


def _derive_et(volume, level_heights, fillvalue=-999.0, threshold_dbz=18.0, return_topped=False):
    """
    Derive echo-top height from a reflectivity volume.

    The implementation follows the standard threshold-crossing idea used by
    operational ET products: locate the highest level meeting the threshold and
    linearly interpolate to the crossing height when the next level drops below
    threshold. When the highest available level still exceeds the threshold, a
    topped flag is reported to distinguish "top at last sampled level" from a
    fully resolved echo top.
    """
    return _derive_et_impl(
        volume,
        level_heights,
        fillvalue=fillvalue,
        threshold_dbz=threshold_dbz,
        return_topped=return_topped,
    )


def _get_sweep_beam_widths(prd, sweep_count):
    try:
        beam_widths = np.asarray(prd.scan_info["beam_width"].values, dtype=np.float64)
    except Exception:
        beam_widths = np.full(int(sweep_count), 1.0, dtype=np.float64)
    if beam_widths.size == 0:
        return np.full(int(sweep_count), 1.0, dtype=np.float64)
    if beam_widths.size < int(sweep_count):
        pad_value = float(beam_widths[-1])
        beam_widths = np.pad(beam_widths, (0, int(sweep_count) - beam_widths.size), constant_values=pad_value)
    return np.where(np.isfinite(beam_widths[:int(sweep_count)]), beam_widths[:int(sweep_count)], 1.0)


def _compute_observability_mask(
    vol_range,
    fix_elevation,
    radar_height,
    grid_x,
    grid_y,
    level_heights,
    radar_x,
    radar_y,
    beam_widths=None,
    effective_earth_radius=None,
):
    local_x = np.asarray(grid_x, dtype=np.float64) - float(radar_x)
    local_y = np.asarray(grid_y, dtype=np.float64) - float(radar_y)
    fix_elevation = np.asarray(fix_elevation, dtype=np.float64)
    level_heights = np.asarray(level_heights, dtype=np.float64)
    mask = np.zeros((level_heights.size,) + local_x.shape, dtype=bool)
    if fix_elevation.size == 0:
        return mask
    if beam_widths is None:
        beam_widths = np.full(fix_elevation.size, 1.0, dtype=np.float64)
    beam_widths = np.asarray(beam_widths, dtype=np.float64)
    if beam_widths.size < fix_elevation.size:
        beam_widths = np.pad(beam_widths, (0, fix_elevation.size - beam_widths.size), constant_values=float(beam_widths[-1]))
    for ilevel, level_height in enumerate(level_heights):
        _, ranges, elevation = cartesian_xyz_to_antenna(
            local_x,
            local_y,
            float(level_height),
            float(radar_height),
            effective_earth_radius=effective_earth_radius,
        )
        inside = np.isfinite(ranges) & np.isfinite(elevation)
        if fix_elevation.size == 1:
            beam_half = max(float(beam_widths[0]) * 0.5, 0.1)
            first_gate = float(vol_range[0][0])
            last_gate = float(vol_range[0][-1])
            mask[ilevel] = (
                inside
                & (np.abs(elevation - fix_elevation[0]) <= beam_half)
                & (ranges >= first_gate)
                & (ranges <= last_gate)
            )
            continue
        inside &= elevation >= fix_elevation[0]
        inside &= elevation <= fix_elevation[-1]
        if not np.any(inside):
            continue
        upper_index = np.searchsorted(fix_elevation, elevation, side="right")
        upper_index = np.clip(upper_index, 1, fix_elevation.size - 1)
        lower_index = upper_index - 1
        first_gate = np.full_like(ranges, np.nan, dtype=np.float64)
        last_gate = np.full_like(ranges, np.nan, dtype=np.float64)
        for sweep_index in range(fix_elevation.size - 1):
            sweep_mask = inside & (lower_index == sweep_index)
            if not np.any(sweep_mask):
                continue
            first_gate[sweep_mask] = max(float(vol_range[sweep_index][0]), float(vol_range[sweep_index + 1][0]))
            last_gate[sweep_mask] = min(float(vol_range[sweep_index][-1]), float(vol_range[sweep_index + 1][-1]))
        mask[ilevel] = inside & np.isfinite(first_gate) & np.isfinite(last_gate) & (ranges >= first_gate) & (ranges <= last_gate)
    return mask


def _compute_quality_weight_volume(
    vol_range,
    fix_elevation,
    radar_height,
    grid_x,
    grid_y,
    level_heights,
    radar_x,
    radar_y,
    beam_widths=None,
    effective_earth_radius=None,
    quality_weight_terms=None,
    range_decay_m=230000.0,
    beam_radius_ref_m=3000.0,
    vertical_gap_ref_deg=0.5,
):
    """
    Build per-radar geometric quality weights for network compositing.

    This is a pycwr quality-weight implementation inspired by the quality-index
    family used in operational mosaics. It is intentionally more explicit than
    the Chapter 4 experiments summarized in
    "吴翀;双偏振雷达的资料质量分析,相态识別及组网应用[D];南京信息工程大学;2018年":
    the thesis discusses exponential and quality-based overlap handling, while
    pycwr combines range, beam broadening, vertical representativeness,
    optional QC masks, and optional blockage weights into one reusable weight
    volume.
    """
    quality_weight_terms = set(quality_weight_terms or [])
    local_x = np.asarray(grid_x, dtype=np.float64) - float(radar_x)
    local_y = np.asarray(grid_y, dtype=np.float64) - float(radar_y)
    fix_elevation = np.asarray(fix_elevation, dtype=np.float64)
    level_heights = np.asarray(level_heights, dtype=np.float64)
    weight = np.zeros((level_heights.size,) + local_x.shape, dtype=np.float64)
    if fix_elevation.size == 0:
        return weight
    if beam_widths is None:
        beam_widths = np.full(fix_elevation.size, 1.0, dtype=np.float64)
    beam_widths = np.asarray(beam_widths, dtype=np.float64)
    if beam_widths.size < fix_elevation.size:
        beam_widths = np.pad(beam_widths, (0, fix_elevation.size - beam_widths.size), constant_values=float(beam_widths[-1]))
    for ilevel, level_height in enumerate(level_heights):
        _, ranges, elevation = cartesian_xyz_to_antenna(
            local_x,
            local_y,
            float(level_height),
            float(radar_height),
            effective_earth_radius=effective_earth_radius,
        )
        observable = np.isfinite(ranges) & np.isfinite(elevation)
        if fix_elevation.size == 1:
            beam_half = max(float(beam_widths[0]) * 0.5, 0.1)
            observable &= np.abs(elevation - fix_elevation[0]) <= beam_half
            observable &= ranges >= float(vol_range[0][0])
            observable &= ranges <= float(vol_range[0][-1])
            representative_beam = np.full_like(ranges, float(beam_widths[0]), dtype=np.float64)
            gap_half = np.zeros_like(ranges, dtype=np.float64)
        else:
            observable &= elevation >= fix_elevation[0]
            observable &= elevation <= fix_elevation[-1]
            upper_index = np.searchsorted(fix_elevation, elevation, side="right")
            upper_index = np.clip(upper_index, 1, fix_elevation.size - 1)
            lower_index = upper_index - 1
            first_gate = np.full_like(ranges, np.nan, dtype=np.float64)
            last_gate = np.full_like(ranges, np.nan, dtype=np.float64)
            representative_beam = np.full_like(ranges, np.nan, dtype=np.float64)
            gap_half = np.full_like(ranges, np.nan, dtype=np.float64)
            for sweep_index in range(fix_elevation.size - 1):
                sweep_mask = observable & (lower_index == sweep_index)
                if not np.any(sweep_mask):
                    continue
                first_gate[sweep_mask] = max(float(vol_range[sweep_index][0]), float(vol_range[sweep_index + 1][0]))
                last_gate[sweep_mask] = min(float(vol_range[sweep_index][-1]), float(vol_range[sweep_index + 1][-1]))
                representative_beam[sweep_mask] = 0.5 * (
                    float(beam_widths[sweep_index]) + float(beam_widths[sweep_index + 1])
                )
                gap_half[sweep_mask] = 0.5 * (
                    float(fix_elevation[sweep_index + 1]) - float(fix_elevation[sweep_index])
                )
            observable &= np.isfinite(first_gate) & np.isfinite(last_gate)
            observable &= ranges >= first_gate
            observable &= ranges <= last_gate
        if not np.any(observable):
            continue
        local_weight = np.ones_like(ranges, dtype=np.float64)
        if "range" in quality_weight_terms:
            local_weight *= np.exp(-np.square(ranges / float(range_decay_m)))
        if "beam" in quality_weight_terms:
            beam_radius = ranges * np.tan(np.deg2rad(representative_beam) * 0.5)
            local_weight *= 1.0 / (1.0 + np.square(beam_radius / float(beam_radius_ref_m)))
        if "vertical" in quality_weight_terms and fix_elevation.size > 1:
            local_weight *= 1.0 / (1.0 + np.square(gap_half / float(vertical_gap_ref_deg)))
        weight[ilevel] = np.where(observable, local_weight, 0.0)
    return weight


def _circle_lonlat(lon, lat, radius_m, npts=180):
    geod = pyproj.Geod(ellps="WGS84")
    az = np.linspace(0.0, 360.0, int(npts), endpoint=False)
    lons, lats, _ = geod.fwd(
        np.full_like(az, float(lon)),
        np.full_like(az, float(lat)),
        az,
        np.full_like(az, float(radius_m)),
    )
    return lons, lats


def _save_figure(fig, output_path, dpi=180, bbox_inches="tight"):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    save_kwargs = {"dpi": dpi}
    if bbox_inches is not None:
        save_kwargs["bbox_inches"] = bbox_inches
    fig.savefig(str(output_path), **save_kwargs)
    plt.close(fig)
    return str(output_path)


def _get_reference_font():
    available = {font.name for font in font_manager.fontManager.ttflist}
    family = next((name for name in _REFERENCE_FONT_CANDIDATES if name in available), "DejaVu Sans")
    return font_manager.FontProperties(family=family)


def _build_reference_cr_colormap():
    cmap = ListedColormap(_REFERENCE_CR_COLORS, name="reference_cr")
    cmap.set_bad((1.0, 1.0, 1.0, 0.0))
    norm = BoundaryNorm(_REFERENCE_CR_BOUNDS, cmap.N, clip=True)
    return cmap, norm


def _resolve_reference_canvas(plot_canvas_px=None):
    if plot_canvas_px is None:
        return _REFERENCE_CANVAS_PX
    width_px, height_px = plot_canvas_px
    return int(width_px), int(height_px)


def _validate_reference_plot_resolution(lon_res_deg, lat_res_deg):
    if float(lon_res_deg) > 0.015 or float(lat_res_deg) > 0.015:
        raise ValueError("reference plot style requires lon/lat resolution no coarser than 0.01 degree.")


def _format_reference_time(target_time):
    if target_time.tzinfo is None:
        target_time = target_time.replace(tzinfo=timezone.utc)
    if ZoneInfo is not None:
        local_time = target_time.astimezone(ZoneInfo("Asia/Shanghai"))
    else:  # pragma: no cover
        local_time = target_time.astimezone(timezone(timedelta(hours=8)))
    return local_time.strftime("%Y-%m-%d %H:%M")


def _draw_reference_colorbar(fig, cmap, norm, font_properties):
    nd_ax = fig.add_axes([0.045, 0.565, 0.02, 0.04])
    nd_ax.set_facecolor("white")
    for spine in nd_ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_edgecolor("black")
    nd_ax.set_xticks([])
    nd_ax.set_yticks([])
    nd_ax.text(
        1.45,
        0.5,
        "ND",
        transform=nd_ax.transAxes,
        ha="left",
        va="center",
        fontsize=10,
        fontproperties=font_properties,
    )

    cax = fig.add_axes([0.045, 0.145, 0.02, 0.42])
    cbar = ColorbarBase(
        cax,
        cmap=cmap,
        norm=norm,
        boundaries=_REFERENCE_CR_BOUNDS,
        ticks=_REFERENCE_CR_TICKS,
        orientation="vertical",
    )
    cbar.ax.invert_yaxis()
    cbar.ax.tick_params(length=0, labelsize=10)
    cbar.ax.set_yticklabels([str(int(value)) for value in _REFERENCE_CR_TICKS], fontproperties=font_properties)
    return cbar


def _plot_reference_cr_field(lon, lat, data, target_time, output_path, plot_canvas_px=None):
    from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

    from ..configure.location_config import get_cn_shp_info
    from ..draw._plot_core import ccrs, cfeature

    if ccrs is None:
        raise ImportError("reference plot style requires cartopy.")

    width_px, height_px = _resolve_reference_canvas(plot_canvas_px=plot_canvas_px)
    font_properties = _get_reference_font()
    cmap, norm = _build_reference_cr_colormap()
    dpi = 100
    fig = plt.figure(dpi=dpi, facecolor="white")
    fig.set_size_inches((width_px + 1e-3) / dpi, (height_px + 1e-3) / dpi, forward=True)
    ax = fig.add_axes([0.11, 0.08, 0.80, 0.84], projection=ccrs.PlateCarree())
    ax.set_facecolor("white")

    masked = np.ma.masked_invalid(np.asarray(data, dtype=np.float64))
    mesh = ax.pcolormesh(
        lon,
        lat,
        masked,
        cmap=cmap,
        norm=norm,
        shading="auto",
        transform=ccrs.PlateCarree(),
        zorder=2,
    )

    extent = [float(np.min(lon)), float(np.max(lon)), float(np.min(lat)), float(np.max(lat))]
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    shp_info = get_cn_shp_info()
    ax.add_feature(
        cfeature.ShapelyFeature(
            shp_info.geometries(),
            ccrs.PlateCarree(),
            edgecolor="#303030",
            facecolor="none",
        ),
        linewidth=1.0,
        zorder=4,
        alpha=0.9,
    )
    grid_x = np.arange(np.floor(extent[0]), np.floor(extent[1]) + 1.0, 1.0)
    grid_y = np.arange(np.floor(extent[2]), np.floor(extent[3]) + 1.0, 1.0)
    gridliner = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=False,
        linewidth=0.8,
        color="#8A8A8A",
        linestyle=(0, (5, 5)),
        alpha=0.9,
        xlocs=grid_x,
        ylocs=grid_y,
        zorder=1,
    )
    gridliner.xlines = True
    gridliner.ylines = True
    ax.set_xticks(grid_x, crs=ccrs.PlateCarree())
    ax.set_yticks(grid_y, crs=ccrs.PlateCarree())
    ax.xaxis.set_major_formatter(LongitudeFormatter(number_format=".0f"))
    ax.yaxis.set_major_formatter(LatitudeFormatter(number_format=".0f"))
    ax.tick_params(labelsize=10)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontproperties(font_properties)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
        spine.set_edgecolor("#303030")
    ax.set_xlabel("")
    ax.set_ylabel("")

    _draw_reference_colorbar(fig, cmap, norm, font_properties)
    fig.text(
        0.96,
        0.955,
        "产品名:组合反射率",
        ha="right",
        va="top",
        fontsize=17,
        fontproperties=font_properties,
        color="#222222",
        weight="bold",
    )
    fig.text(
        0.96,
        0.905,
        _format_reference_time(target_time),
        ha="right",
        va="top",
        fontsize=17,
        fontproperties=font_properties,
        color="#222222",
    )
    fig.text(
        0.045,
        0.095,
        "单位:dBZ",
        ha="left",
        va="bottom",
        fontsize=12,
        fontproperties=font_properties,
        color="#222222",
    )
    mesh.set_rasterized(True)
    return _save_figure(fig, output_path, dpi=dpi, bbox_inches=None)


def plot_network_overview(selected_radars, worker_results, lon_min, lon_max, lat_min, lat_max, output_path, field_name="dBZ"):
    fig, ax = plt.subplots(figsize=(8, 7))
    for radar_item, result in zip(selected_radars, worker_results):
        ax.scatter(result["radar_lon"], result["radar_lat"], s=32, c="black")
        ax.text(result["radar_lon"], result["radar_lat"], " %s" % result["radar_id"], fontsize=8, va="bottom")
        coverage_km = result["actual_max_range_km_by_field"].get(field_name)
        if coverage_km is not None:
            circle_lon, circle_lat = _circle_lonlat(result["radar_lon"], result["radar_lat"], coverage_km * 1000.0)
            ax.plot(circle_lon, circle_lat, color="#1f77b4", linewidth=0.8, alpha=0.8)
        first_gate_m = result.get("first_gate_m_by_field", {}).get(field_name)
        if first_gate_m is not None and first_gate_m > 0:
            blind_lon, blind_lat = _circle_lonlat(result["radar_lon"], result["radar_lat"], first_gate_m)
            ax.plot(blind_lon, blind_lat, color="#d62728", linewidth=0.8, linestyle="--", alpha=0.8)
    lon_min = float(lon_min)
    lon_max = float(lon_max)
    lat_min = float(lat_min)
    lat_max = float(lat_max)
    if lon_min == lon_max:
        lon_min -= 0.01
        lon_max += 0.01
    if lat_min == lat_max:
        lat_min -= 0.01
        lat_max += 0.01
    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title("Radar Coverage Overview")
    ax.grid(True, alpha=0.25, linestyle=":")
    return _save_figure(fig, output_path)


def _plot_2d_field(lon, lat, data, title, output_path, cmap="viridis", vmin=None, vmax=None):
    fig, ax = plt.subplots(figsize=(8, 6))
    mesh = ax.pcolormesh(lon, lat, np.asarray(data, dtype=np.float64), shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    fig.colorbar(mesh, ax=ax, shrink=0.9)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(title)
    ax.grid(True, alpha=0.2, linestyle=":")
    return _save_figure(fig, output_path)


def plot_network_level_diagnostic(lon, lat, level_height, cappi, coverage_count, blind_mask, output_path, field_title="CAPPI"):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    dbz_range = CINRAD_field_normvar.get("dBZ", (0.0, 75.0))
    mesh = axes[0].pcolormesh(lon, lat, np.asarray(cappi, dtype=np.float64), shading="auto",
                              cmap=CINRAD_COLORMAP.get("dBZ", "turbo"), vmin=dbz_range[0], vmax=dbz_range[1])
    fig.colorbar(mesh, ax=axes[0], shrink=0.9, label="dBZ")
    axes[0].set_title("%s %.0f m" % (field_title, level_height))
    blind_numeric = np.where(np.asarray(blind_mask, dtype=bool), 0.0, np.asarray(coverage_count, dtype=np.float64))
    mesh_cov = axes[1].pcolormesh(lon, lat, blind_numeric, shading="auto", cmap="viridis", vmin=0.0)
    fig.colorbar(mesh_cov, ax=axes[1], shrink=0.9, label="Radar count")
    if np.asarray(blind_mask).shape[0] >= 2 and np.asarray(blind_mask).shape[1] >= 2:
        axes[1].contour(lon, lat, np.asarray(blind_mask, dtype=np.float64), levels=[0.5], colors="red", linewidths=0.8)
    axes[1].set_title("Coverage / Blind Mask %.0f m" % level_height)
    for ax in axes:
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.grid(True, alpha=0.2, linestyle=":")
    return _save_figure(fig, output_path)


def grid_single_radar_to_latlon_3d(
    prd,
    field_name,
    grid_x,
    grid_y,
    level_heights,
    radar_x,
    radar_y,
    fillvalue=-999.0,
    effective_earth_radius=None,
    range_mode=None,
    max_range_km=None,
    blind_method="mask",
    return_metadata=False,
):
    """Grid a single radar volume onto the shared Cartesian grid."""
    if hasattr(prd, "_resolve_field_range_mode"):
        range_mode = prd._resolve_field_range_mode(field_name, range_mode=range_mode)
    prd.get_vol_data(
        field_name=field_name,
        fillvalue=fillvalue,
        range_mode=range_mode,
        max_range_km=max_range_km,
    )
    vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, _, _ = prd.vol
    local_x = np.asarray(grid_x, dtype=np.float64) - float(radar_x)
    local_y = np.asarray(grid_y, dtype=np.float64) - float(radar_y)
    grid_value = get_CAPPI_3d(
        vol_azimuth,
        vol_range,
        np.asarray(fix_elevation, dtype=np.float64),
        vol_value,
        float(radar_height),
        local_x,
        local_y,
        np.asarray(level_heights, dtype=np.float64),
        float(fillvalue),
        effective_earth_radius=effective_earth_radius,
        blind_method=blind_method,
    )
    valid_source = [values[np.isfinite(values) & (values != fillvalue)] for values in vol_value]
    valid_source = [values for values in valid_source if values.size]
    if not valid_source:
        return np.full_like(grid_value, float(fillvalue), dtype=np.float64)
    source_min = min(float(np.min(values)) for values in valid_source)
    source_max = max(float(np.max(values)) for values in valid_source)
    tolerance = max(1e-6, 1e-6 * max(abs(source_min), abs(source_max), 1.0))
    invalid_mask = np.isfinite(grid_value) & (
        (grid_value < (source_min - tolerance)) | (grid_value > (source_max + tolerance))
    )
    if np.any(invalid_mask):
        grid_value = np.where(invalid_mask, float(fillvalue), grid_value)
    if return_metadata:
        actual_max_range_km = max(float(ranges[-1]) for ranges in vol_range) / 1000.0
        actual_max_range_km_by_sweep = [float(ranges[-1]) / 1000.0 for ranges in vol_range]
        return grid_value, {
            "range_mode": range_mode,
            "actual_max_range_km": actual_max_range_km,
            "actual_max_range_km_by_sweep": actual_max_range_km_by_sweep,
        }
    return grid_value


def grid_single_radar_cr_to_latlon(
    prd,
    field_name,
    grid_x,
    grid_y,
    radar_x,
    radar_y,
    fillvalue=-999.0,
    effective_earth_radius=None,
    range_mode=None,
    max_range_km=None,
):
    """Grid a single radar composite reflectivity onto the shared Cartesian grid."""
    if hasattr(prd, "_resolve_field_range_mode"):
        range_mode = prd._resolve_field_range_mode(field_name, range_mode=range_mode)
    prd.get_vol_data(
        field_name=field_name,
        fillvalue=fillvalue,
        range_mode=range_mode,
        max_range_km=max_range_km,
    )
    vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, _, _ = prd.vol
    local_x = np.asarray(grid_x, dtype=np.float64) - float(radar_x)
    local_y = np.asarray(grid_y, dtype=np.float64) - float(radar_y)
    from ..core.RadarGrid import get_CR_xy as py_get_cr_xy

    try:
        from ..core.RadarGridC import get_CR_xy as cy_get_cr_xy
        get_cr_xy = cy_get_cr_xy
    except ImportError:
        get_cr_xy = py_get_cr_xy
    grid_value = get_cr_xy(
        vol_azimuth,
        vol_range,
        np.asarray(fix_elevation, dtype=np.float64),
        vol_value,
        float(radar_height),
        local_x,
        local_y,
        float(fillvalue),
        effective_earth_radius=effective_earth_radius,
    )
    return np.asarray(grid_value, dtype=np.float64)


def _prepare_radar_grid_worker(
    radar_item,
    station_lon,
    station_lat,
    station_alt,
    field_names,
    field_range_mode,
    grid_x,
    grid_y,
    level_heights,
    product_level_heights,
    lon_0,
    lat_0,
    fillvalue,
    effective_earth_radius,
    max_range_km,
    blind_method,
    use_qc,
    qc_method,
    qc_band,
    qc_use_existing_kdp,
    qc_fallback,
    qc_clear_air_mode,
    qc_clear_air_max_ref,
    qc_clear_air_max_rhohv,
    qc_clear_air_max_phidp_texture,
    qc_clear_air_max_snr,
    quality_weight_terms,
    range_decay_m,
    beam_radius_ref_m,
    vertical_gap_ref_deg,
    qc_weight_field,
    blockage_config,
):
    """Read and grid a single radar volume for all requested fields."""
    prd = read_auto(
        radar_item["path"],
        station_lon=station_lon,
        station_lat=station_lat,
        station_alt=station_alt,
        effective_earth_radius=effective_earth_radius,
    )
    radar_lon = float(station_lon)
    radar_lat = float(station_lat)
    proj = pyproj.Proj({"proj": "aeqd", "lon_0": lon_0, "lat_0": lat_0})
    radar_x, radar_y = proj(radar_lon, radar_lat, inverse=False)
    qc_applied = False
    qc_error = None
    if use_qc:
        if qc_method != "dualpol":
            raise ValueError("Unsupported qc_method: %s" % qc_method)
        if _radar_supports_dualpol_qc(prd):
            try:
                prd = prd.apply_dualpol_qc(
                    inplace=True,
                    band=qc_band,
                    use_existing_kdp=qc_use_existing_kdp,
                    clear_air_mode=qc_clear_air_mode,
                    clear_air_max_ref=qc_clear_air_max_ref,
                    clear_air_max_rhohv=qc_clear_air_max_rhohv,
                    clear_air_max_phidp_texture=qc_clear_air_max_phidp_texture,
                    clear_air_max_snr=qc_clear_air_max_snr,
                )
                qc_applied = True
            except Exception as exc:
                if qc_fallback != "original":
                    raise
                qc_error = str(exc)
        else:
            if qc_fallback != "original":
                raise ValueError("Radar %s does not have the moments required for dual-pol QC." % radar_item["radar_id"])
            qc_error = "missing_dualpol_fields"

    field_volumes = {}
    actual_max_range_km_by_field = {}
    actual_max_range_km_by_sweep_by_field = {}
    source_field_by_field = {}
    first_gate_m_by_field = {}
    field_quality_weights = {}
    cr_grid = None
    observability_mask = None
    qc_weight_grid = None
    blockage_weight = None
    product_reflectivity_grid = None
    product_reflectivity_quality = None
    product_qc_weight_grid = None
    reflectivity_source_field = _resolve_qc_field_name("dBZ", qc_applied)
    for field_name in field_names:
        source_field_name = _resolve_qc_field_name(field_name, qc_applied)
        try:
            radar_grid, metadata = grid_single_radar_to_latlon_3d(
                prd,
                source_field_name,
                grid_x,
                grid_y,
                level_heights,
                radar_x,
                radar_y,
                fillvalue=fillvalue,
                effective_earth_radius=getattr(prd, "effective_earth_radius", effective_earth_radius),
                range_mode=field_range_mode[field_name],
                max_range_km=max_range_km,
                blind_method=blind_method,
                return_metadata=True,
            )
            actual_max_range_km_by_field[field_name] = metadata["actual_max_range_km"]
            actual_max_range_km_by_sweep_by_field[field_name] = metadata.get("actual_max_range_km_by_sweep")
            prd.get_vol_data(
                field_name=source_field_name,
                fillvalue=fillvalue,
                range_mode=field_range_mode[field_name],
                max_range_km=max_range_km,
            )
            _, vol_range, fix_elevation, _, radar_height, _, _ = prd.vol
            beam_widths = _get_sweep_beam_widths(prd, len(fix_elevation))
            first_gate_m_by_field[field_name] = min(float(ranges[0]) for ranges in vol_range)
            geometry_quality = _compute_quality_weight_volume(
                vol_range,
                fix_elevation,
                radar_height,
                grid_x,
                grid_y,
                level_heights,
                radar_x,
                radar_y,
                beam_widths=beam_widths,
                effective_earth_radius=getattr(prd, "effective_earth_radius", effective_earth_radius),
                quality_weight_terms=quality_weight_terms,
                range_decay_m=range_decay_m,
                beam_radius_ref_m=beam_radius_ref_m,
                vertical_gap_ref_deg=vertical_gap_ref_deg,
            )
            if source_field_name == reflectivity_source_field:
                cr_grid = grid_single_radar_cr_to_latlon(
                    prd,
                    source_field_name,
                    grid_x,
                    grid_y,
                    radar_x,
                    radar_y,
                    fillvalue=fillvalue,
                    effective_earth_radius=getattr(prd, "effective_earth_radius", effective_earth_radius),
                    range_mode=field_range_mode[field_name],
                    max_range_km=max_range_km,
                )
                observability_mask = _compute_observability_mask(
                    vol_range,
                    fix_elevation,
                    radar_height,
                    grid_x,
                    grid_y,
                    level_heights,
                    radar_x,
                    radar_y,
                    beam_widths=beam_widths,
                    effective_earth_radius=getattr(prd, "effective_earth_radius", effective_earth_radius),
                )
                if qc_applied and ("qc" in quality_weight_terms) and qc_weight_field:
                    try:
                        qc_grid = grid_single_radar_to_latlon_3d(
                            prd,
                            qc_weight_field,
                            grid_x,
                            grid_y,
                            level_heights,
                            radar_x,
                            radar_y,
                            fillvalue=fillvalue,
                            effective_earth_radius=getattr(prd, "effective_earth_radius", effective_earth_radius),
                            range_mode="aligned",
                            max_range_km=max_range_km,
                            blind_method="mask",
                            return_metadata=False,
                        )
                        qc_weight_grid = np.clip(
                            np.where(
                                np.isfinite(qc_grid) & (qc_grid != fillvalue),
                                qc_grid,
                                0.0,
                            ),
                            0.0,
                            1.0,
                        )
                    except KeyError:
                        qc_weight_grid = None
                blockage_value = 1.0
                if blockage_config:
                    blockage_value = blockage_config.get(radar_item["radar_id"], 1.0)
                blockage_array = np.asarray(blockage_value, dtype=np.float64)
                if blockage_array.ndim == 0:
                    blockage_weight = np.full(grid_x.shape, float(blockage_array), dtype=np.float64)
                else:
                    if blockage_array.shape != np.asarray(grid_x).shape:
                        raise ValueError(
                            "blockage_config for %s must match the target grid shape." % radar_item["radar_id"]
                        )
                    blockage_weight = np.clip(blockage_array, 0.0, 1.0)
                if product_level_heights is not None and len(product_level_heights) > 0:
                    product_level_heights = np.asarray(product_level_heights, dtype=np.float64)
                    if not np.array_equal(product_level_heights, np.asarray(level_heights, dtype=np.float64)):
                        product_reflectivity_grid = grid_single_radar_to_latlon_3d(
                            prd,
                            source_field_name,
                            grid_x,
                            grid_y,
                            product_level_heights,
                            radar_x,
                            radar_y,
                            fillvalue=fillvalue,
                            effective_earth_radius=getattr(prd, "effective_earth_radius", effective_earth_radius),
                            range_mode=field_range_mode[field_name],
                            max_range_km=max_range_km,
                            blind_method=blind_method,
                            return_metadata=False,
                        )
                        product_reflectivity_quality = _compute_quality_weight_volume(
                            vol_range,
                            fix_elevation,
                            radar_height,
                            grid_x,
                            grid_y,
                            product_level_heights,
                            radar_x,
                            radar_y,
                            beam_widths=beam_widths,
                            effective_earth_radius=getattr(prd, "effective_earth_radius", effective_earth_radius),
                            quality_weight_terms=quality_weight_terms,
                            range_decay_m=range_decay_m,
                            beam_radius_ref_m=beam_radius_ref_m,
                            vertical_gap_ref_deg=vertical_gap_ref_deg,
                        )
                        if qc_applied and ("qc" in quality_weight_terms) and qc_weight_field:
                            if product_qc_weight_grid is None:
                                try:
                                    product_qc_grid = grid_single_radar_to_latlon_3d(
                                        prd,
                                        qc_weight_field,
                                        grid_x,
                                        grid_y,
                                        product_level_heights,
                                        radar_x,
                                        radar_y,
                                        fillvalue=fillvalue,
                                        effective_earth_radius=getattr(prd, "effective_earth_radius", effective_earth_radius),
                                        range_mode="aligned",
                                        max_range_km=max_range_km,
                                        blind_method="mask",
                                        return_metadata=False,
                                    )
                                    product_qc_weight_grid = np.clip(
                                        np.where(
                                            np.isfinite(product_qc_grid) & (product_qc_grid != fillvalue),
                                            product_qc_grid,
                                            0.0,
                                        ),
                                        0.0,
                                        1.0,
                                    )
                                except KeyError:
                                    product_qc_weight_grid = None
                        if product_qc_weight_grid is not None and ("qc" in quality_weight_terms):
                            product_reflectivity_quality = product_reflectivity_quality * product_qc_weight_grid
                        if blockage_weight is not None and ("blockage" in quality_weight_terms):
                            product_reflectivity_quality = (
                                product_reflectivity_quality * blockage_weight[np.newaxis, :, :]
                            )
            quality_weight = geometry_quality
            if qc_weight_grid is not None and ("qc" in quality_weight_terms):
                quality_weight = quality_weight * qc_weight_grid
            if blockage_weight is not None and ("blockage" in quality_weight_terms):
                quality_weight = quality_weight * blockage_weight[np.newaxis, :, :]
            field_quality_weights[field_name] = np.asarray(quality_weight, dtype=np.float64)
        except KeyError:
            radar_grid = np.full(
                (len(level_heights),) + np.asarray(grid_x).shape,
                float(fillvalue),
                dtype=np.float64,
            )
            actual_max_range_km_by_field[field_name] = None
            first_gate_m_by_field[field_name] = None
            field_quality_weights[field_name] = np.zeros_like(radar_grid, dtype=np.float64)
            actual_max_range_km_by_sweep_by_field[field_name] = None
        field_volumes[field_name] = radar_grid
        source_field_by_field[field_name] = source_field_name

    return {
        "radar_id": radar_item["radar_id"],
        "path": radar_item["path"],
        "scan_time": radar_item["scan_time"].isoformat(),
        "radar_lon": radar_lon,
        "radar_lat": radar_lat,
        "radar_x": float(radar_x),
        "radar_y": float(radar_y),
        "field_volumes": field_volumes,
        "actual_max_range_km_by_field": actual_max_range_km_by_field,
        "actual_max_range_km_by_sweep_by_field": actual_max_range_km_by_sweep_by_field,
        "first_gate_m_by_field": first_gate_m_by_field,
        "source_field_by_field": source_field_by_field,
        "field_quality_weights": field_quality_weights,
        "product_reflectivity_grid": product_reflectivity_grid,
        "product_reflectivity_quality": product_reflectivity_quality,
        "cr_grid": cr_grid,
        "observability_mask": observability_mask,
        "qc_applied": qc_applied,
        "qc_error": qc_error,
    }


def _build_composite_weight_cache(radar_sites_xy, grid_x, grid_y, method="exp_weighted", influence_radius_m=300000.0):
    """Precompute per-radar 2D distance or weight fields reused across multiple variables."""
    method = str(method)
    radar_sites_xy = np.asarray(radar_sites_xy, dtype=np.float64)
    grid_x = np.asarray(grid_x, dtype=np.float64)
    grid_y = np.asarray(grid_y, dtype=np.float64)
    cache = []
    for radar_x, radar_y in radar_sites_xy:
        distance2 = (grid_x - radar_x) ** 2 + (grid_y - radar_y) ** 2
        if method in ("exp_weighted", "quality_weighted"):
            cache.append(np.exp(-4.0 * distance2 / (float(influence_radius_m) ** 2)))
        elif method == "nearest":
            cache.append(distance2)
        elif method == "max":
            return None
        else:
            raise ValueError("Unsupported composite_method: %s." % method)
    return cache


def compose_network_volume(
    radar_volumes,
    radar_sites_xy,
    grid_x,
    grid_y,
    fillvalue=-999.0,
    method="exp_weighted",
    influence_radius_m=300000.0,
    weight_cache=None,
    quality_weights=None,
):
    """
    Composite multiple single-radar 3D volumes onto the shared grid.

    Supported methods map to the overlap strategies summarized in
    "吴翀;双偏振雷达的资料质量分析,相态识別及组网应用[D];南京信息工程大学;2018年"
    (Chapter 4): nearest, max, exponential weighting, and a pycwr
    quality-weighted extension. The quality-weighted branch is also aligned in
    spirit with MRMS/RQI-style quality-aware radar merging, but it is not a
    verbatim reimplementation of one published national mosaic algorithm.
    """
    if not radar_volumes:
        raise ValueError("radar_volumes must contain at least one volume.")
    method = str(method)
    if method not in ("exp_weighted", "nearest", "max", "quality_weighted"):
        raise ValueError("Unsupported composite_method: %s." % method)
    influence_radius_m = float(influence_radius_m)
    if influence_radius_m <= 0:
        raise ValueError("influence_radius_m must be positive.")

    radar_sites_xy = np.asarray(radar_sites_xy, dtype=np.float64)
    if radar_sites_xy.shape != (len(radar_volumes), 2):
        raise ValueError("radar_sites_xy must have shape (nradar, 2).")

    grid_x = np.asarray(grid_x, dtype=np.float64)
    grid_y = np.asarray(grid_y, dtype=np.float64)
    radar_volumes = [np.asarray(volume, dtype=np.float64) for volume in radar_volumes]
    composite = np.full(radar_volumes[0].shape, float(fillvalue), dtype=np.float64)
    for volume in radar_volumes[1:]:
        if volume.shape != composite.shape:
            raise ValueError("All radar volumes must share the same shape.")

    valid_stack = np.stack(
        [np.isfinite(volume) & (volume != fillvalue) for volume in radar_volumes],
        axis=0,
    )
    valid_count = np.sum(valid_stack, axis=0, dtype=np.int16)
    volume_stack = np.stack(radar_volumes, axis=0)

    if method == "max":
        masked = np.where(valid_stack, volume_stack, -np.inf)
        reduced = np.max(masked, axis=0)
        has_data = valid_count > 0
        composite[has_data] = reduced[has_data]
        return composite, valid_count

    if weight_cache is None:
        weight_cache = _build_composite_weight_cache(
            radar_sites_xy,
            grid_x,
            grid_y,
            method=method,
            influence_radius_m=influence_radius_m,
        )
    weight_stack_2d = np.stack(weight_cache, axis=0)

    if method == "nearest":
        distance_stack = np.where(valid_stack, weight_stack_2d[:, np.newaxis, :, :], np.inf)
        nearest_index = np.argmin(distance_stack, axis=0)
        nearest_distance = np.take_along_axis(distance_stack, nearest_index[np.newaxis, :, :, :], axis=0)[0]
        chosen_value = np.take_along_axis(volume_stack, nearest_index[np.newaxis, :, :, :], axis=0)[0]
        has_data = np.isfinite(nearest_distance)
        composite[has_data] = chosen_value[has_data]
        return composite, valid_count

    weight_stack = weight_stack_2d[:, np.newaxis, :, :]
    if method == "quality_weighted":
        if quality_weights is None:
            quality_stack = np.ones_like(volume_stack, dtype=np.float64)
        else:
            quality_stack = np.stack([np.asarray(weight, dtype=np.float64) for weight in quality_weights], axis=0)
            if quality_stack.shape != volume_stack.shape:
                raise ValueError("quality_weights must match radar_volumes in shape.")
        weight_stack = weight_stack * quality_stack
    weighted_sum = np.sum(np.where(valid_stack, volume_stack * weight_stack, 0.0), axis=0)
    total_weight = np.sum(np.where(valid_stack, weight_stack, 0.0), axis=0)
    has_data = total_weight > 0
    composite[has_data] = weighted_sum[has_data] / total_weight[has_data]
    return composite, valid_count


def build_network_dataset(
    data_vars,
    lon,
    lat,
    level_heights,
    selected_radars,
    target_time,
    composite_method,
    influence_radius_m,
    fillvalue,
    effective_earth_radius,
    projection_origin,
    field_range_mode,
    blind_method,
    requested_max_range_km,
    actual_max_range_km_by_field,
    actual_max_range_km_by_sweep_by_field,
    use_qc,
    qc_method,
    qc_clear_air_mode,
    qc_clear_air_max_ref,
    qc_clear_air_max_rhohv,
    qc_clear_air_max_phidp_texture,
    qc_clear_air_max_snr,
    qc_applied_by_radar,
    source_field_by_radar,
    quality_weight_terms,
    config_path=None,
):
    """Create the xarray dataset for the 3D radar network product."""
    dataset = xr.Dataset(coords={"z": level_heights, "lat": lat, "lon": lon})
    corrected_fields = {"dBZ": "Zc", "ZDR": "ZDRc", "KDP": "KDPc"}
    for field_name, values in data_vars.items():
        radar_sources = source_field_by_radar[field_name]
        field_qc_flags = []
        for radar_id, source_field_name in radar_sources.items():
            expected_corrected = corrected_fields.get(field_name)
            field_qc_flags.append(
                bool(qc_applied_by_radar.get(radar_id, {}).get("applied")) and expected_corrected == source_field_name
            )
        if field_qc_flags and all(field_qc_flags):
            qc_attr = "true"
        elif any(field_qc_flags):
            qc_attr = "partial"
        else:
            qc_attr = "false"
        dataset[field_name] = (("z", "lat", "lon"), np.where(values == fillvalue, np.nan, values))
        dataset[field_name].attrs = {
            "standard_name": "radar_volume_mosaic",
            "long_name": "%s 3D radar network mosaic" % field_name,
            "coordinates": "z lat lon",
            "range_mode": field_range_mode[field_name],
            "requested_max_range_km": (
                "full"
                if requested_max_range_km is None
                else float(requested_max_range_km)
            ),
            "actual_max_range_km_by_radar": json.dumps(actual_max_range_km_by_field[field_name], ensure_ascii=False),
            "actual_max_range_km_by_radar_and_sweep": json.dumps(
                actual_max_range_km_by_sweep_by_field[field_name],
                ensure_ascii=False,
            ),
            "source_field_name_by_radar": json.dumps(source_field_by_radar[field_name], ensure_ascii=False),
            "qc_applied": qc_attr,
        }
    dataset.coords["z"].attrs = {"units": "m", "long_name": "height_above_mean_sea_level"}
    dataset.coords["lat"].attrs = {"units": "degree_north", "long_name": "latitude"}
    dataset.coords["lon"].attrs = {"units": "degree_east", "long_name": "longitude"}
    dataset.attrs = {
        "product_type": "radar_network_3d",
        "target_time": target_time.isoformat(),
        "radar_count": len(selected_radars),
        "radar_files": "|".join(item["path"] for item in selected_radars),
        "composite_method": composite_method,
        "influence_radius_m": float(influence_radius_m),
        "fillvalue": float(fillvalue),
        "projection": "aeqd",
        "projection_lon_0": float(projection_origin[0]),
        "projection_lat_0": float(projection_origin[1]),
        "field_range_mode": json.dumps(field_range_mode, ensure_ascii=False),
        "blind_method": blind_method,
        "requested_max_range_km": "full" if requested_max_range_km is None else float(requested_max_range_km),
        "actual_max_range_km_by_field": json.dumps(actual_max_range_km_by_field, ensure_ascii=False),
        "use_qc": "true" if use_qc else "false",
        "qc_method": qc_method if use_qc else "none",
        "qc_clear_air_mode": qc_clear_air_mode if use_qc else "ignore",
        "qc_clear_air_max_ref": float(qc_clear_air_max_ref),
        "qc_clear_air_max_rhohv": float(qc_clear_air_max_rhohv),
        "qc_clear_air_max_phidp_texture": float(qc_clear_air_max_phidp_texture),
        "qc_clear_air_max_snr": float(qc_clear_air_max_snr),
        "qc_applied_by_radar": json.dumps(qc_applied_by_radar, ensure_ascii=False),
        "quality_weight_terms": json.dumps(quality_weight_terms),
        "algorithm_references": json.dumps(_PRODUCT_REFERENCE_NOTES["MOSAIC"], ensure_ascii=False),
        "effective_earth_radius": (
            "default"
            if effective_earth_radius is None
            else float(effective_earth_radius)
        ),
        "network_config_path": "" if config_path is None else str(config_path),
    }
    return dataset


def radar_network_3d_to_netcdf(dataset, output_path):
    """Write the 3D radar network dataset to a netCDF file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    dataset.to_netcdf(str(output_path))
    return str(output_path)


def _default_network_output_path(output_dir, target_time, field_names):
    """Build a default netCDF output path under the configured directory."""
    stem = "radar_network_3d_%s_%s.nc" % (
        target_time.strftime("%Y%m%d%H%M%S"),
        "_".join(field_names),
    )
    return str(Path(output_dir) / stem)


def run_radar_network_3d(
    target_time,
    config_path=None,
    radar_dirs=None,
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    lon_res_deg=None,
    lat_res_deg=None,
    level_heights=None,
    product_level_heights=None,
    field_names=None,
    output_products=None,
    output_path=None,
    time_tolerance_minutes=None,
    composite_method=None,
    influence_radius_m=None,
    fillvalue=None,
    effective_earth_radius=None,
    field_range_mode=None,
    max_range_km=None,
    parallel=None,
    max_workers=None,
    blind_method=None,
    use_qc=None,
    qc_method=None,
    qc_band=None,
    qc_use_existing_kdp=None,
    qc_fallback=None,
    qc_clear_air_mode=None,
    qc_clear_air_max_ref=None,
    qc_clear_air_max_rhohv=None,
    qc_clear_air_max_phidp_texture=None,
    qc_clear_air_max_snr=None,
    plot_overview=None,
    plot_style=None,
    plot_canvas_px=None,
    plot_height_levels=None,
    plot_output_dir=None,
    plot_format=None,
    plot_product_for_levels=None,
    coverage_policy=None,
    quality_weight_terms=None,
    range_decay_m=None,
    beam_radius_ref_m=None,
    vertical_gap_ref_deg=None,
    qc_weight_field=None,
    blockage_config=None,
    vil_min_dbz=None,
    vil_max_dbz_cap=None,
    et_threshold_dbz=None,
    pattern=None,
    lon_0=None,
    lat_0=None,
):
    """Run the full 3D radar network workflow and optionally write netCDF output."""
    if isinstance(target_time, str):
        target_time = datetime.fromisoformat(target_time)
    if target_time is None:
        raise ValueError("target_time is required.")

    network_cfg = load_network_config(config_path) if config_path is not None else getattr(cfg, "network", None)
    radar_dirs = _resolve_network_value("radar_dirs", radar_dirs, network_cfg=network_cfg)
    lon_min = _resolve_network_value("lon_min", lon_min, network_cfg=network_cfg)
    lon_max = _resolve_network_value("lon_max", lon_max, network_cfg=network_cfg)
    lat_min = _resolve_network_value("lat_min", lat_min, network_cfg=network_cfg)
    lat_max = _resolve_network_value("lat_max", lat_max, network_cfg=network_cfg)
    lon_res_deg = _resolve_network_value("lon_res_deg", lon_res_deg, network_cfg=network_cfg)
    lat_res_deg = _resolve_network_value("lat_res_deg", lat_res_deg, network_cfg=network_cfg)
    level_heights = _resolve_network_value("level_heights", level_heights, network_cfg=network_cfg)
    time_tolerance_minutes = _resolve_network_value("time_tolerance_minutes", time_tolerance_minutes, network_cfg=network_cfg)
    composite_method = _resolve_network_value("composite_method", composite_method, network_cfg=network_cfg) or "quality_weighted"
    influence_radius_m = _resolve_network_value("influence_radius_m", influence_radius_m, network_cfg=network_cfg) or 300000.0
    fillvalue = _resolve_network_value("fillvalue", fillvalue, network_cfg=network_cfg)
    effective_earth_radius = _resolve_network_value("effective_earth_radius", effective_earth_radius, network_cfg=network_cfg)
    max_range_km = _resolve_network_value("max_range_km", max_range_km, network_cfg=network_cfg)
    pattern = _resolve_network_value("file_pattern", pattern, network_cfg=network_cfg) or "*.bin*"
    parallel = _resolve_network_value("parallel", parallel, network_cfg=network_cfg)
    max_workers = _resolve_network_value("max_workers", max_workers, network_cfg=network_cfg)
    blind_method = _resolve_network_value("blind_method", blind_method, network_cfg=network_cfg)
    use_qc = _resolve_network_value("use_qc", use_qc, network_cfg=network_cfg)
    qc_method = _resolve_network_value("qc_method", qc_method, network_cfg=network_cfg) or "dualpol"
    qc_band = _resolve_network_value("qc_band", qc_band, network_cfg=network_cfg) or "C"
    qc_use_existing_kdp = _resolve_network_value("qc_use_existing_kdp", qc_use_existing_kdp, network_cfg=network_cfg)
    qc_fallback = _resolve_network_value("qc_fallback", qc_fallback, network_cfg=network_cfg) or "original"
    qc_clear_air_mode = _resolve_network_value("qc_clear_air_mode", qc_clear_air_mode, network_cfg=network_cfg) or "label"
    qc_clear_air_max_ref = _resolve_network_value("qc_clear_air_max_ref", qc_clear_air_max_ref, network_cfg=network_cfg)
    qc_clear_air_max_rhohv = _resolve_network_value("qc_clear_air_max_rhohv", qc_clear_air_max_rhohv, network_cfg=network_cfg)
    qc_clear_air_max_phidp_texture = _resolve_network_value(
        "qc_clear_air_max_phidp_texture",
        qc_clear_air_max_phidp_texture,
        network_cfg=network_cfg,
    )
    qc_clear_air_max_snr = _resolve_network_value("qc_clear_air_max_snr", qc_clear_air_max_snr, network_cfg=network_cfg)
    output_products = _coerce_output_products(output_products, network_cfg=network_cfg)
    plot_overview = _resolve_network_value("plot_overview", plot_overview, network_cfg=network_cfg)
    plot_style = _resolve_network_value("plot_style", plot_style, network_cfg=network_cfg) or "reference"
    plot_canvas_px = _resolve_network_value("plot_canvas_px", plot_canvas_px, network_cfg=network_cfg)
    plot_height_levels = _coerce_plot_height_levels(level_heights, plot_height_levels, network_cfg=network_cfg)
    product_level_heights = _coerce_product_level_heights(level_heights, product_level_heights, network_cfg=network_cfg)
    plot_output_dir = _resolve_network_value("plot_output_dir", plot_output_dir, network_cfg=network_cfg)
    plot_format = _resolve_network_value("plot_format", plot_format, network_cfg=network_cfg) or "png"
    plot_product_for_levels = _resolve_network_value("plot_product_for_levels", plot_product_for_levels, network_cfg=network_cfg) or "CAPPI"
    coverage_policy = _resolve_network_value("coverage_policy", coverage_policy, network_cfg=network_cfg) or "observable_only"
    quality_weight_terms = _coerce_quality_weight_terms(quality_weight_terms, network_cfg=network_cfg)
    range_decay_m = _resolve_network_value("range_decay_m", range_decay_m, network_cfg=network_cfg)
    beam_radius_ref_m = _resolve_network_value("beam_radius_ref_m", beam_radius_ref_m, network_cfg=network_cfg)
    vertical_gap_ref_deg = _resolve_network_value("vertical_gap_ref_deg", vertical_gap_ref_deg, network_cfg=network_cfg)
    qc_weight_field = _resolve_network_value("qc_weight_field", qc_weight_field, network_cfg=network_cfg)
    blockage_config = _resolve_network_value("blockage_config", blockage_config, network_cfg=network_cfg)
    vil_min_dbz = _resolve_network_value("vil_min_dbz", vil_min_dbz, network_cfg=network_cfg)
    vil_max_dbz_cap = _resolve_network_value("vil_max_dbz_cap", vil_max_dbz_cap, network_cfg=network_cfg)
    et_threshold_dbz = _resolve_network_value("et_threshold_dbz", et_threshold_dbz, network_cfg=network_cfg)
    field_names = _coerce_field_names(field_names, network_cfg=network_cfg)
    _validate_network_field_names(field_names)
    field_range_mode = _coerce_field_range_mode(field_names, field_range_mode, network_cfg=network_cfg)
    output_dir = _resolve_network_value("output_dir", None, network_cfg=network_cfg)

    if radar_dirs is None:
        raise ValueError("radar_dirs are required.")
    if level_heights is None:
        raise ValueError("level_heights are required.")
    if None in (lon_min, lon_max, lat_min, lat_max, lon_res_deg, lat_res_deg):
        raise ValueError("lon/lat bounds and resolution are required.")
    if fillvalue is None:
        fillvalue = -999.0
    if qc_use_existing_kdp is None:
        qc_use_existing_kdp = True
    if qc_clear_air_max_ref is None:
        qc_clear_air_max_ref = 15.0
    if qc_clear_air_max_rhohv is None:
        qc_clear_air_max_rhohv = 0.97
    if qc_clear_air_max_phidp_texture is None:
        qc_clear_air_max_phidp_texture = 10.0
    if qc_clear_air_max_snr is None:
        qc_clear_air_max_snr = 20.0
    if use_qc is None:
        use_qc = False
    if plot_overview is None:
        plot_overview = True
    if coverage_policy != "observable_only":
        raise ValueError("Only coverage_policy='observable_only' is currently supported.")
    if vil_min_dbz is None:
        vil_min_dbz = 18.0
    if vil_max_dbz_cap is None:
        vil_max_dbz_cap = 56.0
    if et_threshold_dbz is None:
        et_threshold_dbz = 18.0
    blind_method = _normalize_blind_method(blind_method if blind_method is not None else "hybrid")
    plot_style = str(plot_style).lower()
    if plot_style not in {"reference", "simple"}:
        raise ValueError("plot_style must be 'reference' or 'simple'.")
    if ("CR" in output_products or "VIL" in output_products or "ET" in output_products or plot_height_levels.size > 0) and "dBZ" not in field_names:
        field_names = list(field_names) + ["dBZ"]
        field_range_mode = _coerce_field_range_mode(field_names, field_range_mode, network_cfg=network_cfg)
    if plot_style == "reference" and ("CR" in output_products):
        _validate_reference_plot_resolution(lon_res_deg, lat_res_deg)

    if range_decay_m is None:
        range_decay_m = 230000.0
    if beam_radius_ref_m is None:
        beam_radius_ref_m = 3000.0
    if vertical_gap_ref_deg is None:
        vertical_gap_ref_deg = 0.5
    if qc_weight_field is None:
        qc_weight_field = "QC_MASK"

    selected_radars = select_radar_files(
        radar_dirs,
        target_time,
        tolerance_minutes=time_tolerance_minutes if time_tolerance_minutes is not None else 10,
        pattern=pattern,
    )
    station_info = []
    for item in selected_radars:
        try:
            station_lat, station_lon, station_alt, _ = get_radar_info(item["path"])
            station_info.append((float(station_lat), float(station_lon), float(station_alt)))
        except Exception:
            prd = read_auto(
                item["path"],
                effective_earth_radius=effective_earth_radius,
            )
            station_info.append(_get_prd_station_info(prd))
    radar_lats = np.array([float(info[0]) for info in station_info], dtype=np.float64)
    radar_lons = np.array([float(info[1]) for info in station_info], dtype=np.float64)
    radar_alts = np.array([float(info[2]) for info in station_info], dtype=np.float64)
    if lon_0 is None:
        lon_0 = float(np.mean(radar_lons))
    if lat_0 is None:
        lat_0 = float(np.mean(radar_lats))

    lon, lat, grid_lon, grid_lat = build_latlon_grid(
        lon_min,
        lon_max,
        lat_min,
        lat_max,
        lon_res_deg,
        lat_res_deg,
    )
    proj = pyproj.Proj({"proj": "aeqd", "lon_0": lon_0, "lat_0": lat_0})
    grid_x, grid_y = proj(grid_lon, grid_lat, inverse=False)
    level_heights = np.asarray(level_heights, dtype=np.float64)

    worker_results = []
    use_parallel = bool(parallel) and len(selected_radars) > 1
    if use_parallel:
        if max_workers is None:
            max_workers = min(len(selected_radars), os.cpu_count() or 1)
        mp_context = mp.get_context("fork") if hasattr(mp, "get_context") and os.name == "posix" else None
        executor_kwargs = {"max_workers": int(max_workers)}
        if mp_context is not None:
            executor_kwargs["mp_context"] = mp_context
        with ProcessPoolExecutor(**executor_kwargs) as executor:
            worker_results = list(
                executor.map(
                    _prepare_radar_grid_worker,
                    selected_radars,
                    radar_lons,
                    radar_lats,
                    radar_alts,
                    [field_names] * len(selected_radars),
                    [field_range_mode] * len(selected_radars),
                    [grid_x] * len(selected_radars),
                    [grid_y] * len(selected_radars),
                    [level_heights] * len(selected_radars),
                    [product_level_heights] * len(selected_radars),
                    [lon_0] * len(selected_radars),
                    [lat_0] * len(selected_radars),
                    [fillvalue] * len(selected_radars),
                    [effective_earth_radius] * len(selected_radars),
                    [max_range_km] * len(selected_radars),
                    [blind_method] * len(selected_radars),
                    [use_qc] * len(selected_radars),
                    [qc_method] * len(selected_radars),
                    [qc_band] * len(selected_radars),
                    [qc_use_existing_kdp] * len(selected_radars),
                    [qc_fallback] * len(selected_radars),
                    [qc_clear_air_mode] * len(selected_radars),
                    [qc_clear_air_max_ref] * len(selected_radars),
                    [qc_clear_air_max_rhohv] * len(selected_radars),
                    [qc_clear_air_max_phidp_texture] * len(selected_radars),
                    [qc_clear_air_max_snr] * len(selected_radars),
                    [quality_weight_terms] * len(selected_radars),
                    [range_decay_m] * len(selected_radars),
                    [beam_radius_ref_m] * len(selected_radars),
                    [vertical_gap_ref_deg] * len(selected_radars),
                    [qc_weight_field] * len(selected_radars),
                    [blockage_config] * len(selected_radars),
                )
            )
    else:
        worker_results = [
            _prepare_radar_grid_worker(
                radar_item,
                station_lon,
                station_lat,
                station_alt,
                field_names,
                field_range_mode,
                grid_x,
                grid_y,
                level_heights,
                product_level_heights,
                lon_0,
                lat_0,
                fillvalue,
                effective_earth_radius,
                max_range_km,
                blind_method,
                use_qc,
                qc_method,
                qc_band,
                qc_use_existing_kdp,
                qc_fallback,
                qc_clear_air_mode,
                qc_clear_air_max_ref,
                qc_clear_air_max_rhohv,
                qc_clear_air_max_phidp_texture,
                qc_clear_air_max_snr,
                quality_weight_terms,
                range_decay_m,
                beam_radius_ref_m,
                vertical_gap_ref_deg,
                qc_weight_field,
                blockage_config,
            )
            for radar_item, station_lon, station_lat, station_alt in zip(
                selected_radars,
                radar_lons,
                radar_lats,
                radar_alts,
            )
        ]

    radar_sites_xy = np.array(
        [[result["radar_x"], result["radar_y"]] for result in worker_results],
        dtype=np.float64,
    )
    actual_max_range_km_by_field = {
        field_name: {
            result["radar_id"]: result["actual_max_range_km_by_field"][field_name]
            for result in worker_results
        }
        for field_name in field_names
    }
    actual_max_range_km_by_sweep_by_field = {
        field_name: {
            result["radar_id"]: result["actual_max_range_km_by_sweep_by_field"][field_name]
            for result in worker_results
        }
        for field_name in field_names
    }
    source_field_by_radar = {
        field_name: {
            result["radar_id"]: result["source_field_by_field"][field_name]
            for result in worker_results
        }
        for field_name in field_names
    }
    qc_applied_by_radar = {
        result["radar_id"]: {
            "applied": bool(result["qc_applied"]),
            "error": result["qc_error"],
        }
        for result in worker_results
    }
    coverage_masks = []
    for result in worker_results:
        if result["observability_mask"] is None:
            coverage_masks.append(np.zeros((len(level_heights),) + grid_x.shape, dtype=bool))
        else:
            coverage_masks.append(np.asarray(result["observability_mask"], dtype=bool))
    coverage_count = np.sum(np.stack(coverage_masks, axis=0), axis=0, dtype=np.int16)
    blind_mask = coverage_count == 0
    data_vars = {}
    weight_cache = _build_composite_weight_cache(
        radar_sites_xy,
        grid_x,
        grid_y,
        method=composite_method,
        influence_radius_m=influence_radius_m,
    )
    for field_name in field_names:
        radar_volumes = [result["field_volumes"][field_name] for result in worker_results]
        field_quality_weights = [result["field_quality_weights"][field_name] for result in worker_results]
        composite, _ = compose_network_volume(
            radar_volumes,
            radar_sites_xy,
            grid_x,
            grid_y,
            fillvalue=fillvalue,
            method=composite_method,
            influence_radius_m=influence_radius_m,
            weight_cache=weight_cache,
            quality_weights=field_quality_weights,
        )
        data_vars[field_name] = composite

    product_dBZ_volume = data_vars.get("dBZ")
    if (
        "dBZ" in field_names
        and product_level_heights.size > 0
        and not np.array_equal(product_level_heights, level_heights)
        and ("VIL" in output_products or "ET" in output_products)
    ):
        product_radar_volumes = []
        product_quality_weights = []
        for result in worker_results:
            if result["product_reflectivity_grid"] is None:
                product_radar_volumes.append(
                    np.full((len(product_level_heights),) + grid_x.shape, float(fillvalue), dtype=np.float64)
                )
                product_quality_weights.append(
                    np.zeros((len(product_level_heights),) + grid_x.shape, dtype=np.float64)
                )
            else:
                product_radar_volumes.append(np.asarray(result["product_reflectivity_grid"], dtype=np.float64))
                product_quality_weights.append(np.asarray(result["product_reflectivity_quality"], dtype=np.float64))
        product_dBZ_volume, _ = compose_network_volume(
            product_radar_volumes,
            radar_sites_xy,
            grid_x,
            grid_y,
            fillvalue=fillvalue,
            method=composite_method,
            influence_radius_m=influence_radius_m,
            weight_cache=weight_cache,
            quality_weights=product_quality_weights,
        )

    dataset = build_network_dataset(
        data_vars,
        lon,
        lat,
        level_heights,
        selected_radars,
        target_time,
        composite_method,
        influence_radius_m,
        fillvalue,
        effective_earth_radius,
        (lon_0, lat_0),
        field_range_mode,
        blind_method,
        max_range_km,
        actual_max_range_km_by_field,
        actual_max_range_km_by_sweep_by_field,
        use_qc,
        qc_method,
        qc_clear_air_mode,
        qc_clear_air_max_ref,
        qc_clear_air_max_rhohv,
        qc_clear_air_max_phidp_texture,
        qc_clear_air_max_snr,
        qc_applied_by_radar,
        source_field_by_radar,
        quality_weight_terms,
        config_path=config_path,
    )
    dataset["coverage_count"] = (("z", "lat", "lon"), coverage_count.astype(np.int16))
    dataset["coverage_count"].attrs = {
        "long_name": "number_of_radars_with_real_observability",
        "coverage_policy": coverage_policy,
        "coordinates": "z lat lon",
    }
    dataset["blind_mask"] = (("z", "lat", "lon"), blind_mask.astype(np.int8))
    dataset["blind_mask"].attrs = {
        "long_name": "blind_zone_mask",
        "flag_values": [0, 1],
        "flag_meanings": "covered blind",
        "coverage_policy": coverage_policy,
        "coordinates": "z lat lon",
    }

    if "dBZ" in data_vars:
        if "CR" in output_products:
            cr_grids = []
            for result in worker_results:
                if result["cr_grid"] is None:
                    cr_grids.append(np.full(grid_x.shape, np.nan, dtype=np.float64))
                else:
                    cr_grids.append(np.asarray(result["cr_grid"], dtype=np.float64))
            cr_stack = np.stack(cr_grids, axis=0)
            cr_valid = np.isfinite(cr_stack) & (cr_stack != fillvalue)
            cr_reduced = np.max(np.where(cr_valid, cr_stack, -np.inf), axis=0)
            cr_composite = np.where(np.any(cr_valid, axis=0), cr_reduced, np.nan)
            dataset["CR"] = (("lat", "lon"), cr_composite)
            dataset["CR"].attrs = {
                "units": "dBZ",
                "long_name": "network_composite_reflectivity",
                "composite_method": "max",
            }
        if "VIL" in output_products:
            dataset["VIL"] = (
                ("lat", "lon"),
                _derive_vil(
                    product_dBZ_volume,
                    product_level_heights,
                    fillvalue=fillvalue,
                    min_dbz=vil_min_dbz,
                    max_dbz_cap=vil_max_dbz_cap,
                ),
            )
            dataset["VIL"].attrs = {
                "units": "kg m-2",
                "long_name": "vertically_integrated_liquid",
                "min_dbz": float(vil_min_dbz),
                "max_dbz_cap": float(vil_max_dbz_cap),
                "source_levels": json.dumps(product_level_heights.tolist()),
                "references": json.dumps(_PRODUCT_REFERENCE_NOTES["VIL"], ensure_ascii=False),
                "implementation_note": (
                    "Uses Greene and Clark (1972) liquid-water relation and a 56 dBZ cap; "
                    "the low-reflectivity cutoff is a configurable pycwr workflow parameter."
                ),
            }
        if "ET" in output_products:
            et_values, et_topped = _derive_et(
                product_dBZ_volume,
                product_level_heights,
                fillvalue=fillvalue,
                threshold_dbz=et_threshold_dbz,
                return_topped=True,
            )
            dataset["ET"] = (
                ("lat", "lon"),
                et_values,
            )
            dataset["ET"].attrs = {
                "units": "m",
                "long_name": "echo_top_height",
                "threshold_dbz": float(et_threshold_dbz),
                "source_levels": json.dumps(product_level_heights.tolist()),
                "references": json.dumps(_PRODUCT_REFERENCE_NOTES["ET"], ensure_ascii=False),
                "implementation_note": (
                    "Highest threshold-exceeding layer plus linear interpolation to the crossing height; "
                    "see ET_TOPPED for columns that still exceed the threshold at the highest sampled level."
                ),
            }
            dataset["ET_TOPPED"] = (("lat", "lon"), np.asarray(et_topped, dtype=np.uint8))
            dataset["ET_TOPPED"].attrs = {
                "units": "1",
                "long_name": "echo_top_topped_flag",
                "flag_values": np.array([0, 1], dtype=np.uint8),
                "flag_meanings": "resolved topped_at_highest_sampled_level",
                "threshold_dbz": float(et_threshold_dbz),
                "source_levels": json.dumps(product_level_heights.tolist()),
                "references": json.dumps(_PRODUCT_REFERENCE_NOTES["ET"], ensure_ascii=False),
            }

    plot_files = {}
    plot_dir = _resolve_plot_output_dir(plot_output_dir, output_path, output_dir)
    if plot_dir is not None:
        plot_dir = str(plot_dir)
        plot_suffix = plot_format.lstrip(".")
        if plot_overview:
            plot_files["overview"] = plot_network_overview(
                selected_radars,
                worker_results,
                lon_min,
                lon_max,
                lat_min,
                lat_max,
                Path(plot_dir) / ("network_overview_%s.%s" % (target_time.strftime("%Y%m%d%H%M%S"), plot_suffix)),
            )
        for field_name in ("CR", "VIL", "ET"):
            if field_name in dataset:
                plot_path = Path(plot_dir) / ("%s_%s.%s" % (field_name.lower(), target_time.strftime("%Y%m%d%H%M%S"), plot_suffix))
                if field_name == "CR" and plot_style == "reference":
                    plot_files[field_name] = _plot_reference_cr_field(
                        lon,
                        lat,
                        dataset[field_name].values,
                        target_time,
                        plot_path,
                        plot_canvas_px=plot_canvas_px,
                    )
                else:
                    cmap = CINRAD_COLORMAP.get("dBZ", "turbo") if field_name == "CR" else "viridis"
                    value_range = CINRAD_field_normvar.get("dBZ", (0.0, 75.0)) if field_name == "CR" else (None, None)
                    plot_files[field_name] = _plot_2d_field(
                        lon,
                        lat,
                        dataset[field_name].values,
                        field_name,
                        plot_path,
                        cmap=cmap,
                        vmin=value_range[0],
                        vmax=value_range[1],
                    )
        for level_height in plot_height_levels:
            if level_heights.size == 0:
                break
            level_index = int(np.argmin(np.abs(level_heights - float(level_height))))
            level_key = "%dm" % int(round(float(level_heights[level_index])))
            level_field_name = str(plot_product_for_levels).upper()
            if level_field_name == "CAPPI" or level_field_name not in dataset:
                level_field_data = dataset["dBZ"].isel(z=level_index).values if "dBZ" in dataset else np.full(grid_x.shape, np.nan)
                level_field_title = "CAPPI"
            else:
                level_field_data = dataset[level_field_name].values
                level_field_title = level_field_name
            plot_files["level_%s" % level_key] = plot_network_level_diagnostic(
                lon,
                lat,
                float(level_heights[level_index]),
                level_field_data,
                dataset["coverage_count"].isel(z=level_index).values,
                dataset["blind_mask"].isel(z=level_index).values.astype(bool),
                Path(plot_dir) / ("level_%s_%s.%s" % (level_key, target_time.strftime("%Y%m%d%H%M%S"), plot_suffix)),
                field_title=level_field_title,
            )
    dataset.attrs["output_products"] = json.dumps(output_products)
    dataset.attrs["coverage_policy"] = coverage_policy
    dataset.attrs["product_level_heights"] = json.dumps(product_level_heights.tolist())
    dataset.attrs["plot_style"] = plot_style
    dataset.attrs["plot_canvas_px"] = json.dumps(list(_resolve_reference_canvas(plot_canvas_px=plot_canvas_px)))
    dataset.attrs["plot_files"] = json.dumps(plot_files, ensure_ascii=False)
    if output_path is None and output_dir:
        output_path = _default_network_output_path(output_dir, target_time, field_names)
    if output_path is not None:
        radar_network_3d_to_netcdf(dataset, output_path)
    return dataset
