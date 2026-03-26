"""Interop helpers for exporting PRD objects to Py-ART and xradar workflows."""

from __future__ import annotations

import copy
import importlib
from collections import OrderedDict

import numpy as np
import xarray as xr

from ..configure.default_config import CINRAD_field_mapping, DEFAULT_METADATA
from ..configure.pyart_config import get_fillvalue, get_metadata
from ..io.util import date2num, make_time_unit_str
from .PyartRadar import Radar as InternalRadar


PYART_EXTRA_MESSAGE = (
    "Py-ART export requires the optional dependencies from "
    "`pip install \"pycwr[full]\"`."
)
XRADAR_EXTRA_MESSAGE = (
    "xradar export requires the optional dependencies from "
    "`pip install \"pycwr[full]\"`."
)
XRADAR_FIELD_MAPPING = {
    "dBZ": "DBZ",
    "dBT": "DBTH",
    "V": "VR",
    "W": "WRADH",
    "ZDR": "ZDR",
    "CC": "RHOHV",
    "PhiDP": "PHIDP",
    "KDP": "KDP",
}
_PREFERRED_EXPORT_FIELDS = (
    "dBZ",
    "dBT",
    "V",
    "W",
    "ZDR",
    "CC",
    "PhiDP",
    "KDP",
    "SQI",
    "CPA",
    "LDR",
    "HCL",
    "CF",
    "SNRH",
    "SNRV",
)


def _import_external_pyart_radar_class():
    """Return the upstream Py-ART Radar class when available."""
    for module_name in ("pyart", "arm_pyart"):
        try:
            module = importlib.import_module(module_name)
        except ImportError:
            continue
        core = getattr(module, "core", None)
        radar_class = getattr(core, "Radar", None) if core is not None else None
        if radar_class is not None:
            return radar_class
        try:
            core_module = importlib.import_module("%s.core" % module_name)
        except ImportError:
            continue
        radar_class = getattr(core_module, "Radar", None)
        if radar_class is not None:
            return radar_class
    return None


def _import_datatree_class():
    """Return a DataTree implementation when available."""
    try:
        from xarray import DataTree

        return DataTree
    except ImportError:
        pass
    except AttributeError:
        pass
    try:
        from datatree import DataTree

        return DataTree
    except ImportError:
        return None


def _import_xradar_model():
    """Return the installed xradar.model module when available."""
    try:
        return importlib.import_module("xradar.model")
    except ImportError:
        return None


def resolve_pyart_radar_class(use_external=None, strict=False):
    """
    Resolve the Radar class to use for export.

    Parameters
    ----------
    use_external : bool or None
        ``True`` forces an upstream Py-ART Radar, ``False`` forces the bundled
        compatibility Radar, and ``None`` prefers upstream when installed.
    strict : bool
        When ``True`` and ``use_external`` requests upstream Py-ART, raise a
        clear ``ImportError`` if the optional dependency is missing.
    """
    external_radar_class = _import_external_pyart_radar_class()
    if use_external is False:
        return InternalRadar, False
    if external_radar_class is not None:
        return external_radar_class, True
    if use_external is True and strict:
        raise ImportError(PYART_EXTRA_MESSAGE)
    return InternalRadar, False


def _deepcopy_radar_init_kwargs(radar):
    """Create constructor kwargs from a Radar-like object."""
    return {
        "time": copy.deepcopy(radar.time),
        "_range": copy.deepcopy(radar.range),
        "fields": copy.deepcopy(radar.fields),
        "metadata": copy.deepcopy(radar.metadata),
        "scan_type": copy.deepcopy(radar.scan_type),
        "latitude": copy.deepcopy(radar.latitude),
        "longitude": copy.deepcopy(radar.longitude),
        "altitude": copy.deepcopy(radar.altitude),
        "sweep_number": copy.deepcopy(radar.sweep_number),
        "sweep_mode": copy.deepcopy(radar.sweep_mode),
        "fixed_angle": copy.deepcopy(radar.fixed_angle),
        "sweep_start_ray_index": copy.deepcopy(radar.sweep_start_ray_index),
        "sweep_end_ray_index": copy.deepcopy(radar.sweep_end_ray_index),
        "azimuth": copy.deepcopy(radar.azimuth),
        "elevation": copy.deepcopy(radar.elevation),
        "altitude_agl": copy.deepcopy(radar.altitude_agl),
        "target_scan_rate": copy.deepcopy(radar.target_scan_rate),
        "rays_are_indexed": copy.deepcopy(radar.rays_are_indexed),
        "ray_angle_res": copy.deepcopy(radar.ray_angle_res),
        "scan_rate": copy.deepcopy(radar.scan_rate),
        "antenna_transition": copy.deepcopy(radar.antenna_transition),
        "instrument_parameters": copy.deepcopy(radar.instrument_parameters),
        "radar_calibration": copy.deepcopy(radar.radar_calibration),
        "rotation": copy.deepcopy(radar.rotation),
        "tilt": copy.deepcopy(radar.tilt),
        "roll": copy.deepcopy(radar.roll),
        "drift": copy.deepcopy(radar.drift),
        "heading": copy.deepcopy(radar.heading),
        "pitch": copy.deepcopy(radar.pitch),
        "georefs_applied": copy.deepcopy(radar.georefs_applied),
    }


def clone_radar_to_class(radar, radar_class):
    """Clone a Radar-like object into the target Radar class."""
    if isinstance(radar, radar_class):
        return radar
    return radar_class(**_deepcopy_radar_init_kwargs(radar))


def _invert_field_mapping():
    inverse = {}
    for source_name, target_name in CINRAD_field_mapping.items():
        if target_name is not None and target_name not in inverse:
            inverse[target_name] = source_name
    return inverse


_INVERSE_FIELD_MAPPING = _invert_field_mapping()


def _normalize_prd_field_names(prd, field_names=None):
    """Resolve PRD field names from raw or Py-ART-style names."""
    if hasattr(prd, "available_fields"):
        available = set(prd.available_fields(range_mode=None))
    else:
        available = {name for sweep in prd.fields for name in sweep.data_vars}
    def _resolve(name):
        requested = str(name)
        if hasattr(prd, "resolve_field_name"):
            resolved = prd.resolve_field_name(requested, range_mode=None, required=False)
            if resolved is not None:
                return requested, resolved
        if requested in available:
            return requested, requested
        mapped = _INVERSE_FIELD_MAPPING.get(requested)
        if mapped is None:
            return None, None
        if hasattr(prd, "resolve_field_name"):
            resolved = prd.resolve_field_name(mapped, range_mode=None, required=False)
            if resolved is not None:
                return mapped, resolved
        if mapped in available:
            return mapped, mapped
        return None, None
    if field_names is None:
        ordered = []
        consumed = set()
        for logical_name in _PREFERRED_EXPORT_FIELDS:
            requested, resolved = _resolve(logical_name)
            if requested is None or requested in ordered:
                continue
            ordered.append(requested)
            consumed.add(resolved)
        for sweep in prd.fields:
            for name in sweep.data_vars:
                if name in available and name not in ordered and name not in consumed:
                    ordered.append(name)
        return ordered

    resolved = []
    for name in field_names:
        requested, source = _resolve(name)
        if requested is not None:
            resolved.append(requested)
            continue
        raise KeyError(name)
    return resolved


def _resolve_export_source_field(prd, sweep, field_name, range_mode=None):
    """Resolve the actual PRD field used to satisfy one export field."""
    if hasattr(prd, "resolve_field_name"):
        return prd.resolve_field_name(field_name, sweep=sweep, range_mode=range_mode, required=False)
    available = set(prd.available_fields(sweep=sweep, range_mode=range_mode))
    return str(field_name) if str(field_name) in available else None


def _canonical_field_name(field_name):
    """Return the Py-ART field name for a PRD field."""
    return CINRAD_field_mapping.get(field_name) or field_name


def _xradar_field_name(field_name):
    """Return the xradar/CfRadial-style moment name for a PRD field."""
    return XRADAR_FIELD_MAPPING.get(field_name, _canonical_field_name(field_name))


def _isoformat_datetime(value):
    """Return an ISO-8601 string for a scalar datetime-like value."""
    array = np.asarray(value)
    if np.issubdtype(array.dtype, np.datetime64):
        return np.datetime_as_string(array.astype("datetime64[ns]"), unit="s") + "Z"
    return str(value)


def _xradar_coordinate_attrs(prd, common_range):
    """Return xradar-aligned attrs for standard sweep coordinates."""
    model = _import_xradar_model()
    first_gate = float(common_range[0]) if common_range.size else 0.0
    if common_range.size > 1:
        gate_spacing = float(common_range[1] - common_range[0])
    else:
        gate_spacing = first_gate
    if model is None:
        range_attrs = copy.deepcopy(DEFAULT_METADATA["range"])
        azimuth_attrs = copy.deepcopy(DEFAULT_METADATA["azimuth"])
        elevation_attrs = copy.deepcopy(DEFAULT_METADATA["elevation"])
        latitude_attrs = copy.deepcopy(DEFAULT_METADATA["latitude"])
        longitude_attrs = copy.deepcopy(DEFAULT_METADATA["longitude"])
        altitude_attrs = copy.deepcopy(DEFAULT_METADATA["altitude"])
        time_attrs = copy.deepcopy(DEFAULT_METADATA["time"])
        time_attrs["units"] = "seconds since %s" % _isoformat_datetime(prd.scan_info["start_time"].values)
        range_attrs["meters_to_center_of_first_gate"] = first_gate
        range_attrs["meters_between_gates"] = gate_spacing
        return range_attrs, time_attrs, azimuth_attrs, elevation_attrs, latitude_attrs, longitude_attrs, altitude_attrs
    range_attrs = model.get_range_attrs()
    range_attrs["meters_to_center_of_first_gate"] = first_gate
    range_attrs["meters_between_gates"] = gate_spacing
    time_attrs = model.get_time_attrs(_isoformat_datetime(prd.scan_info["start_time"].values))
    return (
        range_attrs,
        time_attrs,
        model.get_azimuth_attrs(),
        model.get_elevation_attrs(),
        model.get_latitude_attrs(),
        model.get_longitude_attrs(),
        model.get_altitude_attrs(),
    )


def _xradar_moment_attrs(field_name, fallback_attrs):
    """Return xradar-aligned attrs for a radar moment field."""
    model = _import_xradar_model()
    target_name = _xradar_field_name(field_name)
    if model is not None and target_name in getattr(model, "sweep_vars_mapping", {}):
        attrs = copy.deepcopy(model.sweep_vars_mapping[target_name])
        attrs["coordinates"] = "elevation azimuth range time"
        return attrs
    attrs = copy.deepcopy(fallback_attrs)
    attrs["coordinates"] = "elevation azimuth range time"
    return attrs


def _to_python_datetimes(values):
    """Convert a datetime coordinate array to Python datetimes."""
    array = np.asarray(values)
    if np.issubdtype(array.dtype, np.datetime64):
        return array.astype("datetime64[us]").astype(object)
    return array.astype(object)


def _resolve_range_index(common_range, field_range):
    """Return the common-range indices matching a field range vector."""
    common = np.asarray(common_range, dtype=np.float64)
    field = np.asarray(field_range, dtype=np.float64)
    if field.size == 0:
        return np.array([], dtype=np.intp)
    if common.size < field.size:
        raise ValueError("common_range is shorter than field_range")
    if np.allclose(common[: field.size], field, rtol=0.0, atol=1.0e-6):
        return np.arange(field.size, dtype=np.intp)
    indices = np.searchsorted(common, field)
    if np.any(indices >= common.size):
        raise ValueError("field_range extends beyond the common Py-ART range grid")
    if not np.allclose(common[indices], field, rtol=0.0, atol=1.0e-6):
        raise ValueError("field_range is incompatible with the common Py-ART range grid")
    return indices.astype(np.intp)


def _resolve_common_range(prd, field_names, range_mode):
    """Return a shared range vector suitable for Py-ART-style export."""
    best_range = None
    for sweep in prd.scan_info.sweep.values:
        for field_name in field_names:
            source_name = _resolve_export_source_field(prd, int(sweep), field_name, range_mode=range_mode)
            if source_name is None:
                continue
            field = prd.get_sweep_field(sweep, source_name, range_mode=range_mode)
            candidate = np.asarray(field.range.values, dtype=np.float32)
            if best_range is None or candidate.size > best_range.size:
                best_range = candidate
    if best_range is None:
        raise ValueError("No exportable fields were found in the PRD volume.")
    return best_range


def _default_scan_metadata(prd):
    """Return common metadata dictionaries built from a PRD object."""
    scan_type = str(np.asarray(prd.scan_info["scan_type"].values).item())
    latitude = get_metadata("latitude")
    longitude = get_metadata("longitude")
    altitude = get_metadata("altitude")
    latitude["data"] = np.array([float(np.asarray(prd.scan_info["latitude"].values))], dtype="float64")
    longitude["data"] = np.array([float(np.asarray(prd.scan_info["longitude"].values))], dtype="float64")
    altitude["data"] = np.array([float(np.asarray(prd.scan_info["altitude"].values))], dtype="float64")

    metadata = get_metadata("metadata")
    metadata["site_name"] = prd.sitename
    metadata["instrument_name"] = prd.sitename
    metadata["source"] = "pycwr.PRD"
    metadata["range_mode"] = "aligned"

    sweep_mode = get_metadata("sweep_mode")
    if scan_type == "ppi":
        sweep_mode["data"] = np.array(prd.nsweeps * ["azimuth_surveillance"], dtype="S")
    elif scan_type == "rhi":
        sweep_mode["data"] = np.array(prd.nsweeps * ["rhi"], dtype="S")
    else:
        sweep_mode["data"] = np.array(prd.nsweeps * ["sector"], dtype="S")

    sweep_number = get_metadata("sweep_number")
    sweep_number["data"] = np.arange(prd.nsweeps, dtype="int32")

    fixed_angle = get_metadata("fixed_angle")
    fixed_angle["data"] = np.asarray(prd.scan_info["fixed_angle"].values, dtype="float32")

    return scan_type, latitude, longitude, altitude, metadata, sweep_mode, sweep_number, fixed_angle


def _build_instrument_parameters(prd, template_radar=None):
    """Create instrument parameter dictionaries for a PRD export."""
    if template_radar is not None and getattr(template_radar, "instrument_parameters", None) is not None:
        return copy.deepcopy(template_radar.instrument_parameters)

    rays_per_sweep = np.asarray(prd.scan_info["rays_per_sweep"].values, dtype=np.int32)
    nyquist_velocity = get_metadata("nyquist_velocity")
    nyquist_velocity["data"] = np.repeat(
        np.asarray(prd.scan_info["nyquist_velocity"].values, dtype=np.float32),
        rays_per_sweep,
    )
    unambiguous_range = get_metadata("unambiguous_range")
    unambiguous_range["data"] = np.repeat(
        np.asarray(prd.scan_info["unambiguous_range"].values, dtype=np.float32),
        rays_per_sweep,
    )
    frequency = get_metadata("frequency")
    frequency["data"] = np.array([float(np.asarray(prd.scan_info["frequency"].values))], dtype=np.float32)
    return {
        "nyquist_velocity": nyquist_velocity,
        "unambiguous_range": unambiguous_range,
        "frequency": frequency,
    }


def build_radar_from_prd(prd, range_mode=None, field_names=None, radar_class=None, template_radar=None):
    """Build a Radar-like object directly from a PRD volume."""
    field_names = _normalize_prd_field_names(prd, field_names)
    if hasattr(prd, "_resolve_export_range_mode"):
        range_mode = prd._resolve_export_range_mode(field_names=field_names, range_mode=range_mode)
    common_range = _resolve_common_range(prd, field_names, range_mode)

    field_buffers = {}
    sweep_starts = []
    sweep_ends = []
    time_chunks = []
    azimuth_chunks = []
    elevation_chunks = []
    ray_offset = 0

    for field_name in field_names:
        field_buffers[_canonical_field_name(field_name)] = np.full(
            (prd.nrays, common_range.size),
            np.nan,
            dtype=np.float32,
        )

    for sweep in prd.scan_info.sweep.values:
        sweep_dataset = prd.fields[int(sweep)]
        sweep_size = int(sweep_dataset.sizes["time"])
        sweep_starts.append(ray_offset)
        sweep_ends.append(ray_offset + sweep_size - 1)
        time_chunks.append(_to_python_datetimes(sweep_dataset["time"].values))
        azimuth_chunks.append(np.asarray(sweep_dataset["azimuth"].values, dtype=np.float32))
        elevation_chunks.append(np.asarray(sweep_dataset["elevation"].values, dtype=np.float32))

        for field_name in field_names:
            source_name = _resolve_export_source_field(prd, int(sweep), field_name, range_mode=range_mode)
            if source_name is None:
                continue
            field = prd.get_sweep_field(sweep, source_name, range_mode=range_mode)
            target_name = _canonical_field_name(field_name)
            range_index = _resolve_range_index(common_range, field.range.values)
            field_buffers[target_name][ray_offset : ray_offset + sweep_size, range_index] = np.asarray(
                field.values,
                dtype=np.float32,
            )
        ray_offset += sweep_size

    dts = np.concatenate(time_chunks)
    units = make_time_unit_str(min(dts))
    time = get_metadata("time")
    time["units"] = units
    time["data"] = date2num(dts, units).astype("float32")

    _range = get_metadata("range")
    _range["data"] = common_range.astype("float32")
    if common_range.size:
        _range["meters_to_center_of_first_gate"] = float(common_range[0])
    if common_range.size > 1:
        _range["meters_between_gates"] = float(common_range[1] - common_range[0])

    azimuth = get_metadata("azimuth")
    azimuth["data"] = np.concatenate(azimuth_chunks).astype("float32")
    elevation = get_metadata("elevation")
    elevation["data"] = np.concatenate(elevation_chunks).astype("float32")

    sweep_start_ray_index = get_metadata("sweep_start_ray_index")
    sweep_start_ray_index["data"] = np.asarray(sweep_starts, dtype="int32")
    sweep_end_ray_index = get_metadata("sweep_end_ray_index")
    sweep_end_ray_index["data"] = np.asarray(sweep_ends, dtype="int32")

    scan_type, latitude, longitude, altitude, metadata, sweep_mode, sweep_number, fixed_angle = _default_scan_metadata(prd)
    metadata["range_mode"] = range_mode

    fields = {}
    fill_value = get_fillvalue()
    for field_name in field_names:
        target_name = _canonical_field_name(field_name)
        field_metadata = get_metadata(target_name)
        field_metadata["data"] = np.ma.masked_array(
            field_buffers[target_name],
            mask=np.isnan(field_buffers[target_name]),
            fill_value=fill_value,
        )
        field_metadata["_FillValue"] = fill_value
        fields[target_name] = field_metadata

    instrument_parameters = _build_instrument_parameters(prd, template_radar=template_radar)
    radar_class = InternalRadar if radar_class is None else radar_class
    return radar_class(
        time=time,
        _range=_range,
        fields=fields,
        metadata=metadata,
        scan_type=scan_type,
        latitude=latitude,
        longitude=longitude,
        altitude=altitude,
        sweep_number=sweep_number,
        sweep_mode=sweep_mode,
        fixed_angle=fixed_angle,
        sweep_start_ray_index=sweep_start_ray_index,
        sweep_end_ray_index=sweep_end_ray_index,
        azimuth=azimuth,
        elevation=elevation,
        altitude_agl=copy.deepcopy(getattr(template_radar, "altitude_agl", None)),
        target_scan_rate=copy.deepcopy(getattr(template_radar, "target_scan_rate", None)),
        rays_are_indexed=copy.deepcopy(getattr(template_radar, "rays_are_indexed", None)),
        ray_angle_res=copy.deepcopy(getattr(template_radar, "ray_angle_res", None)),
        scan_rate=copy.deepcopy(getattr(template_radar, "scan_rate", None)),
        antenna_transition=copy.deepcopy(getattr(template_radar, "antenna_transition", None)),
        instrument_parameters=instrument_parameters,
        radar_calibration=copy.deepcopy(getattr(template_radar, "radar_calibration", None)),
        rotation=copy.deepcopy(getattr(template_radar, "rotation", None)),
        tilt=copy.deepcopy(getattr(template_radar, "tilt", None)),
        roll=copy.deepcopy(getattr(template_radar, "roll", None)),
        drift=copy.deepcopy(getattr(template_radar, "drift", None)),
        heading=copy.deepcopy(getattr(template_radar, "heading", None)),
        pitch=copy.deepcopy(getattr(template_radar, "pitch", None)),
        georefs_applied=copy.deepcopy(getattr(template_radar, "georefs_applied", None)),
    )


def export_pyart_radar(
    prd,
    existing_radar=None,
    range_mode=None,
    field_names=None,
    use_external=None,
    strict=False,
):
    """Export a PRD volume to an upstream or bundled Radar object."""
    if hasattr(prd, "_resolve_export_range_mode"):
        range_mode = prd._resolve_export_range_mode(field_names=field_names, range_mode=range_mode)
    radar_class, is_external = resolve_pyart_radar_class(use_external=use_external, strict=strict)
    if range_mode == "aligned" and field_names is None and existing_radar is not None:
        return clone_radar_to_class(existing_radar, radar_class)
    return build_radar_from_prd(
        prd,
        range_mode=range_mode,
        field_names=field_names,
        radar_class=radar_class,
        template_radar=existing_radar,
    )


def build_xradar_sweep_datasets(prd, range_mode=None, field_names=None):
    """Build sweep-level xarray datasets for xradar-style workflows."""
    field_names = _normalize_prd_field_names(prd, field_names)
    if hasattr(prd, "_resolve_export_range_mode"):
        range_mode = prd._resolve_export_range_mode(field_names=field_names, range_mode=range_mode)
    datasets = OrderedDict()
    model = _import_xradar_model()

    for sweep in prd.scan_info.sweep.values:
        sweep = int(sweep)
        candidate_fields = []
        active_field_names = []
        for field_name in field_names:
            source_name = _resolve_export_source_field(prd, sweep, field_name, range_mode=range_mode)
            if source_name is None:
                continue
            candidate_fields.append(prd.get_sweep_field(sweep, source_name, range_mode=range_mode))
            active_field_names.append(field_name)
        if candidate_fields:
            reference_field = max(candidate_fields, key=lambda field: field.range.size)
        else:
            fallback_name = next(iter(prd.fields[sweep].data_vars))
            reference_field = prd.get_sweep_field(sweep, fallback_name, range_mode=range_mode)
        common_range = np.asarray(reference_field.range.values, dtype=np.float32)
        (
            range_attrs,
            time_attrs,
            azimuth_attrs,
            elevation_attrs,
            latitude_attrs,
            longitude_attrs,
            altitude_attrs,
        ) = _xradar_coordinate_attrs(prd, common_range)
        dataset = xr.Dataset(
            coords={
                "time": ("time", np.asarray(reference_field.time.values)),
                "azimuth": ("time", np.asarray(reference_field.azimuth.values, dtype=np.float32)),
                "elevation": ("time", np.asarray(reference_field.elevation.values, dtype=np.float32)),
                "range": ("range", common_range),
                "frequency": ("frequency", np.array([float(np.asarray(prd.scan_info["frequency"].values))], dtype=np.float32)),
            },
            attrs={
                "site_name": prd.metadata["site_name"],
                "range_mode": range_mode,
            },
        )
        dataset = dataset.assign_coords(
            latitude=xr.DataArray(float(np.asarray(prd.scan_info["latitude"].values))),
            longitude=xr.DataArray(float(np.asarray(prd.scan_info["longitude"].values))),
            altitude=xr.DataArray(float(np.asarray(prd.scan_info["altitude"].values))),
        )
        dataset["sweep_number"] = xr.DataArray(np.int32(sweep))
        dataset["sweep_mode"] = xr.DataArray(str(np.asarray(prd.scan_info["sweep_mode"].values[sweep])))
        dataset["follow_mode"] = xr.DataArray(str(np.asarray(prd.scan_info["follow_mode"].values[sweep])))
        dataset["prt_mode"] = xr.DataArray(str(np.asarray(prd.scan_info["prt_mode"].values[sweep])))
        dataset["sweep_fixed_angle"] = xr.DataArray(float(np.asarray(prd.scan_info["fixed_angle"].values[sweep])), attrs=copy.deepcopy(DEFAULT_METADATA["fixed_angle"]))
        dataset["nyquist_velocity"] = xr.DataArray(float(np.asarray(prd.scan_info["nyquist_velocity"].values[sweep])), attrs=copy.deepcopy(DEFAULT_METADATA["nyquist_velocity"]))
        dataset["unambiguous_range"] = xr.DataArray(float(np.asarray(prd.scan_info["unambiguous_range"].values[sweep])), attrs=copy.deepcopy(DEFAULT_METADATA["unambiguous_range"]))
        dataset["latitude"].attrs = latitude_attrs
        dataset["longitude"].attrs = longitude_attrs
        dataset["altitude"].attrs = altitude_attrs
        dataset["time"].attrs = time_attrs
        dataset["azimuth"].attrs = azimuth_attrs
        dataset["elevation"].attrs = elevation_attrs
        dataset["range"].attrs = range_attrs
        dataset["frequency"].attrs = {"standard_name": "", "units": "s-1"}
        for field_name, field in zip(active_field_names, candidate_fields):
            field_data = np.full(reference_field.shape, np.nan, dtype=np.float32)
            range_index = _resolve_range_index(common_range, field.range.values)
            field_data[:, range_index] = np.asarray(field.values, dtype=np.float32)
            dataset[_xradar_field_name(field_name)] = xr.DataArray(
                field_data,
                dims=("time", "range"),
                attrs=_xradar_moment_attrs(field_name, field.attrs),
            )
            dataset[_xradar_field_name(field_name)].encoding["_FillValue"] = get_fillvalue()
        if model is not None:
            dataset = model.conform_cfradial2_sweep_group(dataset, optional=True)
            dataset.attrs["site_name"] = prd.metadata["site_name"]
            dataset.attrs["range_mode"] = range_mode
        datasets["sweep_%d" % sweep] = dataset
    return datasets


def export_xradar_tree(prd, range_mode=None, field_names=None, strict=True):
    """Export a PRD volume as a DataTree suitable for xradar-style workflows."""
    if hasattr(prd, "_resolve_export_range_mode"):
        range_mode = prd._resolve_export_range_mode(field_names=field_names, range_mode=range_mode)
    datatree_class = _import_datatree_class()
    if datatree_class is None:
        if strict:
            raise ImportError(XRADAR_EXTRA_MESSAGE)
        return build_xradar_sweep_datasets(prd, range_mode=range_mode, field_names=field_names)

    root_attrs = {
        "Conventions": prd.metadata["Conventions"],
        "version": prd.metadata["version"],
        "title": prd.metadata["title"],
        "instrument_name": prd.metadata["instrument_name"],
        "institution": prd.metadata["institution"],
        "references": prd.metadata["references"],
        "source": prd.metadata["source"],
        "history": prd.metadata["history"],
        "comment": prd.metadata["comment"],
        "platform_is_mobile": prd.metadata["platform_is_mobile"],
        "site_name": prd.metadata["site_name"],
        "scan_name": prd.metadata["scan_name"],
        "scan_id": prd.metadata["scan_id"],
        "ray_times_increase": prd.metadata["ray_times_increase"],
        "simulated": prd.metadata["simulated"],
    }
    root = datatree_class(
        xr.Dataset(
            data_vars={
                "volume_number": xr.DataArray(np.int32(prd.metadata["volume_number"])),
                "time_coverage_start": xr.DataArray(_isoformat_datetime(prd.scan_info["start_time"].values)),
                "time_coverage_end": xr.DataArray(_isoformat_datetime(prd.scan_info["end_time"].values)),
                "latitude": xr.DataArray(float(np.asarray(prd.scan_info["latitude"].values))),
                "longitude": xr.DataArray(float(np.asarray(prd.scan_info["longitude"].values))),
                "altitude": xr.DataArray(float(np.asarray(prd.scan_info["altitude"].values))),
                "sweep_group_name": xr.DataArray(np.asarray(prd.scan_info["sweep_group_name"].values, dtype=object), dims=("sweep",)),
                "sweep_fixed_angle": xr.DataArray(np.asarray(prd.scan_info["fixed_angle"].values, dtype=np.float32), dims=("sweep",)),
                "platform_type": xr.DataArray(str(prd.metadata["platform_type"])),
                "instrument_type": xr.DataArray(str(prd.metadata["instrument_type"])),
                "primary_axis": xr.DataArray(str(prd.metadata["primary_axis"])),
                "frequency": xr.DataArray(np.array([float(np.asarray(prd.scan_info["frequency"].values))], dtype=np.float32), dims=("frequency",)),
            },
            coords={"sweep": np.arange(prd.nsweeps, dtype=np.int32)},
            attrs={
                **root_attrs,
                "range_mode": range_mode,
            }
        ),
        name="root",
    )
    root.dataset["latitude"].attrs = _xradar_coordinate_attrs(prd, np.array([0.0], dtype=np.float32))[4]
    root.dataset["longitude"].attrs = _xradar_coordinate_attrs(prd, np.array([0.0], dtype=np.float32))[5]
    root.dataset["altitude"].attrs = _xradar_coordinate_attrs(prd, np.array([0.0], dtype=np.float32))[6]
    root.dataset["sweep_fixed_angle"].attrs = copy.deepcopy(DEFAULT_METADATA["fixed_angle"])
    for name, dataset in build_xradar_sweep_datasets(prd, range_mode=range_mode, field_names=field_names).items():
        root[name] = datatree_class(dataset, name=name)
    return root
