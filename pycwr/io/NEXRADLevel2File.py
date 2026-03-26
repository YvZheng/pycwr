# -*- coding: utf-8 -*-
"""Write PRD volumes to NEXRAD Level II archive formats."""

from __future__ import annotations

import datetime
import os
import re
import struct
from collections import OrderedDict

import numpy as np


RECORD_SIZE = 2432
MSG_HEADER_SIZE = 16
MSG31_HEADER_SIZE = 72
VOLUME_HEADER_SIZE = 24
COMPRESSION_RECORD_SIZE = 12

VOLUME_HEADER_FMT = ">9s3sII4s"
MSG1_HEADER_FMT = (
    ">I H h H H H H H "
    "H H H H H H H f "
    "H H H H H "
    "8s 2s 2s 2s "
    "h h h H "
    "32s"
)
MSG1_HEADER_SIZE = struct.calcsize(MSG1_HEADER_FMT)


NEXRAD_MSG31_FIELD_SPECS = OrderedDict(
    [
        ("dBZ", {"moment": b"REF", "scale": 2.0, "offset": 66.0, "word_size": 8}),
        ("V", {"moment": b"VEL", "scale": 2.0, "offset": 129.0, "word_size": 8}),
        ("W", {"moment": b"SW ", "scale": 2.0, "offset": 129.0, "word_size": 8}),
        ("ZDR", {"moment": b"ZDR", "scale": 16.0, "offset": 128.0, "word_size": 8}),
        ("PhiDP", {"moment": b"PHI", "scale": 100.0, "offset": 2.0, "word_size": 16}),
        ("CC", {"moment": b"RHO", "scale": 200.0, "offset": 5.0, "word_size": 8}),
    ]
)
NEXRAD_MSG1_FIELD_SPECS = OrderedDict(
    [
        ("dBZ", {"scale": 2.0, "offset": 66.0}),
        ("V", {"scale": 2.0, "offset": 129.0}),
        ("W", {"scale": 2.0, "offset": 129.0}),
    ]
)


def _ensure_ppi_volume(prd):
    scan_type = str(np.asarray(prd.scan_info["scan_type"].values).item())
    if scan_type != "ppi":
        raise ValueError("NEXRAD Level II export only supports ppi PRD volumes.")


def _ensure_output_path(filename, overwrite=False):
    path = os.path.abspath(filename)
    if os.path.exists(path) and not overwrite:
        raise FileExistsError("Output file already exists: %s" % path)
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    return path


def _python_datetime(value):
    array = np.asarray(value)
    if np.issubdtype(array.dtype, np.datetime64):
        return array.astype("datetime64[ms]").tolist()
    if isinstance(value, datetime.datetime):
        return value
    raise TypeError("Unable to convert %r to datetime." % (value,))


def _compute_mjd(dt):
    return int((dt.date() - datetime.date(1970, 1, 1)).days + 1)


def _compute_milliseconds(dt):
    midnight = datetime.datetime(dt.year, dt.month, dt.day, tzinfo=dt.tzinfo)
    return int((dt - midnight).total_seconds() * 1000)


def _normalize_icao(icao):
    if icao is None:
        return b"PYCW"
    text = str(icao).strip().upper()
    if not text:
        return b"PYCW"
    return text[:4].ljust(4).encode("ascii", "ignore")


def _infer_site_code(prd, fallback="PYCW"):
    sitename = str(getattr(prd, "sitename", "") or "")
    match = re.search(r"([A-Z0-9]{4,8})$", sitename.upper())
    if match is not None:
        return match.group(1)
    return fallback


def _normalize_station_bytes(text, length, fallback):
    raw = (str(text).strip() if text is not None else "") or fallback
    return raw[:length].ljust(length).encode("ascii", "ignore")


def _available_supported_fields(prd, supported_specs, field_names=None, strict=True, required_fields=None):
    available = set(prd.available_fields(range_mode=None))
    def _is_resolvable(field_name):
        if hasattr(prd, "resolve_field_name"):
            return prd.resolve_field_name(field_name, range_mode=None, required=False) is not None
        return field_name in available
    required = tuple(required_fields or ())
    if field_names is None:
        selected = [name for name in supported_specs if _is_resolvable(name)]
    else:
        selected = []
        for name in field_names:
            if name not in supported_specs:
                if strict:
                    raise ValueError("Field %s is not supported by this export format." % name)
                continue
            if not _is_resolvable(name):
                if strict:
                    raise KeyError(name)
                continue
            selected.append(name)
    for name in required:
        if name not in selected:
            if field_names is None:
                raise ValueError("Field %s is required by this export format." % name)
            raise ValueError("Field %s is required by this export format." % name)
    if not selected:
        raise ValueError("No exportable fields are available for this format.")
    return selected


def _sweep_state(isweep, iray, rays_per_sweep, nsweeps):
    if isweep == 0 and iray == 0:
        return 3
    if isweep == nsweeps - 1 and iray == rays_per_sweep - 1:
        return 4
    if iray == 0:
        return 0
    if iray == rays_per_sweep - 1:
        return 2
    return 1


def _get_prd_ray(prd, sweep, field_name, iray):
    source_name = field_name
    if hasattr(prd, "resolve_field_name"):
        source_name = prd.resolve_field_name(field_name, sweep=sweep, range_mode=None, required=False)
    if source_name is None:
        raise KeyError(field_name)
    field = prd.get_sweep_field(sweep, source_name, range_mode=None)
    values = np.asarray(field.values[iray], dtype=np.float32)
    ranges = np.asarray(field["range"].values, dtype=np.float64)
    return values, ranges


def _range_geometry(range_values):
    range_values = np.asarray(range_values, dtype=np.float64)
    if range_values.size == 0:
        raise ValueError("Export range vector cannot be empty.")
    if range_values.size > 1:
        spacing = float(range_values[1] - range_values[0])
        if spacing <= 0.0:
            raise ValueError("Range spacing must be positive.")
    else:
        spacing = float(range_values[0]) if float(range_values[0]) > 0.0 else 1.0
    return int(round(float(range_values[0]))), int(round(spacing))


def _encode_quantized(values, scale, offset, max_code, missing_code=0, word_size=8):
    codes = np.round(np.asarray(values, dtype=np.float64) * float(scale) + float(offset))
    codes[~np.isfinite(values)] = missing_code
    codes = np.clip(codes.astype(np.int64), 0, max_code)
    if word_size == 8:
        return codes.astype(np.uint8).tobytes()
    if word_size == 16:
        return codes.astype(">u2").tobytes()
    raise ValueError("Unsupported word size: %s" % word_size)


def _pack_msg31_moment_block(moment_name, values, first_gate, gate_spacing, scale, offset, word_size):
    encoded = _encode_quantized(values, scale, offset, 255 if word_size == 8 else 65535, word_size=word_size)
    header = struct.pack(
        ">1s3sI H h h h h B B f f",
        b"B",
        moment_name,
        0,
        len(values),
        first_gate,
        gate_spacing,
        0,
        0,
        0,
        word_size,
        float(scale),
        float(offset),
    )
    block = header + encoded
    if len(block) % 2:
        block += b"\x00"
    return block


def _pack_msg31_radial(prd, sweep, iray, field_names, seq_id):
    azimuth = float(prd.fields[sweep].azimuth.values[iray])
    elevation = float(prd.fields[sweep].elevation.values[iray])
    dt = _python_datetime(prd.fields[sweep].time.values[iray])
    mjd = _compute_mjd(dt)
    milliseconds = _compute_milliseconds(dt)
    rays_per_sweep = int(prd.scan_info["rays_per_sweep"].values[sweep])
    radial_state = _sweep_state(sweep, iray, rays_per_sweep, int(prd.nsweeps))
    nyquist = float(prd.scan_info["nyquist_velocity"].values[sweep]) if "nyquist_velocity" in prd.scan_info else 0.0
    altitude = int(round(float(prd.scan_info["altitude"].values)))
    lat = float(prd.scan_info["latitude"].values)
    lon = float(prd.scan_info["longitude"].values)
    vcp = int(np.asarray(prd.metadata.get("volume_number", 0)))

    moment_blocks = []
    max_range_m = 0.0
    for field_name in field_names:
        spec = NEXRAD_MSG31_FIELD_SPECS[field_name]
        values, ranges = _get_prd_ray(prd, sweep, field_name, iray)
        first_gate, gate_spacing = _range_geometry(ranges)
        if ranges.size:
            max_range_m = max(max_range_m, float(ranges[-1]))
        moment_blocks.append(
            _pack_msg31_moment_block(
                spec["moment"],
                values,
                first_gate,
                gate_spacing,
                spec["scale"],
                spec["offset"],
                spec["word_size"],
            )
        )

    msg_header = struct.pack(
        ">HBBHHIHH",
        0,
        1,
        31,
        seq_id & 0x7FFF,
        mjd,
        milliseconds,
        1,
        1,
    )

    ptr_vol = MSG31_HEADER_SIZE
    ptr_elv = ptr_vol + 44
    ptr_rad = ptr_elv + 12
    pointers = [ptr_vol, ptr_elv, ptr_rad]
    next_offset = ptr_rad + 20
    for block in moment_blocks:
        pointers.append(next_offset)
        next_offset += len(block)
    full_pointers = pointers[:10] + [0] * max(0, 10 - len(pointers))

    radial_length = 44 + 12 + 20 + sum(len(block) for block in moment_blocks)
    msg31_header = struct.pack(
        ">4s I H H f B B H B B B B f B b H " + "I" * 10,
        b"RAD1",
        milliseconds,
        mjd,
        iray + 1,
        azimuth,
        0,
        0,
        radial_length,
        0,
        0,
        sweep + 1,
        radial_state,
        elevation,
        0,
        0,
        3 + len(moment_blocks),
        *full_pointers[:10],
    )

    vol_block = struct.pack(
        ">1s3sHBBffhhfffffH2s",
        b"V",
        b"VOL",
        40,
        1,
        0,
        lat,
        lon,
        altitude,
        0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        max(0, min(65535, vcp)),
        b"\x00\x00",
    )
    elv_block = struct.pack(">1s3sHhf", b"E", b"ELV", 12, 0, elevation)
    radial_block = struct.pack(
        ">1s3sHhffh2s",
        b"R",
        b"RAD",
        20,
        int(round(max_range_m / 100.0)),
        0.0,
        0.0,
        int(round(nyquist * 100.0)),
        b"\x00\x00",
    )

    radial = bytearray()
    radial += msg_header
    radial += msg31_header
    radial += vol_block
    radial += elv_block
    radial += radial_block
    for block in moment_blocks:
        radial += block
    radial += struct.pack(">I", 0xFFFFFFFF)
    radial_size_halfwords = (len(radial) - MSG_HEADER_SIZE + 4) // 2
    radial[0:2] = struct.pack(">H", radial_size_halfwords)
    return bytes(radial), (seq_id + 1) & 0x7FFF


def _pack_msg5(prd, icao):
    now = _python_datetime(prd.scan_info["start_time"].values)
    mjd = _compute_mjd(now)
    milliseconds = _compute_milliseconds(now)
    msg_header = struct.pack(
        ">HBBHHIHH",
        0,
        1,
        5,
        0,
        mjd,
        milliseconds,
        1,
        1,
    )
    msg5_main = struct.pack(
        ">H H H H H B B 10s",
        0,
        1,
        0,
        int(prd.nsweeps),
        0,
        2,
        2,
        b"\x00" * 10,
    )
    elev_blocks = bytearray()
    fixed_angles = np.asarray(prd.scan_info["fixed_angle"].values, dtype=np.float64)
    for angle in fixed_angles:
        elev_blocks += struct.pack(
            ">HBBBBHHhhhhhhHHH2sHHH2sHHH2s",
            int(round(float(angle) * 65536.0 / 360.0)),
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            b"\x00\x00",
            0,
            0,
            0,
            b"\x00\x00",
            0,
            0,
            0,
            b"\x00\x00",
        )
    msg5 = bytearray()
    msg5 += msg_header
    msg5 += msg5_main
    msg5 += elev_blocks
    msg5[0:2] = struct.pack(">H", (len(msg5) - 2) // 2)
    if len(msg5) > RECORD_SIZE:
        raise ValueError("MSG5 record exceeds the fixed NEXRAD record size.")
    msg5 += b"\x00" * (RECORD_SIZE - len(msg5))
    return bytes(msg5)


def write_nexrad_level2_msg31(prd, filename, field_names=None, strict=True, overwrite=False, icao=None):
    _ensure_ppi_volume(prd)
    path = _ensure_output_path(filename, overwrite=overwrite)
    selected_fields = _available_supported_fields(prd, NEXRAD_MSG31_FIELD_SPECS, field_names=field_names, strict=strict)
    start_time = _python_datetime(prd.scan_info["start_time"].values)
    volume_header = struct.pack(
        VOLUME_HEADER_FMT,
        b"AR2V0006.",
        b"001",
        _compute_mjd(start_time),
        _compute_milliseconds(start_time),
        _normalize_icao(icao),
    )
    seq_id = 0
    with open(path, "wb") as handle:
        handle.write(volume_header)
        handle.write(b"\x00" * COMPRESSION_RECORD_SIZE)
        handle.write(_pack_msg5(prd, icao))
        for sweep in range(int(prd.nsweeps)):
            rays_per_sweep = int(prd.scan_info["rays_per_sweep"].values[sweep])
            for iray in range(rays_per_sweep):
                radial, seq_id = _pack_msg31_radial(prd, sweep, iray, selected_fields, seq_id)
                handle.write(radial)
    return path


def _pack_msg1_radial(prd, sweep, iray, seq_id):
    rays_per_sweep = int(prd.scan_info["rays_per_sweep"].values[sweep])
    radial_state = _sweep_state(sweep, iray, rays_per_sweep, int(prd.nsweeps))
    dt = _python_datetime(prd.fields[sweep].time.values[iray])
    mjd = _compute_mjd(dt)
    milliseconds = _compute_milliseconds(dt)
    nyquist = float(prd.scan_info["nyquist_velocity"].values[sweep]) if "nyquist_velocity" in prd.scan_info else 0.0
    unambiguous_range = float(prd.scan_info["unambiguous_range"].values[sweep]) if "unambiguous_range" in prd.scan_info else 0.0
    azimuth = float(prd.fields[sweep].azimuth.values[iray])
    elevation = float(prd.fields[sweep].elevation.values[iray])

    ref_values, ref_ranges = _get_prd_ray(prd, sweep, "dBZ", iray)
    vel_values, vel_ranges = _get_prd_ray(prd, sweep, "V", iray)
    sw_values, sw_ranges = _get_prd_ray(prd, sweep, "W", iray)
    first_gate, gate_spacing = _range_geometry(ref_ranges)
    max_gates = (RECORD_SIZE - MSG_HEADER_SIZE - MSG1_HEADER_SIZE - 4) // 3
    ngates = min(len(ref_values), len(vel_values), len(sw_values), max_gates)
    ref_values = ref_values[:ngates]
    vel_values = vel_values[:ngates]
    sw_values = sw_values[:ngates]

    msg_header = struct.pack(
        ">HBBHHIHH",
        1214,
        1,
        1,
        seq_id & 0x7FFF,
        mjd,
        milliseconds,
        1,
        1,
    )
    angle_scale = 180.0 / (4096.0 * 8.0)
    msg1_header = struct.pack(
        MSG1_HEADER_FMT,
        milliseconds,
        mjd,
        int(round(unambiguous_range / 100.0)),
        max(0, min(0xFFFF, int(round(azimuth / angle_scale)))),
        iray + 1,
        radial_state,
        max(0, min(0xFFFF, int(round(elevation / angle_scale)))),
        sweep + 1,
        first_gate,
        first_gate,
        gate_spacing,
        gate_spacing,
        ngates,
        ngates,
        0,
        0.0,
        MSG1_HEADER_SIZE,
        MSG1_HEADER_SIZE + ngates,
        MSG1_HEADER_SIZE + 2 * ngates,
        2,
        0,
        b"\x00" * 8,
        b"\x00" * 2,
        b"\x00" * 2,
        b"\x00" * 2,
        int(round(nyquist * 100.0)),
        0,
        0,
        0,
        b"\x00" * 32,
    )
    radial = bytearray()
    radial += msg_header
    radial += msg1_header
    radial += _encode_quantized(ref_values, 2.0, 66.0, 255, word_size=8)
    radial += _encode_quantized(vel_values, 2.0, 129.0, 255, word_size=8)
    radial += _encode_quantized(sw_values, 2.0, 129.0, 255, word_size=8)
    radial += struct.pack(">I", 0xFFFFFFFF)
    if len(radial) > RECORD_SIZE:
        raise ValueError("MSG1 radial exceeds fixed record size.")
    radial += b"\x00" * (RECORD_SIZE - len(radial))
    return bytes(radial), (seq_id + 1) & 0x7FFF


def write_nexrad_level2_msg1(prd, filename, field_names=None, strict=True, overwrite=False, icao=None):
    _ensure_ppi_volume(prd)
    _available_supported_fields(
        prd,
        NEXRAD_MSG1_FIELD_SPECS,
        field_names=field_names,
        strict=strict,
        required_fields=("dBZ", "V", "W"),
    )
    path = _ensure_output_path(filename, overwrite=overwrite)
    start_time = _python_datetime(prd.scan_info["start_time"].values)
    volume_header = struct.pack(
        VOLUME_HEADER_FMT,
        b"ARCHIVE2 ",
        b"   ",
        _compute_mjd(start_time),
        _compute_milliseconds(start_time),
        _normalize_icao(icao),
    )
    compression_record = struct.pack(">I", 0) + struct.pack(">H", RECORD_SIZE) + (b"\x00" * 6)
    seq_id = 0
    with open(path, "wb") as handle:
        handle.write(volume_header)
        handle.write(compression_record)
        for sweep in range(int(prd.nsweeps)):
            rays_per_sweep = int(prd.scan_info["rays_per_sweep"].values[sweep])
            for iray in range(rays_per_sweep):
                radial, seq_id = _pack_msg1_radial(prd, sweep, iray, seq_id)
                handle.write(radial)
    return path
