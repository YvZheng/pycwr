# -*- coding: utf-8 -*-
"""Lightweight local Flask app replacing the legacy Qt radar viewer."""

from __future__ import annotations

import io
import json
import os
import re
import secrets
import threading
import webbrowser
from collections import OrderedDict
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from flask import Flask, jsonify, render_template, request, send_file

from ..configure.default_config import HYDROMETEOR_CLASS_NAMES_ZH
from ..configure.location_config import last_open_dir
from ..draw._plot_core import (
    CartesianReferenceOptions,
    MapOptions,
    ccrs,
    default_sweep_title,
    plot_cartesian,
    plot_map,
    plot_vertical_section,
    resolve_field_data,
)
from ..interp import parse_radar_time_from_filename
from ..io import read_auto
from ..io.util import radar_format
from .web_colors import SPECIAL_COLORS, build_web_style


SUPPORTED_PATTERNS = (".a", ".v", ".bin", ".ar2", ".bin.bz2", ".ar2.bz2", ".gz")
BLOCKED_ARCHIVE_PATTERNS = (".tar.gz", ".tgz", ".zip", ".whl")
PREFERRED_FIELDS = ("dBZ", "V", "W", "ZDR", "KDP", "CC", "SNRH", "SNRV", "SQI", "RR", "LDR", "PhiDP")
IGNORED_DIRECTORIES = {"__pycache__", ".git", ".hg", ".svn", "dist", "build", ".venv", "venv", "node_modules"}
IGNORED_SUFFIXES = {".py", ".pyc", ".pyo", ".so", ".pyd", ".dll", ".dylib", ".md", ".txt", ".json", ".yaml", ".yml", ".toml", ".ini", ".cfg", ".log", ".html", ".css", ".js"}
_STATION_CODE_RE = re.compile(r"(Z[A-Z]?\d{3,4})", re.IGNORECASE)


class RadarFileCache(object):
    def __init__(self, max_items=6):
        self.max_items = max_items
        self._items = OrderedDict()

    def get(self, file_path):
        file_path = os.path.abspath(file_path)
        mtime = os.path.getmtime(file_path)
        cached = self._items.get(file_path)
        if cached is not None and cached["mtime"] == mtime:
            self._items.move_to_end(file_path)
            return cached["radar"]
        radar = read_auto(file_path)
        self._items[file_path] = {"mtime": mtime, "radar": radar}
        self._items.move_to_end(file_path)
        while len(self._items) > self.max_items:
            self._items.popitem(last=False)
        return radar


def _is_relative_to(path, root):
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def _normalize_allowed_roots(allowed_roots):
    if allowed_roots is None:
        return None
    roots = []
    for item in allowed_roots or []:
        if not item:
            continue
        root = Path(item).expanduser().resolve(strict=False)
        if root.is_dir() and root not in roots:
            roots.append(root)
    if not roots:
        roots.append(Path(_default_directory()).expanduser().resolve(strict=False))
    return roots


def _resolve_within_roots(raw_path, allowed_roots, expect_directory=False):
    if not raw_path:
        raise PermissionError("Path access is not allowed.")
    candidate = Path(raw_path).expanduser().resolve(strict=False)
    if allowed_roots is not None and not any(_is_relative_to(candidate, root) for root in allowed_roots):
        raise PermissionError("Requested path is outside the allowed roots.")
    if expect_directory:
        if not candidate.is_dir():
            raise FileNotFoundError("Directory does not exist.")
    else:
        if not candidate.is_file():
            raise FileNotFoundError("File does not exist.")
    return str(candidate)


def _default_directory():
    try:
        with open(last_open_dir, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
            candidate = payload.get("lastOpenDir")
            if candidate and os.path.isdir(candidate):
                return candidate
    except Exception:
        pass
    cwd = os.getcwd()
    for candidate in (cwd, os.path.dirname(cwd)):
        try:
            if _scan_tree(candidate).get("file_count", 0) > 0:
                return candidate
        except Exception:
            continue
    return cwd


def _scan_directory(directory):
    directory = os.path.abspath(directory)
    if not os.path.isdir(directory):
        raise FileNotFoundError("Directory does not exist.")
    files = []
    for item in sorted(Path(directory).iterdir(), key=lambda path: path.name.lower()):
        if not item.is_file():
            continue
        info = _supported_file_info(item)
        if info is not None:
            files.append(info)
    return files


def _supported_file_info(item):
    if not _matches_supported_filename(item.name):
        return None
    try:
        detected_format = radar_format(str(item))
    except Exception:
        detected_format = None
    return {
        "type": "file",
        "name": item.name,
        "path": str(item.resolve()),
        "size": item.stat().st_size,
        "format": detected_format,
    }


def _matches_supported_filename(name):
    lower_name = str(name).lower()
    if any(lower_name.endswith(pattern) for pattern in BLOCKED_ARCHIVE_PATTERNS):
        return False
    if Path(lower_name).suffix in IGNORED_SUFFIXES:
        return False
    return any(lower_name.endswith(pattern) for pattern in SUPPORTED_PATTERNS)


def _scan_tree(directory):
    directory = os.path.abspath(directory)
    if not os.path.isdir(directory):
        raise FileNotFoundError("Directory does not exist.")
    node = _tree_directory_node(Path(directory))
    return node or {
        "type": "directory",
        "name": Path(directory).name or Path(directory).anchor,
        "path": str(Path(directory).resolve()),
        "children": [],
        "file_count": 0,
    }


def _tree_directory_node(directory):
    children = []
    file_count = 0
    entries = sorted(
        Path(directory).iterdir(),
        key=lambda path: (0 if path.is_dir() else 1, path.name.lower()),
    )
    for item in entries:
        if item.name.startswith(".") or item.name in IGNORED_DIRECTORIES:
            continue
        if item.is_dir():
            child = _tree_directory_node(item)
            if child is None:
                continue
            children.append(child)
            file_count += int(child.get("file_count", 0))
            continue
        if not item.is_file():
            continue
        info = _supported_file_info(item)
        if info is None:
            continue
        children.append(info)
        file_count += 1
    if not children:
        return None
    return {
        "type": "directory",
        "name": directory.name or directory.anchor,
        "path": str(directory.resolve()),
        "children": children,
        "file_count": file_count,
    }


def _flatten_file_nodes(node, files=None):
    files = [] if files is None else files
    if not node:
        return files
    if node.get("type") == "file":
        files.append(node)
        return files
    for child in node.get("children", []):
        _flatten_file_nodes(child, files)
    return files


def _extract_station_code(path):
    match = _STATION_CODE_RE.search(Path(path).name)
    if match is None:
        return "UNKNOWN"
    return match.group(1).upper()


def _safe_parse_scan_time(path):
    try:
        return parse_radar_time_from_filename(path)
    except Exception:
        return None


def _build_catalog(tree):
    stations = {}
    for item in _flatten_file_nodes(tree):
        station_id = _extract_station_code(item["path"])
        scan_time = _safe_parse_scan_time(item["path"])
        station = stations.setdefault(
            station_id,
            {
                "station_id": station_id,
                "station_name": station_id,
                "files": [],
            },
        )
        station["files"].append(
            {
                "name": item["name"],
                "path": item["path"],
                "format": item.get("format"),
                "scan_time": scan_time.isoformat() if scan_time is not None else None,
            }
        )
    ordered_stations = []
    for station_id in sorted(stations):
        files = sorted(
            stations[station_id]["files"],
            key=lambda item: (item["scan_time"] is None, item["scan_time"] or "", item["name"]),
        )
        latest_time = next((item["scan_time"] for item in reversed(files) if item["scan_time"] is not None), None)
        ordered_stations.append(
            {
                "station_id": station_id,
                "station_name": stations[station_id]["station_name"],
                "file_count": len(files),
                "latest_time": latest_time,
                "files": files,
            }
        )
    return {"stations": ordered_stations}


def _available_fields(radar):
    names = []
    for sweep in radar.fields:
        names.extend(list(sweep.data_vars))
    return _ordered_fields(names)


def _ordered_fields(names):
    seen = set()
    ordered = []
    for preferred in PREFERRED_FIELDS:
        if preferred in names and preferred not in seen:
            ordered.append(preferred)
            seen.add(preferred)
    for name in names:
        if name not in seen:
            ordered.append(name)
            seen.add(name)
    return ordered


def _field_sweeps(radar, fields):
    support = {}
    for field_name in fields:
        support[field_name] = [
            int(index)
            for index, sweep in enumerate(radar.fields)
            if field_name in sweep.data_vars
        ]
    return support


def _field_legends(fields):
    legends = {}
    if "HCL" in fields:
        legends["HCL"] = [
            {"value": index + 1, "label_zh": label}
            for index, label in enumerate(HYDROMETEOR_CLASS_NAMES_ZH)
        ]
    return legends


def _native_field_support(radar, fields):
    support = {}
    for field_name in fields:
        if not hasattr(radar, "has_extended_field"):
            support[field_name] = []
            continue
        support[field_name] = [int(sweep) for sweep in range(int(radar.nsweeps)) if radar.has_extended_field(sweep, field_name)]
    return support


def _field_range_km(field):
    ranges = np.asarray(field["range"].values, dtype=float)
    if ranges.size == 0:
        return 0.0
    return round(float(ranges[-1]) / 1000.0, 3)


def _sweep_profiles(radar):
    profiles = []
    for sweep_index, sweep in enumerate(radar.fields):
        field_names = _ordered_fields(list(sweep.data_vars))
        field_ranges = {}
        for field_name in field_names:
            aligned_range = _field_range_km(sweep[field_name])
            native_range = aligned_range
            if hasattr(radar, "has_extended_field") and hasattr(radar, "get_sweep_field") and radar.has_extended_field(sweep_index, field_name):
                native_range = _field_range_km(radar.get_sweep_field(sweep_index, field_name, range_mode="native"))
            field_ranges[field_name] = {
                "aligned": aligned_range,
                "native": native_range,
            }
        profiles.append(
            {
                "sweep": int(sweep_index),
                "fixed_angle": float(radar.scan_info.fixed_angle.values[sweep_index]),
                "fields": field_names,
                "field_ranges_km": field_ranges,
                "aligned_max_range_km": max((item["aligned"] for item in field_ranges.values()), default=0.0),
                "native_max_range_km": max((item["native"] for item in field_ranges.values()), default=0.0),
            }
        )
    return profiles


def _parse_bool(name, default=False):
    raw = request.args.get(name)
    if raw is None:
        return default
    return raw.lower() in ("1", "true", "yes", "on")


def _parse_int(name, default, minimum=None, maximum=None):
    value = request.args.get(name, default)
    value = int(value)
    if minimum is not None:
        value = max(minimum, value)
    if maximum is not None:
        value = min(maximum, value)
    return value


def _parse_range_mode(default="native"):
    value = request.args.get("range_mode", default).strip().lower()
    return "native" if value == "native" else "aligned"


def _parse_optional_float(name, minimum=None, maximum=None):
    raw = request.args.get(name)
    if raw is None or raw == "" or raw.lower() == "auto":
        return None
    value = float(raw)
    if minimum is not None:
        value = max(minimum, value)
    if maximum is not None:
        value = min(maximum, value)
    return value


def _parse_float(name, minimum=None, maximum=None):
    raw = request.args.get(name)
    if raw is None or raw == "":
        raise ValueError("Missing `%s` parameter." % name)
    value = float(raw)
    if minimum is not None:
        value = max(minimum, value)
    if maximum is not None:
        value = min(maximum, value)
    return value


def _parse_preset(default="standard"):
    value = request.args.get("preset", default).strip().lower()
    if value not in ("standard", "fine_detail", "weak_echo"):
        return default
    return value


def _clip_field_range(field, max_range_km):
    if max_range_km is None or "range" not in field.coords:
        return field
    ranges = np.asarray(field["range"].values, dtype=float)
    if ranges.size == 0:
        return field
    limit_m = float(max_range_km) * 1000.0
    mask = ranges <= (limit_m + 1e-6)
    if not np.any(mask):
        return field.isel(range=slice(0, 1))
    return field.isel(range=slice(0, int(mask.sum())))


def _ring_spacing_km(max_range_km):
    if max_range_km is None:
        return 25.0
    if max_range_km <= 60.0:
        return 10.0
    if max_range_km <= 120.0:
        return 20.0
    if max_range_km <= 260.0:
        return 25.0
    return 50.0


def _map_tick_step(max_range_km):
    if max_range_km is None:
        return 0.5
    if max_range_km <= 80.0:
        return 0.2
    if max_range_km <= 180.0:
        return 0.5
    return 1.0


def _theme_axes(fig, ax, artist):
    background = SPECIAL_COLORS["Background1"]
    axis_color = SPECIAL_COLORS["Axis"]
    text_color = SPECIAL_COLORS["Text1"]
    fig.patch.set_facecolor(background)
    ax.set_facecolor(background)
    ax.tick_params(colors=axis_color)
    ax.xaxis.label.set_color(text_color)
    ax.yaxis.label.set_color(text_color)
    ax.title.set_color(text_color)
    for spine in getattr(ax, "spines", {}).values():
        spine.set_color(axis_color)
    colorbar = getattr(artist, "pycwr_colorbar", None)
    if colorbar is not None:
        colorbar.outline.set_edgecolor(axis_color)
        colorbar.ax.tick_params(color=axis_color, labelcolor=text_color)
        colorbar.ax.yaxis.label.set_color(text_color)
        colorbar.ax.xaxis.label.set_color(text_color)
        for spine in colorbar.ax.spines.values():
            spine.set_edgecolor(axis_color)


def _error_png(message, width=900, height=800, status=400):
    dpi = 100
    fig = plt.figure(figsize=(width / dpi, height / dpi), dpi=dpi)
    ax = fig.add_subplot(111)
    fig.patch.set_facecolor(SPECIAL_COLORS["Background1"])
    ax.set_facecolor(SPECIAL_COLORS["Background1"])
    ax.axis("off")
    ax.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        color=SPECIAL_COLORS["Text1"],
        fontsize=12,
        wrap=True,
    )
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    buffer.seek(0)
    return send_file(buffer, mimetype="image/png"), status


def _log_and_mask_error(app, exc, message, status):
    app.logger.exception("pycwr web viewer request failed")
    return jsonify({"ok": False, "error": message}), status


def _render_png(radar, field_name, sweep, mapped, continuous, width, height, range_mode, max_range_km, preset):
    field, field_key = resolve_field_data(radar, sweep, field_name, range_mode=range_mode)
    field = _clip_field_range(field, max_range_km)
    style = build_web_style(
        radar,
        sweep,
        field_key,
        field,
        continuous=continuous,
        preset=preset,
        clabel=None,
        colorbar_visible=True,
        colorbar_orientation="vertical",
    )
    dpi = 100
    fig = plt.figure(figsize=(max(width, 320) / dpi, max(height, 320) / dpi), dpi=dpi)
    if mapped:
        if ccrs is None:
            raise RuntimeError("Map plotting requires cartopy.")
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
        artist = plot_map(
            ax,
            field.lon,
            field.lat,
            field,
            style,
            fig=fig,
            map_options=MapOptions(
                data_crs=ccrs.PlateCarree(),
                tick_step_degrees=_map_tick_step(max_range_km),
                province_linewidth=0.8,
                province_alpha=1.0,
            ),
        )
    else:
        ax = fig.add_subplot(111)
        artist = plot_cartesian(
            ax,
            field.x,
            field.y,
            field,
            style,
            fig=fig,
            title=default_sweep_title(radar, sweep),
            reference_options=CartesianReferenceOptions(
                ring_color=SPECIAL_COLORS["Ring"],
                spoke_color=SPECIAL_COLORS["Ring"],
                ring_spacing_km=_ring_spacing_km(max_range_km),
            ),
        )
    if not mapped:
        ax.set_title(default_sweep_title(radar, sweep))
    _theme_axes(fig, ax, artist)
    fig.tight_layout()
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    buffer.seek(0)
    return send_file(buffer, mimetype="image/png")


def _resolve_style_reference(radar, field_name, range_mode):
    for sweep_index in range(int(radar.nsweeps)):
        try:
            field_data, field_key = resolve_field_data(radar, sweep_index, field_name, range_mode=range_mode)
            return sweep_index, field_key, field_data
        except KeyError:
            continue
    raise KeyError("field %s not found in any sweep" % field_name)


def _section_height_limit(section, field_name):
    if field_name not in section:
        return None
    values = np.asarray(section[field_name].values, dtype=float)
    heights = np.asarray(section["z"].values, dtype=float)
    if values.size == 0 or heights.size == 0:
        return None
    valid_heights = heights[np.isfinite(values)]
    if valid_heights.size == 0:
        return None
    top_km = float(np.nanpercentile(valid_heights, 99.5) / 1000.0)
    if not np.isfinite(top_km) or top_km <= 0.0:
        return None
    upper_km = min(18.0, max(3.0, top_km + max(0.5, top_km * 0.08)))
    return (0.0, upper_km)


def _render_section_png(radar, field_name, start_xy_km, end_xy_km, width, height, range_mode, preset):
    start_xy_km = np.asarray(start_xy_km, dtype=np.float64)
    end_xy_km = np.asarray(end_xy_km, dtype=np.float64)
    if start_xy_km.shape != (2,) or end_xy_km.shape != (2,):
        raise ValueError("Section endpoints must contain exactly two coordinates.")
    if np.allclose(start_xy_km, end_xy_km, rtol=0.0, atol=1e-6):
        raise ValueError("Section start and end points must be different.")

    section = radar.extract_section(
        tuple(float(value) for value in start_xy_km),
        tuple(float(value) for value in end_xy_km),
        field_name=field_name,
        point_units="km",
        range_mode=range_mode,
    )
    sweep_index, field_key, field_data = _resolve_style_reference(radar, field_name, range_mode=range_mode)
    style = build_web_style(
        radar,
        sweep_index,
        field_key,
        field_data,
        continuous=False,
        preset=preset,
        clabel=None,
        colorbar_visible=True,
        colorbar_orientation="vertical",
    )

    dpi = 100
    fig = plt.figure(figsize=(max(width, 420) / dpi, max(height, 240) / dpi), dpi=dpi)
    ax = fig.add_subplot(111)
    height_limit = _section_height_limit(section, field_name)
    artist = plot_vertical_section(
        ax,
        section,
        None,
        None,
        style,
        fig=fig,
        height_km=height_limit,
        title="%s Vertical Section" % field_key,
    )
    ax.set_title("%s Vertical Section" % field_key, fontsize=12, pad=16)
    length_km = float(np.hypot(*(end_xy_km - start_xy_km)))
    ax.text(
        0.01,
        1.005,
        "Start (%.1f, %.1f) km | End (%.1f, %.1f) km | Length %.1f km"
        % (start_xy_km[0], start_xy_km[1], end_xy_km[0], end_xy_km[1], length_km),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        color=SPECIAL_COLORS["Text1"],
        fontsize=8.5,
    )
    _theme_axes(fig, ax, artist)
    fig.tight_layout()
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    buffer.seek(0)
    return send_file(buffer, mimetype="image/png")


def create_app(allowed_roots=None, auth_token=None):
    app = Flask(__name__, template_folder="templates", static_folder="static")
    cache = RadarFileCache()
    configured_roots = _normalize_allowed_roots(allowed_roots)
    viewer_token = auth_token or secrets.token_urlsafe(24)
    allowed_roots_label = ", ".join(str(root) for root in configured_roots) if configured_roots else "local filesystem"

    @app.get("/")
    def index():
        return render_template(
            "viewer.html",
            default_directory=_default_directory(),
            special_colors=SPECIAL_COLORS,
            api_token=viewer_token,
            allowed_roots=allowed_roots_label,
        )

    @app.before_request
    def _require_viewer_token():
        if request.endpoint == "index" or request.path.startswith("/static/"):
            return None
        supplied = request.args.get("token") or request.headers.get("X-Pycwr-Token")
        if supplied != viewer_token:
            return jsonify({"ok": False, "error": "Forbidden."}), 403
        return None

    @app.get("/api/files")
    def api_files():
        directory = request.args.get("dir", "").strip() or _default_directory()
        try:
            safe_directory = _resolve_within_roots(directory, configured_roots, expect_directory=True)
            files = _scan_directory(safe_directory)
        except PermissionError:
            return jsonify({"ok": False, "error": "Requested directory is outside the allowed roots."}), 403
        except Exception as exc:
            return _log_and_mask_error(app, exc, "Unable to scan directory.", 400)
        return jsonify({"ok": True, "directory": safe_directory, "files": files})

    @app.get("/api/tree")
    def api_tree():
        directory = request.args.get("dir", "").strip() or _default_directory()
        try:
            safe_directory = _resolve_within_roots(directory, configured_roots, expect_directory=True)
            tree = _scan_tree(safe_directory)
        except PermissionError:
            return jsonify({"ok": False, "error": "Requested directory is outside the allowed roots."}), 403
        except Exception as exc:
            return _log_and_mask_error(app, exc, "Unable to build directory tree.", 400)
        return jsonify({"ok": True, "directory": safe_directory, "tree": tree})

    @app.get("/api/catalog")
    def api_catalog():
        directory = request.args.get("dir", "").strip() or _default_directory()
        try:
            safe_directory = _resolve_within_roots(directory, configured_roots, expect_directory=True)
            tree = _scan_tree(safe_directory)
            catalog = _build_catalog(tree)
        except PermissionError:
            return jsonify({"ok": False, "error": "Requested directory is outside the allowed roots."}), 403
        except Exception as exc:
            return _log_and_mask_error(app, exc, "Unable to build radar catalog.", 400)
        return jsonify({"ok": True, "directory": safe_directory, "catalog": catalog})

    @app.get("/api/metadata")
    def api_metadata():
        file_path = request.args.get("path", "").strip()
        if not file_path:
            return jsonify({"ok": False, "error": "Missing `path` parameter."}), 400
        try:
            safe_file_path = _resolve_within_roots(file_path, configured_roots, expect_directory=False)
            if radar_format(safe_file_path) is None:
                return jsonify({"ok": False, "error": "Unsupported radar file."}), 400
            radar = cache.get(safe_file_path)
        except PermissionError:
            return jsonify({"ok": False, "error": "Requested file is outside the allowed roots."}), 403
        except FileNotFoundError:
            return jsonify({"ok": False, "error": "File does not exist."}), 404
        except ValueError as exc:
            return jsonify({"ok": False, "error": str(exc) or "Unable to read the selected radar file."}), 400
        except Exception as exc:
            return _log_and_mask_error(app, exc, "Unable to read radar metadata.", 500)
        fields = _available_fields(radar)
        sweeps = list(range(int(radar.nsweeps)))
        return jsonify(
            {
                "ok": True,
                "path": os.path.abspath(safe_file_path),
                "name": os.path.basename(safe_file_path),
                "fields": fields,
                "field_legends": _field_legends(fields),
                "field_sweeps": _field_sweeps(radar, fields),
                "sweeps": sweeps,
                "native_support": _native_field_support(radar, fields),
                "sweep_profiles": _sweep_profiles(radar),
                "fixed_angles": [float(value) for value in radar.scan_info.fixed_angle.values],
                "scan_type": str(radar.scan_info.scan_type.values),
                "start_time": str(radar.scan_info.start_time.values),
                "site": {
                    "longitude": float(radar.scan_info.longitude.values),
                    "latitude": float(radar.scan_info.latitude.values),
                    "altitude": float(radar.scan_info.altitude.values),
                },
            }
        )

    def plot_response(mapped):
        file_path = request.args.get("path", "").strip()
        if not file_path:
            return _error_png("Missing `path` parameter.", status=400)
        try:
            safe_file_path = _resolve_within_roots(file_path, configured_roots, expect_directory=False)
            radar = cache.get(safe_file_path)
            available_fields = _available_fields(radar)
            field_name = request.args.get("field", available_fields[0] if available_fields else "dBZ")
            sweep = _parse_int("sweep", 0, minimum=0, maximum=max(int(radar.nsweeps) - 1, 0))
            continuous = _parse_bool("continuous", default=False)
            range_mode = _parse_range_mode(default="native")
            max_range_km = _parse_optional_float("max_range_km", minimum=5.0, maximum=600.0)
            preset = _parse_preset(default="standard")
            width = _parse_int("width", 960, minimum=320, maximum=2200)
            height = _parse_int("height", 960 if not mapped else 820, minimum=320, maximum=2200)
            return _render_png(radar, field_name, sweep, mapped, continuous, width, height, range_mode, max_range_km, preset)
        except PermissionError:
            return _error_png("Requested file is outside the allowed roots.", status=403)
        except FileNotFoundError:
            return _error_png("File does not exist.", status=404)
        except ValueError as exc:
            return _error_png(str(exc) or "Unable to read the selected radar file.", status=400)
        except Exception:
            app.logger.exception("pycwr web viewer plot request failed")
            return _error_png("Unable to render radar plot.", status=500)

    @app.get("/plot/ppi.png")
    def plot_ppi():
        return plot_response(mapped=False)

    @app.get("/plot/ppi_map.png")
    def plot_ppi_map():
        return plot_response(mapped=True)

    @app.get("/plot/section.png")
    def plot_section_png():
        file_path = request.args.get("path", "").strip()
        if not file_path:
            return _error_png("Missing `path` parameter.", status=400)
        try:
            safe_file_path = _resolve_within_roots(file_path, configured_roots, expect_directory=False)
            radar = cache.get(safe_file_path)
            available_fields = _available_fields(radar)
            field_name = request.args.get("field", available_fields[0] if available_fields else "dBZ")
            range_mode = _parse_range_mode(default="native")
            preset = _parse_preset(default="standard")
            width = _parse_int("width", 1200, minimum=420, maximum=2200)
            height = _parse_int("height", 520, minimum=240, maximum=1400)
            start_x_km = _parse_float("start_x_km", minimum=-600.0, maximum=600.0)
            start_y_km = _parse_float("start_y_km", minimum=-600.0, maximum=600.0)
            end_x_km = _parse_float("end_x_km", minimum=-600.0, maximum=600.0)
            end_y_km = _parse_float("end_y_km", minimum=-600.0, maximum=600.0)
            return _render_section_png(
                radar,
                field_name=field_name,
                start_xy_km=(start_x_km, start_y_km),
                end_xy_km=(end_x_km, end_y_km),
                width=width,
                height=height,
                range_mode=range_mode,
                preset=preset,
            )
        except PermissionError:
            return _error_png("Requested file is outside the allowed roots.", status=403)
        except FileNotFoundError:
            return _error_png("File does not exist.", status=404)
        except ValueError as exc:
            return _error_png(str(exc) or "Unable to render the selected section.", status=400)
        except Exception:
            app.logger.exception("pycwr web viewer section request failed")
            return _error_png("Unable to render vertical section.", status=500)

    return app


def launch(host="127.0.0.1", port=8787, open_browser=True):
    normalized_host = str(host).strip()
    if normalized_host not in ("127.0.0.1", "localhost"):
        raise ValueError("The built-in web viewer only supports loopback hosts.")
    app = create_app()
    url = "http://%s:%s/" % (host, port)
    if open_browser:
        threading.Timer(0.8, lambda: webbrowser.open(url)).start()
    print("pycwr web viewer running at %s" % url)
    app.run(host=host, port=port, debug=False, use_reloader=False)


if __name__ == "__main__":
    launch()
