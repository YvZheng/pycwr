from dataclasses import dataclass, field
from importlib import import_module
import math
from typing import Any, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize
from matplotlib.ticker import MaxNLocator, MultipleLocator

from ..configure.default_config import (
    CINRAD_COLORMAP,
    CINRAD_field_bins,
    CINRAD_field_mapping,
    CINRAD_field_normvar,
    DEFAULT_METADATA,
    HYDROMETEOR_CLASS_COLORS,
    HYDROMETEOR_CLASS_SHORT_LABELS_ZH,
)
from ..core.transforms import (
    antenna_vectors_to_cartesian_cwr,
    cartesian_to_geographic_aeqd,
    geographic_to_cartesian_aeqd,
)

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.mpl.geoaxes import GeoAxes
    from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

    from ..configure.location_config import get_cn_shp_info
except ImportError:
    ccrs = None
    cfeature = None
    GeoAxes = None
    get_cn_shp_info = None
    LatitudeFormatter = None
    LongitudeFormatter = None


@dataclass
class ColorbarOptions:
    visible: bool = True
    orientation: str = "vertical"
    label: Optional[str] = None
    ticks: Optional[Sequence[float]] = None
    ticklabels: Optional[Sequence[str]] = None
    cax: Optional[Axes] = None


@dataclass
class PlotStyle:
    cmap: Any
    value_range: Optional[Tuple[float, float]] = None
    bins: Optional[int] = None
    levels: Optional[Sequence[float]] = None
    norm: Any = None
    continuous: bool = False
    colorbar: ColorbarOptions = field(default_factory=ColorbarOptions)


@dataclass
class CartesianReferenceOptions:
    visible: bool = True
    ring_spacing_km: Optional[float] = None
    ring_color: str = "#6A6A6A"
    ring_linewidth: float = 0.6
    ring_linestyle: str = "-"
    spoke_count: int = 12
    spoke_color: str = "#6A6A6A"
    spoke_linewidth: float = 0.6
    spoke_linestyle: str = "-"


@dataclass
class MapOptions:
    data_crs: Any = None
    feature_scale: str = "50m"
    tick_step_degrees: Optional[float] = None
    add_ocean: bool = True
    add_land: bool = True
    add_lakes: bool = True
    add_rivers: bool = True
    add_province_outline: bool = True
    province_linewidth: float = 0.5
    province_alpha: float = 0.8


def require_mpl_axes(ax):
    if not isinstance(ax, Axes):
        raise TypeError("ax should be a matplotlib Axes instance")
    if GeoAxes is not None and isinstance(ax, GeoAxes):
        raise TypeError("ax should be a Cartesian matplotlib Axes, not a cartopy GeoAxes")


def require_cartopy_axes(ax):
    if GeoAxes is None or cfeature is None or get_cn_shp_info is None or ccrs is None:
        raise ImportError(
            "Map plotting requires the full optional stack. "
            "Install it with `pip install \"pycwr[full]\"`."
        )
    get_cn_shp_info()
    if not isinstance(ax, GeoAxes):
        raise TypeError("ax should be a cartopy GeoAxes instance")


def as_numpy(value):
    if value is None:
        return None
    if hasattr(value, "values"):
        value = value.values
    return np.asarray(value)


def as_scalar(value):
    array = as_numpy(value)
    if array is None:
        return None
    if array.shape == ():
        return array.item()
    return array


def ensure_ppi_scan(radar):
    scan_type = as_scalar(radar.scan_info.scan_type)
    if scan_type != "ppi":
        raise ValueError("Only ppi scan plotting is supported by this interface.")


def resolve_field_data(radar, sweep_index, field_name, range_mode=None):
    sweep = radar.fields[sweep_index]
    if field_name in sweep.data_vars:
        field_key = field_name
    else:
        for candidate, mapped_name in CINRAD_field_mapping.items():
            if mapped_name == field_name and candidate in sweep.data_vars:
                field_key = candidate
                break
        else:
            raise KeyError("field %s not found in sweep %s" % (field_name, sweep_index))
    if hasattr(radar, "get_sweep_field"):
        if hasattr(radar, "_resolve_field_range_mode"):
            range_mode = radar._resolve_field_range_mode(field_key, range_mode=range_mode)
        return radar.get_sweep_field(sweep_index, field_key, range_mode=range_mode), field_key
    return sweep[field_key], field_key


def resolve_field_metadata_name(field_key):
    return CINRAD_field_mapping.get(field_key, field_key)


def default_colorbar_label(field_key):
    metadata_name = resolve_field_metadata_name(field_key)
    if field_key == "HCL" or metadata_name == "hydro_class":
        return "Hydrometeor Class / 水凝物分类"
    metadata = DEFAULT_METADATA.get(metadata_name, {})
    units = metadata.get("units")
    if units:
        return "%s (%s)" % (metadata_name, units)
    return metadata_name


def default_sweep_title(radar, sweep_index):
    timestamp = pd.to_datetime(radar.fields[sweep_index].time[0].item()).strftime("UTC %Y-%m-%d %H:%M:%S")
    fixed_angle = float(as_scalar(radar.scan_info.fixed_angle[sweep_index]))
    return "%s Elevation: %.1f deg" % (timestamp, fixed_angle)


def resolve_field_range(radar, sweep_index, field_key, field_data, value_range=None):
    if value_range is not None:
        vmin, vmax = value_range
    else:
        metadata_name = resolve_field_metadata_name(field_key)
        if field_key == "V":
            vmax = float(as_scalar(radar.scan_info.nyquist_velocity[sweep_index]))
            vmin = -1.0 * vmax
        else:
            configured_range = CINRAD_field_normvar.get(metadata_name, -1)
            if configured_range == -1:
                data = as_numpy(field_data).astype(float)
                vmin = float(np.nanmin(data))
                vmax = float(np.nanmax(data))
            else:
                vmin, vmax = configured_range
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        return 0.0, 1.0
    if vmin == vmax:
        delta = max(abs(vmin) * 0.05, 1.0)
        return vmin - delta, vmax + delta
    return float(vmin), float(vmax)


def resolve_field_style(
    radar,
    sweep_index,
    field_key,
    field_data,
    value_range=None,
    cmap=None,
    bins=None,
    continuous=False,
    colorbar_visible=True,
    colorbar_orientation="vertical",
    colorbar_ticks=None,
    colorbar_ticklabels=None,
    colorbar_label=None,
    colorbar_cax=None,
):
    metadata_name = resolve_field_metadata_name(field_key)
    if field_key == "HCL" or metadata_name == "hydro_class":
        tick_values = list(range(1, len(HYDROMETEOR_CLASS_SHORT_LABELS_ZH) + 1))
        return PlotStyle(
            cmap=cmap if cmap is not None else ListedColormap(list(HYDROMETEOR_CLASS_COLORS), name="pycwr_hcl"),
            levels=np.arange(0.5, len(HYDROMETEOR_CLASS_SHORT_LABELS_ZH) + 1.5, 1.0),
            continuous=False,
            colorbar=ColorbarOptions(
                visible=colorbar_visible,
                orientation=colorbar_orientation,
                label=colorbar_label if colorbar_label is not None else default_colorbar_label(field_key),
                ticks=colorbar_ticks if colorbar_ticks is not None else tick_values,
                ticklabels=colorbar_ticklabels if colorbar_ticklabels is not None else list(HYDROMETEOR_CLASS_SHORT_LABELS_ZH),
                cax=colorbar_cax,
            ),
        )
    return PlotStyle(
        cmap=cmap if cmap is not None else CINRAD_COLORMAP.get(metadata_name, "viridis"),
        value_range=resolve_field_range(radar, sweep_index, field_key, field_data, value_range=value_range),
        bins=bins if bins is not None else CINRAD_field_bins.get(metadata_name, 16),
        continuous=continuous,
        colorbar=ColorbarOptions(
            visible=colorbar_visible,
            orientation=colorbar_orientation,
            label=colorbar_label if colorbar_label is not None else default_colorbar_label(field_key),
            ticks=colorbar_ticks,
            ticklabels=colorbar_ticklabels,
            cax=colorbar_cax,
        ),
    )


def normalize_style(style, data=None):
    if style.value_range is None:
        if data is None:
            style.value_range = (0.0, 1.0)
        else:
            array = as_numpy(data).astype(float)
            vmin = float(np.nanmin(array))
            vmax = float(np.nanmax(array))
            if vmin == vmax:
                delta = max(abs(vmin) * 0.05, 1.0)
                style.value_range = (vmin - delta, vmax + delta)
            else:
                style.value_range = (vmin, vmax)
    if style.bins is None:
        style.bins = 16
    return style


def _format_ticks(ticks):
    ticks = np.asarray(ticks)
    if np.allclose(ticks, np.round(ticks), equal_nan=False):
        return ["%2.f" % tick for tick in ticks]
    return ["%.2f" % tick for tick in ticks]


def _resolve_matplotlib_cmap(cmap_spec):
    try:
        return plt.get_cmap(cmap_spec)
    except ValueError:
        # Legacy plotting entry points may bypass pycwr.draw.__getattr__,
        # so ensure pycwr's named colormaps are registered before retrying.
        if isinstance(cmap_spec, str):
            import_module(".colormap", __package__)
            return plt.get_cmap(cmap_spec)
        raise


def build_colormap(style):
    style = normalize_style(style)
    cmap = _resolve_matplotlib_cmap(style.cmap)
    if style.levels is not None:
        levels = np.asarray(style.levels, dtype=float)
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        ticks = list(style.colorbar.ticks) if style.colorbar.ticks is not None else levels
        ticklabels = style.colorbar.ticklabels if style.colorbar.ticklabels is not None else _format_ticks(ticks)
        return cmap, norm, ticks, ticklabels
    if style.norm is not None:
        ticks = list(style.colorbar.ticks) if style.colorbar.ticks is not None else None
        ticklabels = style.colorbar.ticklabels
        return cmap, style.norm, ticks, ticklabels
    vmin, vmax = style.value_range
    if style.continuous:
        norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
        ticks = list(style.colorbar.ticks) if style.colorbar.ticks is not None else None
        ticklabels = style.colorbar.ticklabels
        return cmap, norm, ticks, ticklabels
    levels = MaxNLocator(nbins=style.bins).tick_values(vmin, vmax)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    ticks = list(style.colorbar.ticks) if style.colorbar.ticks is not None else levels
    ticklabels = style.colorbar.ticklabels if style.colorbar.ticklabels is not None else _format_ticks(ticks)
    return cmap, norm, ticks, ticklabels


def create_standard_axes(fig=None, orientation="vertical", projection=None):
    if fig is None:
        fig = plt.figure()
    if orientation == "horizontal":
        ax = fig.add_axes([0.08, 0.22, 0.84, 0.68], projection=projection)
        cax = fig.add_axes([0.08, 0.10, 0.84, 0.05])
    else:
        ax = fig.add_axes([0.08, 0.10, 0.76, 0.82], projection=projection)
        cax = fig.add_axes([0.86, 0.10, 0.03, 0.82])
    return fig, ax, cax


def _edges_from_centers_1d(values):
    values = as_numpy(values).astype(float)
    if values.ndim != 1:
        raise ValueError("values must be 1-D")
    if values.size == 1:
        delta = 0.5
        return np.array([values[0] - delta, values[0] + delta], dtype=float)
    edges = np.empty(values.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (values[:-1] + values[1:])
    edges[0] = values[0] - 0.5 * (values[1] - values[0])
    edges[-1] = values[-1] + 0.5 * (values[-1] - values[-2])
    return edges


def _linear_pad_2d(values):
    values = as_numpy(values).astype(float)
    if values.ndim != 2:
        raise ValueError("values must be 2-D")
    nrows, ncols = values.shape
    padded = np.empty((nrows + 2, ncols + 2), dtype=float)
    padded[1:-1, 1:-1] = values
    padded[0, 1:-1] = values[0, :] if nrows == 1 else 2.0 * values[0, :] - values[1, :]
    padded[-1, 1:-1] = values[-1, :] if nrows == 1 else 2.0 * values[-1, :] - values[-2, :]
    padded[1:-1, 0] = values[:, 0] if ncols == 1 else 2.0 * values[:, 0] - values[:, 1]
    padded[1:-1, -1] = values[:, -1] if ncols == 1 else 2.0 * values[:, -1] - values[:, -2]
    padded[0, 0] = padded[1, 0] + padded[0, 1] - padded[1, 1]
    padded[0, -1] = padded[1, -1] + padded[0, -2] - padded[1, -2]
    padded[-1, 0] = padded[-2, 0] + padded[-1, 1] - padded[-2, 1]
    padded[-1, -1] = padded[-2, -1] + padded[-1, -2] - padded[-2, -2]
    return padded


def _edges_from_centers_2d(values):
    padded = _linear_pad_2d(values)
    return 0.25 * (
        padded[:-1, :-1]
        + padded[1:, :-1]
        + padded[:-1, 1:]
        + padded[1:, 1:]
    )


def pcolormesh_coordinates(x_coords, y_coords, data):
    x_coords = as_numpy(x_coords)
    y_coords = as_numpy(y_coords)
    data = as_numpy(data)
    if data.ndim != 2:
        return x_coords, y_coords
    if x_coords.ndim == 1 and y_coords.ndim == 1:
        if x_coords.size == data.shape[0] and y_coords.size == data.shape[1]:
            return _edges_from_centers_1d(x_coords), _edges_from_centers_1d(y_coords)
        return x_coords, y_coords
    if x_coords.ndim == 2 and y_coords.ndim == 2 and x_coords.shape == data.shape and y_coords.shape == data.shape:
        return _edges_from_centers_2d(x_coords), _edges_from_centers_2d(y_coords)
    return x_coords, y_coords


def ensure_cartesian_coordinates(
    radar_data=None,
    x=None,
    y=None,
    ranges=None,
    azimuth=None,
    elevation=None,
    altitude_m=0.0,
    effective_earth_radius=None,
):
    if x is not None and y is not None:
        return as_numpy(x), as_numpy(y)
    if radar_data is not None and hasattr(radar_data, "x") and hasattr(radar_data, "y"):
        return as_numpy(radar_data.x), as_numpy(radar_data.y)
    if ranges is None or azimuth is None or elevation is None:
        raise ValueError("x/y or range/azimuth/elevation coordinates are required.")
    x, y, _ = antenna_vectors_to_cartesian_cwr(
        as_numpy(ranges),
        as_numpy(azimuth),
        as_numpy(elevation),
        h=altitude_m,
        effective_earth_radius=effective_earth_radius,
    )
    return x, y


def ensure_geographic_coordinates(
    radar_data=None,
    lon=None,
    lat=None,
    x=None,
    y=None,
    ranges=None,
    azimuth=None,
    elevation=None,
    station_lonlat=None,
    altitude_m=0.0,
    effective_earth_radius=None,
):
    if lon is not None and lat is not None:
        return as_numpy(lon), as_numpy(lat)
    if radar_data is not None and hasattr(radar_data, "lon") and hasattr(radar_data, "lat"):
        return as_numpy(radar_data.lon), as_numpy(radar_data.lat)
    if station_lonlat is None:
        raise ValueError("station_lonlat=(lon, lat) is required when lon/lat coordinates are unavailable.")
    x, y = ensure_cartesian_coordinates(
        radar_data=radar_data,
        x=x,
        y=y,
        ranges=ranges,
        azimuth=azimuth,
        elevation=elevation,
        altitude_m=altitude_m,
        effective_earth_radius=effective_earth_radius,
    )
    station_lon, station_lat = station_lonlat
    lon, lat = cartesian_to_geographic_aeqd(x, y, station_lon, station_lat)
    return lon, lat


def _max_cartesian_radius_km(x_m, y_m):
    return float(np.nanmax(np.hypot(as_numpy(x_m), as_numpy(y_m))) / 1000.0)


def cartesian_extent_from_xy(x_m, y_m, extent_km=None):
    if extent_km is not None:
        return extent_km
    max_radius_km = _max_cartesian_radius_km(x_m, y_m)
    return (-max_radius_km, max_radius_km, -max_radius_km, max_radius_km)


def geographic_extent_from_lonlat(lon, lat, extent=None):
    if extent is not None:
        return extent
    return (
        float(np.nanmin(as_numpy(lon))),
        float(np.nanmax(as_numpy(lon))),
        float(np.nanmin(as_numpy(lat))),
        float(np.nanmax(as_numpy(lat))),
    )


def normalize_geographic_extent(extent, min_span_degrees=0.1):
    lon_min, lon_max, lat_min, lat_max = [float(value) for value in extent]
    min_span = float(min_span_degrees)
    lon_span = lon_max - lon_min
    lat_span = lat_max - lat_min
    if lon_span < min_span:
        center = 0.5 * (lon_min + lon_max)
        half_span = 0.5 * min_span
        lon_min = center - half_span
        lon_max = center + half_span
    if lat_span < min_span:
        center = 0.5 * (lat_min + lat_max)
        half_span = 0.5 * min_span
        lat_min = center - half_span
        lat_max = center + half_span
    return lon_min, lon_max, lat_min, lat_max


def nice_step(max_value, target_steps=5):
    if max_value <= 0:
        return 1.0
    rough_step = float(max_value) / float(max(target_steps, 1))
    exponent = math.floor(math.log10(rough_step))
    base = 10 ** exponent
    for factor in (1.0, 2.0, 5.0, 10.0):
        step = factor * base
        if step >= rough_step:
            return step
    return 10.0 * base


def apply_cartesian_axes(ax, extent_km, title=None, labels=None):
    min_x, max_x, min_y, max_y = extent_km
    if labels is None:
        labels = ("East of radar (km)", "North of radar (km)")
    max_extent = max(abs(min_x), abs(max_x), abs(min_y), abs(max_y))
    major_step = nice_step(max_extent, target_steps=4)
    minor_step = max(major_step / 5.0, 1.0)
    ax.set_aspect("equal")
    ax.set_xlim([min_x, max_x])
    ax.set_ylim([min_y, max_y])
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_title(title or "")
    ax.xaxis.set_major_locator(MultipleLocator(major_step))
    ax.xaxis.set_minor_locator(MultipleLocator(minor_step))
    ax.yaxis.set_major_locator(MultipleLocator(major_step))
    ax.yaxis.set_minor_locator(MultipleLocator(minor_step))
    ax.tick_params(axis="both", which="major", direction="in")
    ax.tick_params(axis="both", which="minor", direction="in")


def add_cartesian_reference(ax, extent_km, reference_options=None):
    reference_options = reference_options or CartesianReferenceOptions()
    if not reference_options.visible:
        return
    max_extent = max(abs(extent_km[0]), abs(extent_km[1]), abs(extent_km[2]), abs(extent_km[3]))
    ring_spacing_km = reference_options.ring_spacing_km or nice_step(max_extent, target_steps=5)
    theta = np.linspace(0.0, 2.0 * np.pi, 256)
    rings = np.arange(ring_spacing_km, max_extent + 0.5 * ring_spacing_km, ring_spacing_km)
    for radius in rings:
        x0 = radius * np.cos(theta)
        y0 = radius * np.sin(theta)
        ax.plot(
            x0,
            y0,
            color=reference_options.ring_color,
            linewidth=reference_options.ring_linewidth,
            linestyle=reference_options.ring_linestyle,
            zorder=2,
        )
    if rings.size == 0:
        return
    spoke_angles = np.linspace(0.0, np.pi, reference_options.spoke_count, endpoint=False)
    max_radius = rings[-1]
    for angle in spoke_angles:
        ax.plot(
            [-1.0 * max_radius * np.sin(angle), max_radius * np.sin(angle)],
            [-1.0 * max_radius * np.cos(angle), max_radius * np.cos(angle)],
            color=reference_options.spoke_color,
            linewidth=reference_options.spoke_linewidth,
            linestyle=reference_options.spoke_linestyle,
            zorder=2,
        )


def add_colorbar(ax, mappable, style, fig=None):
    if not style.colorbar.visible:
        return None
    fig = fig if fig is not None else ax.figure
    cmap, _, ticks, ticklabels = build_colormap(style)
    del cmap
    colorbar = fig.colorbar(
        mappable=mappable,
        ax=ax,
        cax=style.colorbar.cax,
        orientation=style.colorbar.orientation,
        ticks=ticks,
    )
    if ticklabels is not None and ticks is not None:
        if style.colorbar.orientation == "vertical":
            colorbar.ax.set_yticklabels(ticklabels)
        else:
            colorbar.ax.set_xticklabels(ticklabels)
    if style.colorbar.label:
        colorbar.set_label(style.colorbar.label)
    mappable.pycwr_colorbar = colorbar
    return colorbar


def plot_cartesian(
    ax,
    x_m,
    y_m,
    data,
    style,
    fig=None,
    extent_km=None,
    title=None,
    labels=None,
    reference_options=None,
    **kwargs
):
    require_mpl_axes(ax)
    cmap, norm, _, _ = build_colormap(style)
    mesh_x, mesh_y = pcolormesh_coordinates(x_m / 1000.0, y_m / 1000.0, data)
    mesh = ax.pcolormesh(
        mesh_x,
        mesh_y,
        as_numpy(data),
        cmap=cmap,
        norm=norm,
        shading="auto",
        **kwargs
    )
    extent_km = cartesian_extent_from_xy(x_m, y_m, extent_km=extent_km)
    add_cartesian_reference(ax, extent_km, reference_options=reference_options)
    apply_cartesian_axes(ax, extent_km, title=title, labels=labels)
    add_colorbar(ax, mesh, style, fig=fig)
    ax.pycwr_last_mappable = mesh
    return mesh


def _default_map_options(map_options):
    map_options = map_options or MapOptions()
    if map_options.data_crs is None:
        map_options.data_crs = ccrs.PlateCarree() if ccrs is not None else None
    return map_options


def _auto_tick_step(extent):
    lon_span = max(abs(extent[1] - extent[0]), 0.1)
    lat_span = max(abs(extent[3] - extent[2]), 0.1)
    return max(nice_step(max(lon_span, lat_span), target_steps=5), 0.1)


def apply_map_axes(ax, extent, map_options, title=None):
    map_options = _default_map_options(map_options)
    data_crs = map_options.data_crs
    extent = normalize_geographic_extent(extent)
    ax.set_extent([extent[0], extent[1], extent[2], extent[3]], crs=data_crs)
    if map_options.add_ocean:
        ax.add_feature(cfeature.OCEAN.with_scale(map_options.feature_scale), zorder=0)
    if map_options.add_land:
        ax.add_feature(
            cfeature.NaturalEarthFeature(
                "physical",
                "land",
                map_options.feature_scale,
                edgecolor="none",
                facecolor="white",
            ),
            zorder=1,
        )
    if map_options.add_lakes:
        ax.add_feature(cfeature.LAKES.with_scale(map_options.feature_scale), zorder=2)
    if map_options.add_rivers:
        ax.add_feature(cfeature.RIVERS.with_scale(map_options.feature_scale), zorder=3)
    if map_options.add_province_outline:
        shp_info = get_cn_shp_info()
        ax.add_feature(
            cfeature.ShapelyFeature(
                shp_info.geometries(),
                data_crs,
                edgecolor="k",
                facecolor="none",
            ),
            linewidth=map_options.province_linewidth,
            linestyle="-",
            zorder=5,
            alpha=map_options.province_alpha,
        )
    ax.set_title(title or "", fontsize=14)
    tick_step = map_options.tick_step_degrees or _auto_tick_step(extent)
    meridians = np.arange(math.floor(extent[0]), math.ceil(extent[1]) + tick_step, tick_step)
    parallels = np.arange(math.floor(extent[2]), math.ceil(extent[3]) + tick_step, tick_step)
    ax.set_xticks(meridians, crs=data_crs)
    ax.set_yticks(parallels, crs=data_crs)
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())


def plot_map(
    ax,
    lon,
    lat,
    data,
    style,
    fig=None,
    extent=None,
    map_options=None,
    **kwargs
):
    require_cartopy_axes(ax)
    map_options = _default_map_options(map_options)
    cmap, norm, _, _ = build_colormap(style)
    extent = geographic_extent_from_lonlat(lon, lat, extent=extent)
    apply_map_axes(ax, extent, map_options=map_options)
    mesh_lon, mesh_lat = pcolormesh_coordinates(lon, lat, data)
    mesh = ax.pcolormesh(
        mesh_lon,
        mesh_lat,
        as_numpy(data),
        cmap=cmap,
        norm=norm,
        transform=map_options.data_crs,
        shading="auto",
        zorder=20,
        **kwargs
    )
    add_colorbar(ax, mesh, style, fig=fig)
    ax.pycwr_last_mappable = mesh
    return mesh


def _trim_section_data(mesh_x, mesh_z, section_data):
    mesh_shape = as_numpy(mesh_x).shape
    array = as_numpy(section_data)
    if array.ndim != 2:
        raise ValueError("section data must be 2-D")
    candidates = [array]
    if array.shape[1] > 1:
        candidates.append(array[:, :-1])
    if array.shape[0] > 1:
        candidates.append(array[:-1, :])
    if array.shape[0] > 1 and array.shape[1] > 1:
        candidates.append(array[:-1, :-1])
    for candidate in candidates:
        if candidate.shape == mesh_shape:
            return candidate
        if candidate.shape[0] + 1 == mesh_shape[0] and candidate.shape[1] + 1 == mesh_shape[1]:
            return candidate
    raise ValueError("section mesh and field shapes do not align")


def _sanitize_section_panel(x_coords, z_coords, section_data):
    """Drop non-finite padded rows/columns before plotting a section panel."""
    x_array = as_numpy(x_coords).astype(float)
    z_array = as_numpy(z_coords).astype(float)
    section_array = as_numpy(section_data).astype(float)
    if x_array.shape != section_array.shape or z_array.shape != section_array.shape:
        return x_array, z_array, section_array
    finite = np.isfinite(x_array) & np.isfinite(z_array)
    if not np.any(finite):
        return x_array[:0, :0], z_array[:0, :0], section_array[:0, :0]
    valid_rows = np.any(finite, axis=1)
    valid_cols = np.all(finite[valid_rows, :], axis=0)
    if not np.any(valid_cols):
        valid_cols = np.any(finite[valid_rows, :], axis=0)
    x_array = x_array[valid_rows][:, valid_cols]
    z_array = z_array[valid_rows][:, valid_cols]
    section_array = section_array[valid_rows][:, valid_cols]
    return x_array, z_array, section_array


def section_dataset_to_arrays(section, field_name=None):
    """Normalize a section dataset into plotting arrays."""
    if field_name is None:
        field_name = section.attrs.get("field_name")
    if not field_name or field_name not in section.data_vars:
        raise KeyError("field_name must reference a data variable in the section dataset")
    if "distance_ground" in section:
        mesh_x = as_numpy(section["distance_ground"])
    elif "distance" in section.coords:
        mesh_x = np.broadcast_to(as_numpy(section["distance"]), as_numpy(section[field_name]).shape)
    else:
        raise KeyError("section dataset must provide distance information")
    return mesh_x, as_numpy(section["z"]), as_numpy(section[field_name]), field_name


def _section_panels(mesh_xy, mesh_z, field_data):
    if hasattr(mesh_xy, "data_vars"):
        dataset_mesh_x, dataset_mesh_z, dataset_field, _ = section_dataset_to_arrays(mesh_xy, field_name=mesh_z)
        return [(dataset_mesh_x, dataset_mesh_z, dataset_field)]
    if isinstance(mesh_xy, (list, tuple)):
        return list(zip(mesh_xy, mesh_z, field_data))
    return [(mesh_xy, mesh_z, field_data)]


def plot_vertical_section(
    ax,
    mesh_xy,
    mesh_z,
    field_data,
    style,
    fig=None,
    height_km=(0.0, 18.0),
    title=None,
    labels=None,
    **kwargs
):
    require_mpl_axes(ax)
    if labels is None:
        labels = ("Distance from section start (km)", "Height (km)")
    cmap, norm, _, _ = build_colormap(style)
    last_mesh = None
    height_limit = height_km
    panels = _section_panels(mesh_xy, mesh_z, field_data)
    if not panels:
        raise ValueError("field_data is empty")
    if height_limit is None:
        max_height = 0.0
        for _, z_coords, section in panels:
            array = as_numpy(section).astype(float)
            if array.size == 0 or np.isnan(array).all():
                continue
            max_height = max(max_height, float(np.nanmax(as_numpy(z_coords)) / 1000.0))
        headroom = 0.5 if max_height <= 0.0 else max(0.5, max_height * 0.05)
        height_limit = (0.0, max_height + headroom)
    for x_coords, z_coords, section in panels:
        array = as_numpy(section).astype(float)
        if array.size == 0 or np.isnan(array).all():
            continue
        panel_x, panel_z, panel_section = _sanitize_section_panel(x_coords, z_coords, section)
        if panel_section.size == 0 or np.isnan(panel_section).all():
            continue
        trimmed = _trim_section_data(panel_x, panel_z, panel_section)
        mesh_x, mesh_z_coords = pcolormesh_coordinates(panel_x / 1000.0, panel_z / 1000.0, trimmed)
        last_mesh = ax.pcolormesh(
            mesh_x,
            mesh_z_coords,
            trimmed,
            cmap=cmap,
            norm=norm,
            shading="auto",
            **kwargs
        )
    if last_mesh is None:
        raise ValueError("field_data is empty")
    ax.set_xlabel(labels[0], fontsize=13)
    ax.set_ylabel(labels[1], fontsize=13)
    ax.set_title(title or "", fontsize=15)
    ax.set_ylim(list(height_limit))
    ax.tick_params(axis="both", which="major", labelsize=11, direction="in")
    ax.grid(True, linestyle=":", linewidth=0.7, alpha=0.5)
    add_colorbar(ax, last_mesh, style, fig=fig)
    ax.pycwr_last_mappable = last_mesh
    return last_mesh


def lonlat_section_points(start_xy_km, end_xy_km, distance_from_start_km, radar_lon, radar_lat):
    x_points_km, y_points_km = interpolate_points_along_line(start_xy_km, end_xy_km, distance_from_start_km)
    lon, lat = cartesian_to_geographic_aeqd(
        x_points_km * 1000.0,
        y_points_km * 1000.0,
        lon_0=radar_lon,
        lat_0=radar_lat,
    )
    return lon, lat


def geographic_line_to_cartesian(start_lonlat, end_lonlat, radar_lon, radar_lat):
    start_x, start_y = geographic_to_cartesian_aeqd(
        lat=start_lonlat[1],
        lon=start_lonlat[0],
        lat_0=radar_lat,
        lon_0=radar_lon,
    )
    end_x, end_y = geographic_to_cartesian_aeqd(
        lat=end_lonlat[1],
        lon=end_lonlat[0],
        lat_0=radar_lat,
        lon_0=radar_lon,
    )
    return (start_x[0], start_y[0]), (end_x[0], end_y[0])


def interpolate_points_along_line(start_points, end_points, distance_from_start):
    start_x, start_y = start_points
    end_x, end_y = end_points
    distance_from_start = as_numpy(distance_from_start)
    line_length = np.sqrt((start_x - end_x) ** 2 + (start_y - end_y) ** 2)
    return (
        distance_from_start / line_length * (end_x - start_x) + start_x,
        distance_from_start / line_length * (end_y - start_y) + start_y,
    )
