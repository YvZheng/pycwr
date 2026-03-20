from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt

from .RadarPlot import Graph, GraphMap
from ._plot_core import ccrs, default_sweep_title, plot_vertical_section, resolve_field_data, resolve_field_style


@dataclass
class EasyPlotResult:
    fig: object
    ax: object
    artist: object


def _finalize_plot(fig, ax, artist, show=False, save=None, dpi=150):
    if save is not None:
        save_path = Path(save)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    return EasyPlotResult(fig=fig, ax=ax, artist=artist)


def _new_cartesian_axes(ax=None, figsize=(8, 8)):
    if ax is not None:
        return ax.figure, ax
    fig, ax = plt.subplots(figsize=figsize)
    return fig, ax


def _new_map_axes(ax=None, figsize=(8, 8), projection=None):
    if ax is not None:
        return ax.figure, ax
    projection = projection or (ccrs.PlateCarree() if ccrs is not None else None)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection=projection)
    return fig, ax


def plot_ppi(
    radar,
    field="dBZ",
    sweep=0,
    ax=None,
    figsize=(8, 8),
    show=False,
    save=None,
    title=None,
    **kwargs
):
    """Beginner-friendly Cartesian PPI plot."""
    fig, ax = _new_cartesian_axes(ax=ax, figsize=figsize)
    artist = Graph(radar).plot_ppi(ax, sweep, field, title=title, **kwargs)
    return _finalize_plot(fig, ax, artist, show=show, save=save)


def plot_ppi_map(
    radar,
    field="dBZ",
    sweep=0,
    ax=None,
    figsize=(8, 8),
    projection=None,
    data_crs=None,
    show=False,
    save=None,
    title=None,
    **kwargs
):
    """Beginner-friendly mapped PPI plot."""
    fig, ax = _new_map_axes(ax=ax, figsize=figsize, projection=projection)
    transform = data_crs or (ccrs.PlateCarree() if ccrs is not None else None)
    GraphMap(radar, transform).plot_ppi_map(ax, sweep, field, title=title, **kwargs)
    artist = getattr(ax, "pycwr_last_mappable", None)
    return _finalize_plot(fig, ax, artist, show=show, save=save)


def plot_section(
    radar,
    start=(-50.0, 0.0),
    end=(50.0, 0.0),
    field="dBZ",
    ax=None,
    figsize=(10, 5),
    show=False,
    save=None,
    title=None,
    point_units="km",
    **kwargs
):
    """Beginner-friendly vertical section plot using Cartesian start/end points."""
    fig, ax = _new_cartesian_axes(ax=ax, figsize=figsize)
    artist = Graph(radar).plot_vcs(
        ax,
        start,
        end,
        field,
        title=title,
        point_units=point_units,
        **kwargs
    )
    return _finalize_plot(fig, ax, artist, show=show, save=save)


def plot_section_lonlat(
    radar,
    start_lonlat,
    end_lonlat,
    field="dBZ",
    ax=None,
    figsize=(10, 5),
    show=False,
    save=None,
    title=None,
    **kwargs
):
    """Beginner-friendly vertical section plot using lon/lat endpoints."""
    fig, ax = _new_cartesian_axes(ax=ax, figsize=figsize)
    artist = GraphMap(radar, ccrs.PlateCarree() if ccrs is not None else None).plot_vcs_map(
        ax,
        start_lonlat,
        end_lonlat,
        field,
        title=title,
        **kwargs
    )
    return _finalize_plot(fig, ax, artist, show=show, save=save)


def plot_rhi(
    radar,
    azimuth,
    field="dBZ",
    ax=None,
    figsize=(10, 5),
    show=False,
    save=None,
    title=None,
    **kwargs
):
    """Beginner-friendly RHI-style vertical plot at a target azimuth."""
    fig, ax = _new_cartesian_axes(ax=ax, figsize=figsize)
    min_max = kwargs.pop("min_max", None)
    cmap = kwargs.pop("cmap", None)
    cmap_bins = kwargs.pop("cmap_bins", None)
    cbar = kwargs.pop("cbar", True)
    orientation = kwargs.pop("orientation", "vertical")
    cbar_ticks = kwargs.pop("cbar_ticks", None)
    cbar_ticklabels = kwargs.pop("cbar_ticklabels", None)
    clabel = kwargs.pop("clabel", None)
    continuously = kwargs.pop("continuously", False)
    height_km = kwargs.pop("height_km", None)
    labels = kwargs.pop("labels", None)
    section = radar.extract_rhi(azimuth=azimuth, field_name=field)
    field_data, field_key = resolve_field_data(radar, 0, field)
    style = resolve_field_style(
        radar,
        0,
        field_key,
        field_data,
        value_range=min_max,
        cmap=cmap,
        bins=cmap_bins,
        continuous=continuously,
        colorbar_visible=cbar,
        colorbar_orientation=orientation,
        colorbar_ticks=cbar_ticks,
        colorbar_ticklabels=cbar_ticklabels,
        colorbar_label=clabel,
    )
    artist = plot_vertical_section(
        ax,
        section,
        None,
        None,
        style,
        fig=fig,
        height_km=height_km,
        title=title or default_sweep_title(radar, 0),
        labels=labels,
        **kwargs
    )
    return _finalize_plot(fig, ax, artist, show=show, save=save)


def plot_wind_profile(
    radar_or_profile,
    ax=None,
    figsize=(6, 8),
    show=False,
    save=None,
    title=None,
    **kwargs
):
    """Beginner-friendly vertical wind profile plot."""
    fig, ax = _new_cartesian_axes(ax=ax, figsize=figsize)
    if hasattr(radar_or_profile, "retrieve_vwp"):
        artist = Graph(radar_or_profile).plot_wind_profile(ax, title=title, **kwargs)
    else:
        artist = Graph(None).plot_wind_profile(ax, profile=radar_or_profile, title=title, **kwargs)
    return _finalize_plot(fig, ax, artist, show=show, save=save)


def plot_vvp(
    radar_or_retrieval,
    sweep=0,
    ax=None,
    figsize=(8, 8),
    show=False,
    save=None,
    title=None,
    **kwargs
):
    """Beginner-friendly local VVP wind-field plot."""
    fig, ax = _new_cartesian_axes(ax=ax, figsize=figsize)
    if hasattr(radar_or_retrieval, "retrieve_vvp"):
        artist = Graph(radar_or_retrieval).plot_vvp(ax, sweep=sweep, title=title, **kwargs)
    else:
        artist = Graph(None).plot_vvp(ax, sweep=sweep, retrieval=radar_or_retrieval, title=title, **kwargs)
    return _finalize_plot(fig, ax, artist, show=show, save=save)


def plot(
    radar,
    kind="ppi",
    field="dBZ",
    sweep=0,
    show=False,
    save=None,
    **kwargs
):
    """Single beginner entry point.

    Supported kinds: ``ppi``, ``ppi_map``, ``section``, ``section_lonlat``.
    """
    if kind == "ppi":
        return plot_ppi(radar, field=field, sweep=sweep, show=show, save=save, **kwargs)
    if kind == "ppi_map":
        return plot_ppi_map(radar, field=field, sweep=sweep, show=show, save=save, **kwargs)
    if kind == "section":
        return plot_section(radar, field=field, show=show, save=save, **kwargs)
    if kind == "section_lonlat":
        return plot_section_lonlat(radar, field=field, show=show, save=save, **kwargs)
    if kind == "rhi":
        return plot_rhi(radar, field=field, show=show, save=save, **kwargs)
    if kind == "wind_profile":
        return plot_wind_profile(radar, show=show, save=save, **kwargs)
    if kind == "vvp":
        return plot_vvp(radar, sweep=sweep, show=show, save=save, **kwargs)
    raise ValueError("unsupported kind: %s" % kind)


class EasyRadarPlotter(object):
    """Small wrapper for users who want `plotter.ppi()` style calls."""

    def __init__(self, radar):
        self.radar = radar

    def ppi(self, **kwargs):
        return plot_ppi(self.radar, **kwargs)

    def ppi_map(self, **kwargs):
        return plot_ppi_map(self.radar, **kwargs)

    def section(self, **kwargs):
        return plot_section(self.radar, **kwargs)

    def section_lonlat(self, **kwargs):
        return plot_section_lonlat(self.radar, **kwargs)

    def rhi(self, **kwargs):
        return plot_rhi(self.radar, **kwargs)

    def wind_profile(self, **kwargs):
        return plot_wind_profile(self.radar, **kwargs)

    def vvp(self, **kwargs):
        return plot_vvp(self.radar, **kwargs)

    def quicklook(self, **kwargs):
        return plot(self.radar, **kwargs)
