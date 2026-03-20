# -*- coding: utf-8 -*-
"""Legacy-compatible single-sweep map plotting helpers."""

import matplotlib.pyplot as plt

from ._plot_core import (
    MapOptions,
    PlotStyle,
    ColorbarOptions,
    as_numpy,
    ccrs,
    create_standard_axes,
    default_sweep_title,
    ensure_geographic_coordinates,
    ensure_ppi_scan,
    plot_map,
    resolve_field_data,
    resolve_field_style,
)

DEFAULT_MAP_CRS = ccrs.PlateCarree() if ccrs is not None else None


class RadarGraphMap(object):
    """Backward-compatible mapped PPI plotting helper."""

    def __init__(self, radar=None, **legacy_kwargs):
        legacy_radar = legacy_kwargs.pop("NuistRadar", None)
        if legacy_kwargs:
            raise TypeError("unexpected keyword arguments: %s" % sorted(legacy_kwargs))
        self.NRadar = radar if radar is not None else legacy_radar

    def plot(
        self,
        sweep,
        field_name,
        range_mode=None,
        normvar=None,
        title=None,
        clabel=None,
        continuously=False,
        extent=None,
        data_crs=None,
        map_options=None,
        colorbar_orientation="vertical",
        colorbar_ticks=None,
        colorbar_ticklabels=None,
        **kwargs
    ):
        ensure_ppi_scan(self.NRadar)
        field, field_key = resolve_field_data(self.NRadar, sweep, field_name, range_mode=range_mode)
        style = resolve_field_style(
            self.NRadar,
            sweep,
            field_key,
            field,
            value_range=normvar,
            continuous=continuously,
            colorbar_orientation=colorbar_orientation,
            colorbar_ticks=colorbar_ticks,
            colorbar_ticklabels=colorbar_ticklabels,
            colorbar_label=clabel,
        )
        if title is None:
            title = default_sweep_title(self.NRadar, sweep)
        local_map_options = map_options or MapOptions(data_crs=data_crs or DEFAULT_MAP_CRS)
        fig, ax, cax = create_standard_axes(orientation=colorbar_orientation, projection=ax_projection(local_map_options))
        style.colorbar.cax = cax
        mesh = plot_map(
            ax,
            field.lon,
            field.lat,
            field,
            style,
            fig=fig,
            extent=extent,
            map_options=local_map_options,
            **kwargs
        )
        ax.set_title(title, fontsize=14)
        return mesh

    @staticmethod
    def GUI_plot(
        NRadar,
        fig,
        ax,
        cax,
        sweep,
        field_name,
        range_mode=None,
        normvar=None,
        title=None,
        clabel=None,
        continuously=False,
        extent=None,
        data_crs=None,
        map_options=None,
        colorbar_orientation="vertical",
        colorbar_ticks=None,
        colorbar_ticklabels=None,
        **kwargs
    ):
        ensure_ppi_scan(NRadar)
        field, field_key = resolve_field_data(NRadar, sweep, field_name, range_mode=range_mode)
        style = resolve_field_style(
            NRadar,
            sweep,
            field_key,
            field,
            value_range=normvar,
            continuous=continuously,
            colorbar_orientation=colorbar_orientation,
            colorbar_ticks=colorbar_ticks,
            colorbar_ticklabels=colorbar_ticklabels,
            colorbar_label=clabel,
            colorbar_cax=cax,
        )
        if title is None:
            title = default_sweep_title(NRadar, sweep)
        local_map_options = map_options or MapOptions(data_crs=data_crs or DEFAULT_MAP_CRS)
        mesh = plot_map(
            ax,
            field.lon,
            field.lat,
            field,
            style,
            fig=fig,
            extent=extent,
            map_options=local_map_options,
            **kwargs
        )
        ax.set_title(title, fontsize=14)
        return mesh

    @staticmethod
    def simple_plot_ppi_map(
        radar_data,
        _range=None,
        azimuth=None,
        elevation=None,
        main_piont=None,
        title=None,
        normvar=None,
        cmap=None,
        cmap_bins=16,
        extend=None,
        projection=DEFAULT_MAP_CRS,
        orient="vertical",
        clabel=None,
        continuously=False,
        main_point=None,
        altitude_m=0.0,
        effective_earth_radius=None,
        map_options=None,
        **kwargs
    ):
        station_lonlat = main_point or main_piont
        lon, lat = ensure_geographic_coordinates(
            radar_data=radar_data,
            ranges=_range,
            azimuth=azimuth,
            elevation=elevation,
            station_lonlat=station_lonlat,
            altitude_m=altitude_m,
            effective_earth_radius=effective_earth_radius,
        )
        local_map_options = map_options or MapOptions(data_crs=projection)
        fig, ax, cax = create_standard_axes(orientation=orient, projection=ax_projection(local_map_options))
        return RadarGraphMap.plot_ppi_map(
            fig,
            ax,
            cax,
            lon,
            lat,
            radar_data,
            title=title,
            normvar=normvar,
            cmap=cmap,
            cmap_bins=cmap_bins,
            extend=extend,
            projection=projection,
            orient=orient,
            clabel=clabel,
            continuously=continuously,
            map_options=local_map_options,
            **kwargs
        )

    @staticmethod
    def simple_plot_ppi_xy_map(
        x,
        y,
        radar_data,
        main_piont,
        title=None,
        normvar=None,
        cmap=None,
        cmap_bins=16,
        extend=None,
        projection=DEFAULT_MAP_CRS,
        orient="vertical",
        clabel=None,
        continuously=False,
        map_options=None,
        **kwargs
    ):
        lon, lat = ensure_geographic_coordinates(
            radar_data=radar_data,
            x=x,
            y=y,
            station_lonlat=main_piont,
        )
        local_map_options = map_options or MapOptions(data_crs=projection)
        fig, ax, cax = create_standard_axes(orientation=orient, projection=ax_projection(local_map_options))
        return RadarGraphMap.plot_ppi_map(
            fig,
            ax,
            cax,
            lon,
            lat,
            radar_data,
            title=title,
            normvar=normvar,
            cmap=cmap,
            cmap_bins=cmap_bins,
            extend=extend,
            projection=projection,
            orient=orient,
            clabel=clabel,
            continuously=continuously,
            map_options=local_map_options,
            **kwargs
        )

    @staticmethod
    def plot_ppi_map(
        fig,
        ax,
        cx,
        lon,
        lat,
        radar_data,
        title=None,
        normvar=None,
        cmap=None,
        cmap_bins=16,
        extend=None,
        projection=DEFAULT_MAP_CRS,
        orient="vertical",
        clabel=None,
        continuously=False,
        map_options=None,
        **kwargs
    ):
        style = PlotStyle(
            cmap=cmap or "viridis",
            value_range=normvar,
            bins=cmap_bins,
            continuous=continuously,
            colorbar=ColorbarOptions(
                visible=True,
                orientation=orient,
                label=clabel,
                cax=cx,
            ),
        )
        local_map_options = map_options or MapOptions(data_crs=projection)
        mesh = plot_map(
            ax,
            lon,
            lat,
            radar_data,
            style,
            fig=fig,
            extent=extend,
            map_options=local_map_options,
            **kwargs
        )
        ax.set_title(title or "", fontsize=14)
        return mesh

    @staticmethod
    def _FixTicks(ticks):
        ticks = as_numpy(ticks)
        if (ticks % 1).sum() == 0:
            return ["%2.f" % i for i in ticks]
        return ["%.2f" % i for i in ticks]


def ax_projection(map_options):
    if map_options is None:
        return DEFAULT_MAP_CRS
    return map_options.data_crs or DEFAULT_MAP_CRS
