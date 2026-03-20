# -*- coding: utf-8 -*-
"""Vertical section plotting helpers with explicit coordinate handling."""

import matplotlib.pyplot as plt

from ._plot_core import (
    ColorbarOptions,
    PlotStyle,
    as_numpy,
    create_standard_axes,
    default_sweep_title,
    ensure_ppi_scan,
    geographic_line_to_cartesian,
    interpolate_points_along_line,
    lonlat_section_points,
    plot_vertical_section,
    resolve_field_data,
    resolve_field_style,
)


class VerticalSection(object):
    """Backward-compatible RHI and VCS plotting helper."""

    def __init__(self, radar=None, **legacy_kwargs):
        legacy_radar = legacy_kwargs.pop("NuistRadar", None)
        if legacy_kwargs:
            raise TypeError("unexpected keyword arguments: %s" % sorted(legacy_kwargs))
        self.NRadar = radar if radar is not None else legacy_radar

    def RHI(
        self,
        azimuth,
        field_name,
        range_mode=None,
        height=(0, 18),
        title=None,
        clabel=None,
        continuously=False,
        normvar=None,
        cmap=None,
        cmap_bins=None,
        colorbar_orientation="horizontal",
        colorbar_ticks=None,
        colorbar_ticklabels=None,
        axis_labels=None,
        **kwargs
    ):
        section = self.NRadar.extract_rhi(azimuth=azimuth, field_name=field_name, range_mode=range_mode)
        field, field_key = resolve_field_data(self.NRadar, 0, field_name, range_mode=range_mode)
        style = resolve_field_style(
            self.NRadar,
            0,
            field_key,
            field,
            value_range=normvar,
            cmap=cmap,
            bins=cmap_bins,
            continuous=continuously,
            colorbar_orientation=colorbar_orientation,
            colorbar_ticks=colorbar_ticks,
            colorbar_ticklabels=colorbar_ticklabels,
            colorbar_label=clabel,
        )
        fig, ax, cax = create_standard_axes(orientation=colorbar_orientation)
        style.colorbar.cax = cax
        if title is None:
            title = default_sweep_title(self.NRadar, 0)
        return plot_vertical_section(
            ax,
            section,
            None,
            None,
            style,
            fig=fig,
            height_km=height,
            title=title,
            labels=axis_labels,
            **kwargs
        )

    def section(
        self,
        start_point,
        end_point,
        field_name,
        range_mode=None,
        height=(0, 18),
        title=None,
        clabel=None,
        continuously=False,
        normvar=None,
        cmap=None,
        cmap_bins=None,
        colorbar_orientation="horizontal",
        colorbar_ticks=None,
        colorbar_ticklabels=None,
        axis_labels=None,
        point_units="km",
        **kwargs
    ):
        ensure_ppi_scan(self.NRadar)
        scale = 1000.0 if point_units == "km" else 1.0
        section = self.NRadar.extract_section(
            (start_point[0] * scale, start_point[1] * scale),
            (end_point[0] * scale, end_point[1] * scale),
            field_name=field_name,
            point_units="m",
            range_mode=range_mode,
        )
        field, field_key = resolve_field_data(self.NRadar, 0, field_name, range_mode=range_mode)
        style = resolve_field_style(
            self.NRadar,
            0,
            field_key,
            field,
            value_range=normvar,
            cmap=cmap,
            bins=cmap_bins,
            continuous=continuously,
            colorbar_orientation=colorbar_orientation,
            colorbar_ticks=colorbar_ticks,
            colorbar_ticklabels=colorbar_ticklabels,
            colorbar_label=clabel,
        )
        fig, ax, cax = create_standard_axes(orientation=colorbar_orientation)
        style.colorbar.cax = cax
        if title is None:
            title = default_sweep_title(self.NRadar, 0)
        return plot_vertical_section(
            ax,
            section,
            None,
            None,
            style,
            fig=fig,
            height_km=height,
            title=title,
            labels=axis_labels,
            **kwargs
        )

    def section_map(
        self,
        start_lonlat,
        end_lonlat,
        field_name,
        range_mode=None,
        height=(0, 18),
        title=None,
        orient="horizontal",
        label=None,
        clabel=None,
        continuously=False,
        normvar=None,
        cmap=None,
        cmap_bins=None,
        colorbar_ticks=None,
        colorbar_ticklabels=None,
        **kwargs
    ):
        ensure_ppi_scan(self.NRadar)
        fig, ax, cax = create_standard_axes(orientation=orient)
        return VerticalSection.GUI_section_map(
            fig,
            ax,
            cax,
            self.NRadar,
            start_lonlat,
            end_lonlat,
            field_name,
            range_mode=range_mode,
            title=title,
            clabel=clabel,
            continuously=continuously,
            normvar=normvar,
            cmap=cmap,
            cmap_bins=cmap_bins,
            colorbar_ticks=colorbar_ticks,
            colorbar_ticklabels=colorbar_ticklabels,
            height=height,
            label=label,
            **kwargs
        )

    @staticmethod
    def GUI_section(
        fig,
        ax,
        cax,
        NRadar,
        start_point,
        end_point,
        field_name,
        range_mode=None,
        title=None,
        clabel=None,
        continuously=False,
        colorbar_orientation="horizontal",
        colorbar_ticks=None,
        colorbar_ticklabels=None,
        normvar=None,
        cmap=None,
        cmap_bins=None,
        height=(0, 18),
        label=None,
        **kwargs
    ):
        ensure_ppi_scan(NRadar)
        section = NRadar.extract_section(
            start_point,
            end_point,
            field_name=field_name,
            point_units="m",
            range_mode=range_mode,
        )
        field, field_key = resolve_field_data(NRadar, 0, field_name, range_mode=range_mode)
        style = resolve_field_style(
            NRadar,
            0,
            field_key,
            field,
            value_range=normvar,
            cmap=cmap,
            bins=cmap_bins,
            continuous=continuously,
            colorbar_orientation=colorbar_orientation,
            colorbar_ticks=colorbar_ticks,
            colorbar_ticklabels=colorbar_ticklabels,
            colorbar_label=clabel,
            colorbar_cax=cax,
        )
        if title is None:
            title = default_sweep_title(NRadar, 0)
        return plot_vertical_section(
            ax,
            section,
            None,
            None,
            style,
            fig=fig,
            height_km=height,
            title=title,
            labels=label,
            **kwargs
        )

    @staticmethod
    def GUI_section_map(
        fig,
        ax,
        cax,
        NRadar,
        start_lonlat,
        end_lonlat,
        field_name,
        range_mode=None,
        title=None,
        clabel=None,
        continuously=False,
        colorbar_ticks=None,
        colorbar_ticklabels=None,
        normvar=None,
        cmap=None,
        cmap_bins=None,
        height=(0, 18),
        label=None,
        orient="horizontal",
        **kwargs
    ):
        ensure_ppi_scan(NRadar)
        field, field_key = resolve_field_data(NRadar, 0, field_name, range_mode=range_mode)
        style = resolve_field_style(
            NRadar,
            0,
            field_key,
            field,
            value_range=normvar,
            cmap=cmap,
            bins=cmap_bins,
            continuous=continuously,
            colorbar_orientation=orient,
            colorbar_ticks=colorbar_ticks,
            colorbar_ticklabels=colorbar_ticklabels,
            colorbar_label=clabel,
            colorbar_cax=cax,
        )
        if title is None:
            title = default_sweep_title(NRadar, 0)
        radar_lon = float(as_numpy(NRadar.scan_info.longitude))
        radar_lat = float(as_numpy(NRadar.scan_info.latitude))
        start_xy, end_xy = geographic_line_to_cartesian(start_lonlat, end_lonlat, radar_lon, radar_lat)
        section = NRadar.extract_section_lonlat(
            start_lonlat,
            end_lonlat,
            field_name=field_name,
            range_mode=range_mode,
        )
        mesh = plot_vertical_section(
            ax,
            section,
            None,
            None,
            style,
            fig=fig,
            height_km=height,
            title=title,
            labels=label,
            **kwargs
        )
        xticks = ax.get_xticks()
        start_xy_km = (start_xy[0] / 1000.0, start_xy[1] / 1000.0)
        end_xy_km = (end_xy[0] / 1000.0, end_xy[1] / 1000.0)
        lon_points, lat_points = lonlat_section_points(start_xy_km, end_xy_km, xticks, radar_lon, radar_lat)
        ax.set_xticklabels(
            ["(%.2f, %.2f)" % (lon_points[i], lat_points[i]) for i in range(len(xticks))],
            rotation=15,
            fontsize=10,
        )
        return mesh

    @staticmethod
    def SectionPlot_VCS(
        fig,
        ax,
        cx,
        mesh_xy,
        mesh_z,
        field_data,
        height=(0, 18),
        title=None,
        normvar=None,
        cmap=None,
        cmap_bins=16,
        orient="horizontal",
        label=None,
        clabel=None,
        continuously=False,
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
        return plot_vertical_section(
            ax,
            mesh_xy,
            mesh_z,
            field_data,
            style,
            fig=fig,
            height_km=height,
            title=title,
            labels=label,
            **kwargs
        )

    @staticmethod
    def SectionPlot_VCS_map(
        fig,
        ax,
        cx,
        start_lonlat,
        end_lonlat,
        field_name,
        NRadar,
        range_mode=None,
        height=(0, 18),
        title=None,
        normvar=None,
        cmap=None,
        cmap_bins=16,
        orient="horizontal",
        label=None,
        clabel=None,
        continuously=False,
        **kwargs
    ):
        mesh = VerticalSection.GUI_section_map(
            fig,
            ax,
            cx,
            NRadar,
            start_lonlat,
            end_lonlat,
            field_name,
            range_mode=range_mode,
            title=title,
            clabel=clabel,
            continuously=continuously,
            height=height,
            label=label,
            orient=orient,
            normvar=normvar,
            cmap=cmap,
            cmap_bins=cmap_bins,
            **kwargs
        )
        return mesh

    @staticmethod
    def _FixTicks(ticks):
        ticks = as_numpy(ticks)
        if (ticks % 1).sum() == 0:
            return ["%2.f" % i for i in ticks]
        return ["%.2f" % i for i in ticks]

    @staticmethod
    def get_points_from_ranges(start_points, end_points, distance_from_start):
        return interpolate_points_along_line(start_points, end_points, distance_from_start)
