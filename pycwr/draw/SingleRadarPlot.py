# -*- coding: utf-8 -*-
"""Legacy-compatible single-sweep Cartesian plotting helpers."""

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from ._plot_core import (
    CartesianReferenceOptions,
    ColorbarOptions,
    PlotStyle,
    add_cartesian_reference,
    as_numpy,
    create_standard_axes,
    default_sweep_title,
    ensure_cartesian_coordinates,
    ensure_ppi_scan,
    plot_cartesian,
    resolve_field_data,
    resolve_field_style,
)


class RadarGraph(object):
    """Backward-compatible PPI plotting helper with explicit unit handling."""

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
        dark=False,
        continuously=False,
        extent_km=None,
        axis_labels=None,
        colorbar_orientation="vertical",
        colorbar_ticks=None,
        colorbar_ticklabels=None,
        reference_options=None,
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
        plot_style = "dark_background" if dark else "default"
        with plt.style.context(plot_style):
            fig, ax, cax = create_standard_axes(orientation=colorbar_orientation)
            style.colorbar.cax = cax
            return plot_cartesian(
                ax,
                field.x,
                field.y,
                field,
                style,
                fig=fig,
                extent_km=extent_km,
                title=title,
                labels=axis_labels,
                reference_options=reference_options,
                **kwargs
            )

    @staticmethod
    def GUI_plot(
        NRadar,
        fig,
        ax,
        cx,
        sweep,
        field_name,
        range_mode=None,
        normvar=None,
        title=None,
        clabel=None,
        continuously=False,
        extent_km=None,
        axis_labels=None,
        colorbar_orientation="vertical",
        colorbar_ticks=None,
        colorbar_ticklabels=None,
        reference_options=None,
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
            colorbar_cax=cx,
        )
        if title is None:
            title = default_sweep_title(NRadar, sweep)
        return plot_cartesian(
            ax,
            field.x,
            field.y,
            field,
            style,
            fig=fig,
            extent_km=extent_km,
            title=title,
            labels=axis_labels,
            reference_options=reference_options,
            **kwargs
        )

    @staticmethod
    def simple_plot_ppi_xy(
        x,
        y,
        radar_data,
        normvar=None,
        cmap=None,
        max_range=None,
        title=None,
        cmap_bins=16,
        orient="vertical",
        label=None,
        clabel=None,
        dark=False,
        continuously=False,
        reference_options=None,
        **kwargs
    ):
        plot_style = "dark_background" if dark else "default"
        with plt.style.context(plot_style):
            fig, ax, cax = create_standard_axes(orientation=orient)
            return RadarGraph.plot_ppi(
                fig,
                ax,
                cax,
                x,
                y,
                radar_data,
                max_range=max_range,
                title=title,
                normvar=normvar,
                cmap=cmap,
                cmap_bins=cmap_bins,
                orient=orient,
                label=label,
                clabel=clabel,
                continuously=continuously,
                reference_options=reference_options,
                **kwargs
            )

    @staticmethod
    def simple_plot_ppi(
        radar_data,
        _range=None,
        azimuth=None,
        elevation=None,
        normvar=None,
        cmap=None,
        max_range=None,
        title=None,
        cmap_bins=16,
        orient="vertical",
        label=None,
        clabel=None,
        dark=False,
        continuously=False,
        altitude_m=0.0,
        effective_earth_radius=None,
        reference_options=None,
        **kwargs
    ):
        x, y = ensure_cartesian_coordinates(
            radar_data=radar_data,
            ranges=_range,
            azimuth=azimuth,
            elevation=elevation,
            altitude_m=altitude_m,
            effective_earth_radius=effective_earth_radius,
        )
        return RadarGraph.simple_plot_ppi_xy(
            x,
            y,
            radar_data,
            normvar=normvar,
            cmap=cmap,
            max_range=max_range,
            title=title,
            cmap_bins=cmap_bins,
            orient=orient,
            label=label,
            clabel=clabel,
            dark=dark,
            continuously=continuously,
            reference_options=reference_options,
            **kwargs
        )

    @staticmethod
    def plot_ppi(
        fig,
        ax,
        cx,
        x,
        y,
        radar_data,
        max_range=None,
        title=None,
        normvar=None,
        cmap=None,
        cmap_bins=16,
        orient="vertical",
        label=None,
        clabel=None,
        continuously=False,
        reference_options=None,
        **kwargs
    ):
        extent_km = max_range
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
        return plot_cartesian(
            ax,
            x,
            y,
            radar_data,
            style,
            fig=fig,
            extent_km=extent_km,
            title=title,
            labels=label,
            reference_options=reference_options,
            **kwargs
        )

    @staticmethod
    def _SetAxis(ax, major, minor, fontsize):
        ax.xaxis.set_major_locator(MultipleLocator(major))
        ax.xaxis.set_minor_locator(MultipleLocator(minor))
        ax.yaxis.set_major_locator(MultipleLocator(major))
        ax.yaxis.set_minor_locator(MultipleLocator(minor))
        ax.tick_params(axis="both", which="major", labelsize=fontsize)

    @staticmethod
    def _SetGrids(ax, max_range, deltaR=50):
        extent = (-max_range, max_range, -max_range, max_range)
        add_cartesian_reference(
            ax,
            extent,
            reference_options=CartesianReferenceOptions(ring_spacing_km=deltaR),
        )

    @staticmethod
    def _FixTicks(ticks):
        ticks = as_numpy(ticks)
        if (ticks % 1).sum() == 0:
            return ["%2.f" % i for i in ticks]
        return ["%.2f" % i for i in ticks]
