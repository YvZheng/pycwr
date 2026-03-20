import matplotlib.pyplot as plt
import numpy as np

from ..configure.default_config import (
    CINRAD_COLORMAP,
    CINRAD_field_bins,
    CINRAD_field_mapping,
    CINRAD_field_normvar,
)
from ..core.transforms import antenna_vectors_to_cartesian, cartesian_to_geographic_aeqd
from ._plot_core import (
    CartesianReferenceOptions,
    ColorbarOptions,
    MapOptions,
    PlotStyle,
    add_cartesian_reference,
    as_numpy,
    default_colorbar_label,
    default_sweep_title,
    geographic_line_to_cartesian,
    lonlat_section_points,
    plot_cartesian,
    plot_map,
    plot_vertical_section,
    require_cartopy_axes,
    require_mpl_axes,
    resolve_field_data,
    resolve_field_style,
)


def _require_mpl_axes(ax):
    require_mpl_axes(ax)


def _require_cartopy_axes(ax):
    require_cartopy_axes(ax)


class Graph(object):
    """Recommended Cartesian plotting interface with unit-explicit options."""

    def __init__(self, NRadar):
        self.Radar = NRadar

    def _field_style(
        self,
        sweep_num,
        field_name,
        field_data,
        cmap=None,
        min_max=None,
        cmap_bins=None,
        cbar=True,
        orientation="vertical",
        cbar_ticks=None,
        cbar_ticklabels=None,
        clabel=None,
        continuously=False,
    ):
        return resolve_field_style(
            self.Radar,
            sweep_num,
            field_name,
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

    def plot_ppi(
        self,
        ax,
        sweep_num,
        field_name,
        range_mode=None,
        cmap=None,
        min_max=None,
        cmap_bins=None,
        cbar=True,
        orientation="vertical",
        cbar_ticks=None,
        cbar_ticklabels=None,
        clabel=None,
        title=None,
        extent_km=None,
        labels=None,
        reference_options=None,
        continuously=False,
        **kwargs
    ):
        _require_mpl_axes(ax)
        field, field_key = resolve_field_data(self.Radar, sweep_num, field_name, range_mode=range_mode)
        style = self._field_style(
            sweep_num,
            field_key,
            field,
            cmap=cmap,
            min_max=min_max,
            cmap_bins=cmap_bins,
            cbar=cbar,
            orientation=orientation,
            cbar_ticks=cbar_ticks,
            cbar_ticklabels=cbar_ticklabels,
            clabel=clabel,
            continuously=continuously,
        )
        if title is None:
            title = default_sweep_title(self.Radar, sweep_num)
        return plot_cartesian(
            ax,
            field.x,
            field.y,
            field,
            style,
            fig=ax.figure,
            extent_km=extent_km,
            title=title,
            labels=labels,
            reference_options=reference_options,
            **kwargs
        )

    def plot_rhi(
        self,
        ax,
        sweep_num,
        field_name,
        range_mode=None,
        cmap=None,
        min_max=None,
        cmap_bins=None,
        cbar=True,
        orientation="vertical",
        cbar_ticks=None,
        cbar_ticklabels=None,
        clabel=None,
        title=None,
        height_km=None,
        labels=None,
        continuously=False,
        **kwargs
    ):
        _require_mpl_axes(ax)
        field, field_key = resolve_field_data(self.Radar, sweep_num, field_name, range_mode=range_mode)
        style = self._field_style(
            sweep_num,
            field_key,
            field,
            cmap=cmap,
            min_max=min_max,
            cmap_bins=cmap_bins,
            cbar=cbar,
            orientation=orientation,
            cbar_ticks=cbar_ticks,
            cbar_ticklabels=cbar_ticklabels,
            clabel=clabel,
            continuously=continuously,
        )
        if title is None:
            title = default_sweep_title(self.Radar, sweep_num)
        mesh_xy = np.hypot(as_numpy(field.x), as_numpy(field.y))
        mesh_z = as_numpy(field.z)
        if height_km is None:
            height_km = (0.0, float(np.nanmax(mesh_z) / 1000.0))
        return plot_vertical_section(
            ax,
            [mesh_xy],
            [mesh_z],
            [field],
            style,
            fig=ax.figure,
            height_km=height_km,
            title=title,
            labels=labels,
            **kwargs
        )

    def plot_vcs(
        self,
        ax,
        start_xy,
        end_xy,
        field_name,
        range_mode=None,
        cmap=None,
        min_max=None,
        cmap_bins=None,
        cbar=True,
        orientation="vertical",
        cbar_ticks=None,
        cbar_ticklabels=None,
        clabel=None,
        title=None,
        height_km=(0, 18),
        labels=None,
        continuously=False,
        point_units="km",
        **kwargs
    ):
        _require_mpl_axes(ax)
        field, field_key = resolve_field_data(self.Radar, 0, field_name, range_mode=range_mode)
        style = self._field_style(
            0,
            field_key,
            field,
            cmap=cmap,
            min_max=min_max,
            cmap_bins=cmap_bins,
            cbar=cbar,
            orientation=orientation,
            cbar_ticks=cbar_ticks,
            cbar_ticklabels=cbar_ticklabels,
            clabel=clabel,
            continuously=continuously,
        )
        scale = 1000.0 if point_units == "km" else 1.0
        section = self.Radar.extract_section(
            (start_xy[0] * scale, start_xy[1] * scale),
            (end_xy[0] * scale, end_xy[1] * scale),
            field_name=field_name,
            point_units="m",
            range_mode=range_mode,
        )
        if title is None:
            title = default_sweep_title(self.Radar, 0)
        return plot_vertical_section(
            ax,
            section,
            None,
            None,
            style,
            fig=ax.figure,
            height_km=height_km,
            title=title,
            labels=labels,
            **kwargs
        )

    def plot_crf(
        self,
        ax,
        range_mode=None,
        cmap=CINRAD_COLORMAP[CINRAD_field_mapping["dBZ"]],
        min_max=CINRAD_field_normvar[CINRAD_field_mapping["dBZ"]],
        cmap_bins=CINRAD_field_bins[CINRAD_field_mapping["dBZ"]],
        cbar=True,
        orientation="vertical",
        cbar_ticks=None,
        cbar_ticklabels=None,
        clabel=None,
        title="Composite Reflectivity",
        extent_km=None,
        labels=None,
        reference_options=None,
        continuously=False,
        **kwargs
    ):
        _require_mpl_axes(ax)
        reference_field = self.Radar.get_sweep_field(0, "dBZ", range_mode=range_mode)
        max_range = int(reference_field.range.max().values)
        x_range = np.arange(-1 * max_range, max_range + 1, 1000.0)
        y_range = x_range
        self.Radar.add_product_CR_xy(x_range, y_range, range_mode=range_mode)
        product_name = "CR" if range_mode == "aligned" else "CR_native"
        radar_data = self.Radar.product[product_name].values
        x, y = np.meshgrid(self.Radar.product["CR"].x_cr.values, self.Radar.product["CR"].y_cr.values, indexing="ij")
        style = PlotStyle(
            cmap=cmap,
            value_range=min_max,
            bins=cmap_bins,
            continuous=continuously,
            colorbar=ColorbarOptions(
                visible=cbar,
                orientation=orientation,
                label=clabel or "CR (dBZ)",
                ticks=cbar_ticks,
                ticklabels=cbar_ticklabels,
            ),
        )
        return plot_cartesian(
            ax,
            x,
            y,
            radar_data,
            style,
            fig=ax.figure,
            extent_km=extent_km,
            title=title,
            labels=labels,
            reference_options=reference_options,
            **kwargs
        )

    def plot_cappi(
        self,
        ax,
        level_height=3000,
        range_mode=None,
        cmap=CINRAD_COLORMAP[CINRAD_field_mapping["dBZ"]],
        min_max=CINRAD_field_normvar[CINRAD_field_mapping["dBZ"]],
        cmap_bins=CINRAD_field_bins[CINRAD_field_mapping["dBZ"]],
        cbar=True,
        orientation="vertical",
        cbar_ticks=None,
        cbar_ticklabels=None,
        clabel=None,
        title=None,
        extent_km=None,
        labels=None,
        reference_options=None,
        continuously=False,
        **kwargs
    ):
        _require_mpl_axes(ax)
        reference_field = self.Radar.get_sweep_field(0, "dBZ", range_mode=range_mode)
        max_range = int(reference_field.range.max().values)
        x_range = np.arange(-1 * max_range, max_range + 1, 1000.0)
        y_range = x_range
        self.Radar.add_product_CAPPI_xy(x_range, y_range, level_height, range_mode=range_mode)
        product_name = "CAPPI_%d" % level_height if range_mode == "aligned" else "CAPPI_%d_native" % level_height
        radar_data = self.Radar.product[product_name].values
        x, y = np.meshgrid(
            self.Radar.product["CAPPI_%d" % level_height]["x_cappi_%d" % level_height].values,
            self.Radar.product["CAPPI_%d" % level_height]["y_cappi_%d" % level_height].values,
            indexing="ij",
        )
        if title is None:
            title = "CAPPI %dm" % level_height
        style = PlotStyle(
            cmap=cmap,
            value_range=min_max,
            bins=cmap_bins,
            continuous=continuously,
            colorbar=ColorbarOptions(
                visible=cbar,
                orientation=orientation,
                label=clabel or "CAPPI (dBZ)",
                ticks=cbar_ticks,
                ticklabels=cbar_ticklabels,
            ),
        )
        return plot_cartesian(
            ax,
            x,
            y,
            radar_data,
            style,
            fig=ax.figure,
            extent_km=extent_km,
            title=title,
            labels=labels,
            reference_options=reference_options,
            **kwargs
        )

    def plot_wind_profile(
        self,
        ax,
        profile=None,
        sweeps=None,
        field_name=None,
        range_mode="aligned",
        title=None,
        color="C0",
        marker="o",
        linewidth=1.5,
        show_barbs=True,
        barb_color="0.2",
        annotate_samples=True,
        **kwargs
    ):
        _require_mpl_axes(ax)
        if profile is None:
            profile = self.Radar.retrieve_vwp(sweeps=sweeps, field_name=field_name, range_mode=range_mode, **kwargs)
        height_km = np.asarray(profile["height"].values, dtype=np.float64) / 1000.0
        speed = np.asarray(profile["wind_speed"].values, dtype=np.float64)
        u = np.asarray(profile["u"].values, dtype=np.float64)
        v = np.asarray(profile["v"].values, dtype=np.float64)
        sample_count = np.asarray(profile["sample_count"].values, dtype=np.float64) if "sample_count" in profile else None
        line = ax.plot(speed, height_km, color=color, marker=marker, linewidth=linewidth)
        ax.set_xlabel("Wind Speed (m/s)")
        ax.set_ylabel("Height (km)")
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
        ax.set_facecolor("#f7f8fa")
        if title is None:
            title = "Vertical Wind Profile"
        ax.set_title(title)
        if sample_count is not None:
            finite_speed = np.isfinite(speed) & np.isfinite(height_km) & np.isfinite(sample_count)
            if np.any(finite_speed):
                ax.scatter(
                    speed[finite_speed],
                    height_km[finite_speed],
                    c=sample_count[finite_speed],
                    s=40,
                    cmap="plasma",
                    edgecolors="white",
                    linewidths=0.6,
                    zorder=3,
                )
        if show_barbs:
            barb_ax = ax.twiny()
            barb_ax.set_xlim(-1.0, 1.0)
            barb_ax.set_xticks([])
            barb_ax.patch.set_alpha(0.0)
            valid = np.isfinite(u) & np.isfinite(v) & np.isfinite(height_km)
            barbs = barb_ax.barbs(
                np.zeros(int(valid.sum()), dtype=np.float64),
                height_km[valid],
                u[valid],
                v[valid],
                length=5.5,
                barbcolor=barb_color,
                flagcolor=barb_color,
            )
            if annotate_samples and sample_count is not None:
                for height_value, count_value in zip(height_km[valid], sample_count[valid]):
                    barb_ax.text(0.12, height_value, str(int(count_value)), fontsize=8, color=barb_color, va="center")
            ax._pycwr_barb_axes = barb_ax
            return barbs
        return line[0]

    def plot_vvp(
        self,
        ax,
        sweep,
        retrieval=None,
        field_name=None,
        range_mode="aligned",
        title=None,
        background="speed",
        background_field=None,
        background_alpha=0.8,
        background_cmap=None,
        background_min_max=None,
        background_cmap_bins=None,
        cbar=True,
        quiver_color="black",
        quiver_scale=None,
        extent_km=None,
        **kwargs
    ):
        _require_mpl_axes(ax)
        if retrieval is None:
            retrieval = self.Radar.retrieve_vvp(sweep, field_name=field_name, range_mode=range_mode, **kwargs)
        x_km = np.asarray(retrieval["x"].values, dtype=np.float64) / 1000.0
        y_km = np.asarray(retrieval["y"].values, dtype=np.float64) / 1000.0
        u = np.asarray(retrieval["u"].values, dtype=np.float64)
        v = np.asarray(retrieval["v"].values, dtype=np.float64)
        speed = np.asarray(retrieval["wind_speed"].values, dtype=np.float64)
        artist = None
        if background_field is not None:
            background_data, background_key = resolve_field_data(
                self.Radar,
                sweep,
                background_field,
                range_mode=range_mode,
            )
            style = self._field_style(
                sweep,
                background_key,
                background_data,
                cmap=background_cmap,
                min_max=background_min_max,
                cmap_bins=background_cmap_bins,
                cbar=False,
            )
            x_bg = np.asarray(background_data["x"].values, dtype=np.float64) / 1000.0
            y_bg = np.asarray(background_data["y"].values, dtype=np.float64) / 1000.0
            value_bg = np.asarray(background_data.values, dtype=np.float64)
            max_range_km = kwargs.get("max_range_km")
            valid_bg = np.isfinite(x_bg) & np.isfinite(y_bg) & np.isfinite(value_bg)
            if max_range_km is not None:
                valid_bg &= (np.hypot(x_bg, y_bg) <= float(max_range_km))
            vmin = None if style.value_range is None else style.value_range[0]
            vmax = None if style.value_range is None else style.value_range[1]
            artist = ax.scatter(
                x_bg[valid_bg],
                y_bg[valid_bg],
                c=value_bg[valid_bg],
                s=8,
                cmap=style.cmap,
                vmin=vmin,
                vmax=vmax,
                linewidths=0.0,
                alpha=background_alpha,
            )
            if cbar:
                label = default_colorbar_label(background_key)
                ax.figure.colorbar(artist, ax=ax, label=label)
        elif background == "speed":
            finite_speed = np.isfinite(speed) & np.isfinite(x_km) & np.isfinite(y_km)
            artist = ax.scatter(
                x_km[finite_speed],
                y_km[finite_speed],
                c=speed[finite_speed],
                s=18,
                cmap="viridis",
                linewidths=0.0,
            )
            if cbar:
                ax.figure.colorbar(artist, ax=ax, label="Wind Speed (m/s)")
        valid = np.isfinite(u) & np.isfinite(v) & np.isfinite(x_km) & np.isfinite(y_km)
        quiver = ax.quiver(
            x_km[valid],
            y_km[valid],
            u[valid],
            v[valid],
            angles="xy",
            scale_units="xy",
            scale=quiver_scale,
            color=quiver_color,
        )
        if title is None:
            title = "Local VVP Wind Field"
        ax.set_title(title)
        ax.set_xlabel("East Distance (km)")
        ax.set_ylabel("North Distance (km)")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
        ax.set_facecolor("#f7f8fa")
        if extent_km is not None:
            ax.set_xlim(extent_km[0], extent_km[1])
            ax.set_ylim(extent_km[2], extent_km[3])
        ax._pycwr_last_mappable = artist
        return quiver

    def add_rings(self, ax, rings, color="#5B5B5B", linestyle="-", linewidth=0.6, **kwargs):
        theta = np.linspace(0, 2 * np.pi, 200)
        for radius in rings:
            x0 = radius * np.cos(theta)
            y0 = radius * np.sin(theta)
            ax.plot(x0, y0, linestyle=linestyle, linewidth=linewidth, color=color, **kwargs)
        for rad in np.arange(0, np.pi, np.pi / 6):
            ax.plot(
                [-1 * rings[-1] * np.sin(rad), rings[-1] * np.sin(rad)],
                [-1 * rings[-1] * np.cos(rad), rings[-1] * np.cos(rad)],
                linestyle=linestyle,
                linewidth=linewidth,
                color=color,
                **kwargs
            )
        return ax

    def add_lines(self, ax, start_xy, end_xy, color="red", marker="x", markersize=12, **kwargs):
        line_x = [start_xy[0], end_xy[0]]
        line_y = [start_xy[1], end_xy[1]]
        return ax.plot(line_x, line_y, color=color, marker=marker, markersize=markersize, **kwargs)


class GraphMap(object):
    """Recommended mapped plotting interface with explicit CRS handling."""

    def __init__(self, NRadar, transform):
        self.Radar = NRadar
        self.transform = transform

    def plot_ppi_map(
        self,
        ax,
        sweep_num,
        field_name,
        range_mode=None,
        extend=None,
        cmap=None,
        min_max=None,
        cmap_bins=None,
        cbar=True,
        orientation="vertical",
        cbar_ticks=None,
        cbar_ticklabels=None,
        clabel=None,
        title=None,
        continuously=False,
        map_options=None,
        **kwargs
    ):
        _require_cartopy_axes(ax)
        field, field_key = resolve_field_data(self.Radar, sweep_num, field_name, range_mode=range_mode)
        style = resolve_field_style(
            self.Radar,
            sweep_num,
            field_key,
            field,
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
        mesh = plot_map(
            ax,
            field.lon,
            field.lat,
            field,
            style,
            fig=ax.figure,
            extent=extend,
            map_options=map_options or MapOptions(data_crs=self.transform),
            **kwargs
        )
        ax.set_title(title or default_sweep_title(self.Radar, sweep_num), fontsize=14)
        return ax

    def plot_cappi_map(
        self,
        ax,
        level_height,
        range_mode=None,
        extend=None,
        cmap=CINRAD_COLORMAP[CINRAD_field_mapping["dBZ"]],
        min_max=CINRAD_field_normvar[CINRAD_field_mapping["dBZ"]],
        cmap_bins=CINRAD_field_bins[CINRAD_field_mapping["dBZ"]],
        cbar=True,
        orientation="vertical",
        cbar_ticks=None,
        cbar_ticklabels=None,
        clabel=None,
        title=None,
        continuously=False,
        map_options=None,
        **kwargs
    ):
        _require_cartopy_axes(ax)
        if extend is None:
            reference_field = self.Radar.get_sweep_field(0, "dBZ", range_mode=range_mode)
            extend = (
                float(np.min(reference_field.lon)),
                float(np.max(reference_field.lon)),
                float(np.min(reference_field.lat)),
                float(np.max(reference_field.lat)),
            )
        x_lon = np.arange(extend[0], extend[1], 0.01)
        y_lat = np.arange(extend[2], extend[3], 0.01)
        self.Radar.add_product_CAPPI_lonlat(x_lon, y_lat, level_height, range_mode=range_mode)
        product_name = "CAPPI_geo_%d" % level_height if range_mode == "aligned" else "CAPPI_geo_%d_native" % level_height
        radar_data = self.Radar.product[product_name].values
        lon, lat = np.meshgrid(x_lon, y_lat, indexing="ij")
        style = PlotStyle(
            cmap=cmap,
            value_range=min_max,
            bins=cmap_bins,
            continuous=continuously,
            colorbar=ColorbarOptions(
                visible=cbar,
                orientation=orientation,
                label=clabel or "CAPPI (dBZ)",
                ticks=cbar_ticks,
                ticklabels=cbar_ticklabels,
            ),
        )
        plot_map(
            ax,
            lon,
            lat,
            radar_data,
            style,
            fig=ax.figure,
            extent=extend,
            map_options=map_options or MapOptions(data_crs=self.transform),
            **kwargs
        )
        ax.set_title(title or "CAPPI %dm" % level_height, fontsize=14)
        return ax

    def plot_crf_map(
        self,
        ax,
        range_mode=None,
        extend=None,
        cmap=CINRAD_COLORMAP[CINRAD_field_mapping["dBZ"]],
        min_max=CINRAD_field_normvar[CINRAD_field_mapping["dBZ"]],
        cmap_bins=CINRAD_field_bins[CINRAD_field_mapping["dBZ"]],
        cbar=True,
        orientation="vertical",
        cbar_ticks=None,
        cbar_ticklabels=None,
        clabel=None,
        title="Composite Reflectivity",
        continuously=False,
        map_options=None,
        **kwargs
    ):
        _require_cartopy_axes(ax)
        if extend is None:
            reference_field = self.Radar.get_sweep_field(0, "dBZ", range_mode=range_mode)
            extend = (
                float(np.min(reference_field.lon)),
                float(np.max(reference_field.lon)),
                float(np.min(reference_field.lat)),
                float(np.max(reference_field.lat)),
            )
        x_lon = np.arange(extend[0], extend[1], 0.01)
        y_lat = np.arange(extend[2], extend[3], 0.01)
        self.Radar.add_product_CR_lonlat(x_lon, y_lat, range_mode=range_mode)
        product_name = "CR_geo" if range_mode == "aligned" else "CR_geo_native"
        radar_data = self.Radar.product[product_name].values
        lon, lat = np.meshgrid(x_lon, y_lat, indexing="ij")
        style = PlotStyle(
            cmap=cmap,
            value_range=min_max,
            bins=cmap_bins,
            continuous=continuously,
            colorbar=ColorbarOptions(
                visible=cbar,
                orientation=orientation,
                label=clabel or "CR (dBZ)",
                ticks=cbar_ticks,
                ticklabels=cbar_ticklabels,
            ),
        )
        plot_map(
            ax,
            lon,
            lat,
            radar_data,
            style,
            fig=ax.figure,
            extent=extend,
            map_options=map_options or MapOptions(data_crs=self.transform),
            **kwargs
        )
        ax.set_title(title, fontsize=14)
        return ax

    def plot_vcs_map(
        self,
        ax,
        start_lonlat,
        end_lonlat,
        field_name,
        range_mode=None,
        cmap=None,
        min_max=None,
        cmap_bins=None,
        cbar=True,
        orientation="vertical",
        cbar_ticks=None,
        cbar_ticklabels=None,
        clabel=None,
        title=None,
        height_km=(0, 18),
        labels=None,
        continuously=False,
        **kwargs
    ):
        _require_mpl_axes(ax)
        field, field_key = resolve_field_data(self.Radar, 0, field_name, range_mode=range_mode)
        style = resolve_field_style(
            self.Radar,
            0,
            field_key,
            field,
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
        radar_lon = float(as_numpy(self.Radar.scan_info.longitude))
        radar_lat = float(as_numpy(self.Radar.scan_info.latitude))
        start_xy, end_xy = geographic_line_to_cartesian(start_lonlat, end_lonlat, radar_lon, radar_lat)
        section = self.Radar.extract_section_lonlat(
            start_lonlat,
            end_lonlat,
            field_name=field_name,
            range_mode=range_mode,
        )
        gci = plot_vertical_section(
            ax,
            section,
            None,
            None,
            style,
            fig=ax.figure,
            height_km=height_km,
            title=title or default_sweep_title(self.Radar, 0),
            labels=labels,
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
        return gci

    def add_lines_map(self, ax, start_lonlat, end_lonlat, color="red", marker="x", **kwargs):
        line_lon = [start_lonlat[0], end_lonlat[0]]
        line_lat = [start_lonlat[1], end_lonlat[1]]
        return ax.plot(line_lon, line_lat, color=color, marker=marker, transform=self.transform, zorder=20, **kwargs)


class RadarDisplay(Graph):
    """Alias of Graph kept as the recommended user-facing display class."""


class RadarMapDisplay(GraphMap):
    """Alias of GraphMap kept as the recommended user-facing map display class."""


def plot_xy(
    ax,
    x,
    y,
    data,
    cmap="CN_ref",
    bounds=np.arange(-5, 76, 5),
    cbar=True,
    orientation="vertical",
    cbar_ticks=None,
    cbar_ticklabels=None,
    clabel=None,
    extent_km=None,
    labels=None,
    reference_options=None,
    **kwargs
):
    _require_mpl_axes(ax)
    style = PlotStyle(
        cmap=cmap,
        levels=bounds,
        colorbar=ColorbarOptions(
            visible=cbar,
            orientation=orientation,
            label=clabel,
            ticks=cbar_ticks,
            ticklabels=cbar_ticklabels,
        ),
    )
    return plot_cartesian(
        ax,
        x,
        y,
        data,
        style,
        fig=ax.figure,
        extent_km=extent_km,
        labels=labels,
        reference_options=reference_options,
        **kwargs
    )


def add_rings(ax, rings, color="#5B5B5B", linestyle="-", linewidth=0.6, **kwargs):
    theta = np.linspace(0, 2 * np.pi, 200)
    for radius in rings:
        x0 = radius * np.cos(theta)
        y0 = radius * np.sin(theta)
        ax.plot(x0, y0, linestyle=linestyle, linewidth=linewidth, color=color, **kwargs)
    for rad in np.arange(0, np.pi, np.pi / 6):
        ax.plot(
            [-1 * rings[-1] * np.sin(rad), rings[-1] * np.sin(rad)],
            [-1 * rings[-1] * np.cos(rad), rings[-1] * np.cos(rad)],
            linestyle=linestyle,
            linewidth=linewidth,
            color=color,
            **kwargs
        )
    return ax


def plot_az_ranges(
    ax,
    _range,
    azimuth,
    elevation,
    data,
    cmap="CN_ref",
    bounds=np.arange(-5, 76, 5),
    cbar=True,
    orientation="vertical",
    cbar_ticks=None,
    cbar_ticklabels=None,
    clabel=None,
    extent_km=None,
    labels=None,
    reference_options=None,
    **kwargs
):
    _require_mpl_axes(ax)
    x, y, _ = antenna_vectors_to_cartesian(_range, azimuth, elevation, edges=True)
    return plot_xy(
        ax,
        x,
        y,
        data,
        cmap=cmap,
        bounds=bounds,
        cbar=cbar,
        orientation=orientation,
        cbar_ticks=cbar_ticks,
        cbar_ticklabels=cbar_ticklabels,
        clabel=clabel,
        extent_km=extent_km,
        labels=labels,
        reference_options=reference_options,
        **kwargs
    )
