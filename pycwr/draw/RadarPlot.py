import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from ..configure.default_config import CINRAD_COLORMAP, CINRAD_field_bins, \
    CINRAD_field_normvar, CINRAD_field_mapping
import numpy as np
from ..configure.location_config import CN_shp_info
import cartopy.feature as cfeature
from ..core.transforms import geographic_to_cartesian_aeqd, cartesian_to_geographic_aeqd, antenna_vectors_to_cartesian
from .VerticalSectionPlot import VerticalSection
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy, matplotlib

class Graph(object):
    """Improved mapping function, cartesian coords, recommended"""
    def __init__(self, NRadar):
        self.Radar = NRadar

    def plot_ppi(self, ax, sweep_num, field_name, cmap=None, min_max=None, cmap_bins=None, cbar=True,
                 orientation="vertical",cbar_ticks=None, cbar_ticklabels=None, clabel=None, **kwargs):
        """
        :param ax: axes.Axes object or array of Axes objects., eg: fig, ax = plt.subplots
        :param sweep_num: The sweep_num volume scan to draw, from 0 start!
        :param field_name: field dict to select data, eg: "dBZ" "V"
        :param cmap: str or Colormap, optional, A Colormap instance or registered colormap name. to see cm.py!
        :param min_max: The colorbar range(vmin, vmax). If None, suitable min/max values are automatically chosen by min max of data!
        :param cmap_bins: bins of colormaps
        :param cbar: if True, plot with colorbar, else not!
        :param orientation: vertical or horizontal, if cbar is True , this is vaild!, colorbar oriention!
        :param cbar_ticks: Set the locations of the tick marks from sequence ticks
        :param cbar_ticklabels: Set the text values of the tick labels.
        :param kwargs: other arguments for pcolormesh!
        :return:
        """
        assert isinstance(ax, matplotlib.axes._axes.Axes), "axes should be matplotlib axes not cartopy axes!"
        if field_name == "V":
            vmax = self.Radar.scan_info.nyquist_velocity[sweep_num].values
            vmin = -1 * vmax
        elif min_max is not None:
            vmin, vmax = min_max
        elif CINRAD_field_normvar[CINRAD_field_mapping[field_name]] == -1:
            vmax = np.nanmax(self.Radar.fields[sweep_num][field_name])
            vmin = np.nanmin(self.Radar.fields[sweep_num][field_name])
        else:
            vmin, vmax = CINRAD_field_normvar[CINRAD_field_mapping[field_name]]
        if cmap is None:
            cmap = CINRAD_COLORMAP[CINRAD_field_mapping[field_name]]
        if cmap_bins is None:
            cmap_bins = CINRAD_field_bins[CINRAD_field_mapping[field_name]]
        ax.set_aspect("equal")
        radar_data = self.Radar.fields[sweep_num][field_name]
        x, y = radar_data.x, radar_data.y
        cmaps = plt.get_cmap(cmap)
        levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmaps.N, clip=True)
        gci = ax.pcolormesh(x / 1000., y / 1000., radar_data, cmap=cmaps, \
                            zorder=0, norm=norm, shading='auto', **kwargs)
        if cbar:
            cb=plt.colorbar(mappable=gci, ax=ax, orientation=orientation)
            if cbar_ticks is None:
                ticks = levels
            else:
                ticks = cbar_ticks
            cb.set_ticks(ticks)

        if cbar_ticklabels is not None:
            if orientation == "vertical":
                cb.ax.set_yticklabels(cbar_ticklabels)
            else:
                cb.ax.set_xticklabels(cbar_ticklabels)

        if clabel is not None:
            cb.set_label(clabel)
        return gci

    def plot_rhi(self, ax, sweep_num, field_name, cmap=None, min_max=None,cmap_bins=None, cbar=True,
                 orientation="vertical",cbar_ticks=None, cbar_ticklabels=None, clabel=None, **kwargs):

        assert isinstance(ax, matplotlib.axes._axes.Axes), "axes should be matplotlib axes not cartopy axes!"
        if field_name == "V":
            vmax = self.Radar.scan_info.nyquist_velocity[0].values
            vmin = -1 * vmax
        elif min_max is not None:
            vmin, vmax = min_max
        elif CINRAD_field_normvar[CINRAD_field_mapping[field_name]] == -1:
            vmax = np.nanmax(self.Radar.fields[0][field_name])
            vmin = np.nanmin(self.Radar.fields[0][field_name])
        else:
            vmin, vmax = CINRAD_field_normvar[CINRAD_field_mapping[field_name]]
        if cmap is None:
            cmap = CINRAD_COLORMAP[CINRAD_field_mapping[field_name]]
        if cmap_bins is None:
            cmap_bins = CINRAD_field_bins[CINRAD_field_mapping[field_name]]

        mesh_xy = (self.Radar.fields[sweep_num].x ** 2 + self.Radar.fields[sweep_num].y ** 2) ** 0.5
        mesh_z = self.Radar.fields[sweep_num].z
        field_data = self.Radar.fields[sweep_num][field_name]
        cmaps = plt.get_cmap(cmap)
        levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmaps.N, clip=True)

        gci = ax.pcolormesh(mesh_xy/1000., mesh_z/1000., field_data, cmap=cmaps,
                                norm=norm, shading='auto', **kwargs)
        if cbar:
            cb = plt.colorbar(mappable=gci, ax=ax, orientation=orientation)
            if cbar_ticks is None:
                ticks = levels
            else:
                ticks = cbar_ticks
            cb.set_ticks(ticks)

        if cbar_ticklabels is not None:
            if orientation == "vertical":
                cb.ax.set_yticklabels(cbar_ticklabels)
            else:
                cb.ax.set_xticklabels(cbar_ticklabels)
        if clabel is not None:
            cb.set_label(clabel)

        return gci

    def plot_vcs(self, ax, start_xy, end_xy, field_name, cmap=None, min_max=None,cmap_bins=None, cbar=True,
                 orientation="vertical",cbar_ticks=None, cbar_ticklabels=None, clabel=None, **kwargs):
        """
        :param ax: axes.Axes object or array of Axes objects., eg: fig, ax = plt.subplots
        :param start_xy: (start_x, start_y) units:km, VCS start position!
        :param end_xy: (end_x, end_y) units:km, VCS end position!
        :param field_name: field dict to select data, eg: "dBZ" "V"
        :param cmap: str or Colormap, optional, A Colormap instance or registered colormap name. to see cm.py!
        :param min_max: The colorbar range(vmin, vmax). If None, suitable min/max values are automatically chosen by min max of data!
        :param cmap_bins: bins of colormap
        :param cbar: bool, if True, plot with colorbar,
        :param orientation: vertical or horizontal, it is vaild when cbar is True
        :param cbar_ticks: Set the locations of the tick marks from sequence ticks
        :param cbar_ticklabels: Set the text values of the tick labels.
        :return:
        """
        assert isinstance(ax, matplotlib.axes._axes.Axes), "axes should be matplotlib axes not cartopy axes!"
        if field_name == "V":
            vmax = self.Radar.scan_info.nyquist_velocity[0].values
            vmin = -1 * vmax
        elif min_max is not None:
            vmin, vmax = min_max
        elif CINRAD_field_normvar[CINRAD_field_mapping[field_name]] == -1:
            vmax = np.nanmax(self.Radar.fields[0][field_name])
            vmin = np.nanmin(self.Radar.fields[0][field_name])
        else:
            vmin, vmax = CINRAD_field_normvar[CINRAD_field_mapping[field_name]]
        if cmap is None:
            cmap = CINRAD_COLORMAP[CINRAD_field_mapping[field_name]]
        if cmap_bins is None:
            cmap_bins = CINRAD_field_bins[CINRAD_field_mapping[field_name]]

        start_point = (start_xy[0] * 1000., start_xy[1] * 1000) ##km to meters
        end_point = (end_xy[0] * 1000., end_xy[1] * 1000) ##km to meters
        mesh_xy, mesh_z, field_data = self.Radar.get_vcs_data(start_point, end_point, field_name)
        cmaps = plt.get_cmap(cmap)
        levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmaps.N, clip=True)
        for isweep, _ in enumerate(mesh_xy):
            gci = ax.pcolormesh(mesh_xy[isweep] / 1000., mesh_z[isweep] / 1000., field_data[isweep], cmap=cmaps,
                                norm=norm, shading='auto', **kwargs)
        if cbar:
            cb = plt.colorbar(mappable=gci, ax=ax, orientation=orientation)
            if cbar_ticks is None:
                ticks = levels
            else:
                ticks = cbar_ticks
            cb.set_ticks(ticks)
            if clabel is not None:
                cb.set_label(clabel)

        if cbar_ticklabels is not None:
            if orientation == "vertical":
                cb.ax.set_yticklabels(cbar_ticklabels)
            else:
                cb.ax.set_xticklabels(cbar_ticklabels)

        return gci

    def plot_crf(self, ax, cmap=CINRAD_COLORMAP[CINRAD_field_mapping["dBZ"]],
                 min_max=CINRAD_field_normvar[CINRAD_field_mapping["dBZ"]], cmap_bins=CINRAD_field_bins[CINRAD_field_mapping["dBZ"]],
                 cbar=True, orientation="vertical",cbar_ticks=None, cbar_ticklabels=None, clabel=None, **kwargs):
        """
        显示组合反射率因子
        :param ax: axes.Axes object or array of Axes objects., eg: fig, ax = plt.subplots
        :param XRange: np.ndarray, 1d, units:meters
        :param YRange: np.ndarray, 1d, units:meters
        :param cmap: str or Colormap, optional, A Colormap instance or registered colormap name. to see cm.py!
        :param min_max: The colorbar range(vmin, vmax). If None, suitable min/max values are automatically chosen by min max of data!
        :param cmap_bins: bins of colormaps
        :param cbar: if True, plot with colorbar, else not!
        :param orientation: vertical or horizontal, if cbar is True , this is vaild!, colorbar oriention!
        :param kwargs: other arguments for pcolormesh!
        :param cbar_ticks: Set the locations of the tick marks from sequence ticks
        :param cbar_ticklabels: Set the text values of the tick labels.
        :return:
        """
        max_range = int(self.Radar.fields[0].range.max().values)
        XRange = np.arange(-1 * max_range, max_range+1, 1000.)
        YRange = XRange
        self.Radar.add_product_CR_xy(XRange, YRange)
        assert isinstance(ax, matplotlib.axes._axes.Axes), "axes should be matplotlib axes not cartopy axes!"

        vmin, vmax = min_max
        ax.set_aspect("equal")
        radar_data = self.Radar.product['CR'].values
        x, y = np.meshgrid(self.Radar.product['CR'].x_cr.values,
                           self.Radar.product['CR'].y_cr.values, indexing="ij")
        cmaps = plt.get_cmap(cmap)
        levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmaps.N, clip=True)
        gci = ax.pcolormesh(x / 1000., y / 1000., radar_data, cmap=cmaps, \
                            zorder=0, norm=norm, shading='auto', **kwargs)
        if cbar:
            cb = plt.colorbar(mappable=gci, ax=ax, orientation=orientation)
            if cbar_ticks is None:
                ticks = levels
            else:
                ticks = cbar_ticks
            cb.set_ticks(ticks)
            if clabel is not None:
                cb.set_label(clabel)

        if cbar_ticklabels is not None:
            if orientation == "vertical":
                cb.ax.set_yticklabels(cbar_ticklabels)
            else:
                cb.ax.set_xticklabels(cbar_ticklabels)

        return gci

    def plot_cappi(self, ax, level_height=3000, cmap=CINRAD_COLORMAP[CINRAD_field_mapping["dBZ"]],
                   min_max=CINRAD_field_normvar[CINRAD_field_mapping["dBZ"]], cmap_bins=CINRAD_field_bins[CINRAD_field_mapping["dBZ"]],
                   cbar=True, orientation="vertical", cbar_ticks=None, cbar_ticklabels=None, clabel=None, **kwargs):
        """
        显示CAPPI图像
        :param ax: axes.Axes object or array of Axes objects., eg: fig, ax = plt.subplots
        :param level_height: height of cappi， units:meters, default, 3000m
        :param cmap: str or Colormap, optional, A Colormap instance or registered colormap name. to see cm.py!
        :param min_max: The colorbar range(vmin, vmax). If None, suitable min/max values are automatically chosen by min max of data!
        :param cmap_bins: bins of colormaps
        :param cbar: if True, plot with colorbar, else not!
        :param orientation: vertical or horizontal, if cbar is True , this is vaild!, colorbar oriention!
        :param cbar_ticks: Set the locations of the tick marks from sequence ticks
        :param cbar_ticklabels: Set the text values of the tick labels.
        :param kwargs: other arguments for pcolormesh!
        :return:
        """
        max_range = int(self.Radar.fields[0].range.max().values)
        XRange = np.arange(-1 * max_range, max_range + 1, 1000.)
        YRange = XRange
        self.Radar.add_product_CAPPI_xy(XRange, YRange, level_height)
        assert isinstance(ax, matplotlib.axes._axes.Axes), "axes should be matplotlib axes not cartopy axes!"

        vmin, vmax = min_max
        ax.set_aspect("equal")
        radar_data = self.Radar.product["CAPPI_%d"%level_height].values
        x, y = np.meshgrid(self.Radar.product["CAPPI_%d"%level_height]['x_cappi_%d'%level_height].values,
                           self.Radar.product["CAPPI_%d"%level_height]['y_cappi_%d'%level_height].values, indexing="ij")
        cmaps = plt.get_cmap(cmap)
        levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmaps.N, clip=True)
        gci = ax.pcolormesh(x / 1000., y / 1000., radar_data, cmap=cmaps, \
                            zorder=0, norm=norm,shading='auto', **kwargs)
        if cbar:
            cb = plt.colorbar(mappable=gci, ax=ax, orientation=orientation)
            if cbar_ticks is None:
                ticks = levels
            else:
                ticks = cbar_ticks
            cb.set_ticks(ticks)
            if clabel is not None:
                cb.set_label(clabel)

        if cbar_ticklabels is not None:
            if orientation == "vertical":
                cb.ax.set_yticklabels(cbar_ticklabels)
            else:
                cb.ax.set_xticklabels(cbar_ticklabels)

        return gci

    def add_rings(self, ax, rings, color="#5B5B5B", linestyle='-', linewidth=0.6, **kwargs):
        """
        :param ax: axes.Axes object or array of Axes objects., eg: fig, ax = plt.subplots
        :param rings: distance from radar (units:km)
        :param color: line color for rings
        :param linestyle: linestyle for rings, {'-', '--', '-.', ':', '', (offset, on-off-seq), ...}
        :param linewidth: linewidth for rings , float
        :param kwargs:  other arguments for ax.plot
        :return:
        """
        theta = np.linspace(0, 2 * np.pi, 200)
        for i in rings:
            x0 = i * np.cos(theta)
            y0 = i * np.sin(theta)
            gci = ax.plot(x0, y0, linestyle=linestyle, linewidth=linewidth, color=color, **kwargs)
        for rad in np.arange(0, np.pi, np.pi / 6):
            gci = ax.plot([-1 * rings[-1] * np.sin(rad), rings[-1] * np.sin(rad)], \
                    [-1 * rings[-1] * np.cos(rad), rings[-1] * np.cos(rad)], \
                    linestyle=linestyle, linewidth=linewidth, color=color, **kwargs)
        return gci

    def add_lines(self, ax, start_xy, end_xy, color='red', marker='x', markersize=12, **kwargs):
        """
        :param ax: ax: axes.Axes object or array of Axes objects., eg: fig, ax = plt.subplots
        :param start_xy: (start_x, start_y) units:km, line start position, units:km
        :param end_xy: (end_x, end_y) units:km, line end position, units:km
        :param color: color for line
        :param marker: 	marker style for line marker
        :param markersize: float, for markersize
        :param kwargs: kwargs are used to specify other properties like a line label in plot!
        :return:
        """
        line_x = [start_xy[0], end_xy[0]]
        line_y = [start_xy[1], end_xy[1]]
        gci = ax.plot(line_x, line_y, color=color, marker=marker,markersize=markersize, **kwargs)
        return gci

class GraphMap(object):
    def __init__(self, NRadar, transform):
        """
        :param NRadar: NRadar object, read from basedata
        :param transform: The transform argument to plotting functions tells Cartopy what coordinate system your data are defined in.
        """
        self.Radar = NRadar
        self.transform = transform
    def plot_ppi_map(self, ax, sweep_num, field_name, extend=None, cmap=None, min_max=None,\
                 cmap_bins=None, cbar=True, orientation="vertical",cbar_ticks=None, cbar_ticklabels=None, \
                     clabel=None, **kwargs):
        """
        :param ax: cartopy.mpl.geoaxes.GeoAxesSubplot, it should get from cartopy, eg:plt.axes(projection=ccrs.PlateCarree())
        :param sweep_num:  The sweep_num volume scan to draw, from 0 start!
        :param field_name: field dict to select data, eg: "dBZ" "V"
        :param extend: (min_lon, max_lon, min_lat, max_lat), Latitude and longitude range, units:degrees
        :param cmap: str or Colormap, optional, A Colormap instance or registered colormap name. to see cm.py!
        :param min_max: The colorbar range(vmin, vmax). If None, suitable min/max values are automatically chosen by min max of data!
        :param cmap_bins:bins of cmaps, int
        :param cbar: bool, if True, plot with colorbar,
        :param orientation: vertical or horizontal, it is vaild when cbar is True
        :param cbar_ticks: Set the locations of the tick marks from sequence ticks
        :param cbar_ticklabels: Set the text values of the tick labels.
        :param kwargs: kwargs: other arguments for pcolormesh!
        :return:
        """
        assert isinstance(ax, cartopy.mpl.geoaxes.GeoAxesSubplot), "axes is not cartopy axes!"
        if field_name == "V":
            vmax = self.Radar.scan_info.nyquist_velocity[sweep_num].values
            vmin = -1 * vmax
        elif min_max is not None:
            vmin, vmax = min_max
        elif CINRAD_field_normvar[CINRAD_field_mapping[field_name]] == -1:
            vmax = np.nanmax(self.Radar.fields[sweep_num][field_name])
            vmin = np.nanmin(self.Radar.fields[sweep_num][field_name])
        else:
            vmin, vmax = CINRAD_field_normvar[CINRAD_field_mapping[field_name]]
        if cmap is None:
            cmap = CINRAD_COLORMAP[CINRAD_field_mapping[field_name]]
        if cmap_bins is None:
            cmap_bins = CINRAD_field_bins[CINRAD_field_mapping[field_name]]
        if extend is None:
            min_lon = np.min(self.Radar.fields[sweep_num].lon)
            max_lon = np.max(self.Radar.fields[sweep_num].lon)
            min_lat = np.min(self.Radar.fields[sweep_num].lat)
            max_lat = np.max(self.Radar.fields[sweep_num].lat)
        else:
            min_lon, max_lon, min_lat, max_lat = extend

        #ax.set_aspect("equal")
        radar_data = self.Radar.fields[sweep_num][field_name]
        lat, lon = radar_data.lat, radar_data.lon
        cmaps = plt.get_cmap(cmap)
        levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmaps.N, clip=True)

        pm = ax.pcolormesh(lon, lat, radar_data, transform=self.transform, cmap=cmap, norm=norm, zorder=4, **kwargs)

        ax.add_feature(cfeature.OCEAN.with_scale('50m'), zorder=0)
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', \
                                                    edgecolor='none', facecolor="white"), zorder=1)
        ax.add_feature(cfeature.LAKES.with_scale('50m'), zorder=2)
        ax.add_feature(cfeature.RIVERS.with_scale('50m'), zorder=3)

        ax.add_feature(cfeature.ShapelyFeature(CN_shp_info.geometries(), self.transform, \
                                               edgecolor='k', facecolor='none'), linewidth=0.5, \
                       linestyle='-', zorder=5, alpha=0.8)
        parallels = np.arange(int(min_lat), np.ceil(max_lat) + 1, 1)
        meridians = np.arange(int(min_lon), np.ceil(max_lon) + 1, 1)
        ax.set_xticks(meridians, crs=self.transform)
        ax.set_yticks(parallels, crs=self.transform)
        lon_formatter = LongitudeFormatter()
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        if cbar:
            cb = plt.colorbar(mappable=pm, ax=ax, orientation=orientation)
            if cbar_ticks is None:
                ticks = levels
            else:
                ticks = cbar_ticks
            cb.set_ticks(ticks)
            if clabel is not None:
                cb.set_label(clabel)

        if cbar_ticklabels is not None:
            if orientation == "vertical":
                cb.ax.set_yticklabels(cbar_ticklabels)
            else:
                cb.ax.set_xticklabels(cbar_ticklabels)

        return ax

    def plot_cappi_map(self,  ax, level_height, extend=None, cmap=CINRAD_COLORMAP[CINRAD_field_mapping['dBZ']],
                    min_max=CINRAD_field_normvar[CINRAD_field_mapping['dBZ']], cmap_bins=CINRAD_field_bins[CINRAD_field_mapping['dBZ']],
                     cbar=True, orientation="vertical",cbar_ticks=None, cbar_ticklabels=None, clabel=None, **kwargs):
        """
        显示CAPPI图像
        :param ax: cartopy.mpl.geoaxes.GeoAxesSubplot, it should get from cartopy, eg:plt.axes(projection=ccrs.PlateCarree())
        :param level_height: height of cappi， units:meters, default, 3000m
        :param extend: (min_lon, max_lon, min_lat, max_lat), Latitude and longitude range, units:degrees
        :param cmap: str or Colormap, optional, A Colormap instance or registered colormap name. to see cm.py!
        :param min_max: The colorbar range(vmin, vmax). If None, suitable min/max values are automatically chosen by min max of data!
        :param cmap_bins:bins of cmaps, int
        :param cbar: bool, if True, plot with colorbar,
        :param orientation: vertical or horizontal, it is vaild when cbar is True
        :param cbar_ticks: Set the locations of the tick marks from sequence ticks
        :param cbar_ticklabels: Set the text values of the tick labels.
        :param kwargs: kwargs: other arguments for pcolormesh!
        :return:
        """
        assert isinstance(ax, cartopy.mpl.geoaxes.GeoAxesSubplot), "axes is not cartopy axes!"
        vmin, vmax = min_max

        if extend is None:
            min_lon = np.min(self.Radar.fields[0].lon)
            max_lon = np.max(self.Radar.fields[0].lon)
            min_lat = np.min(self.Radar.fields[0].lat)
            max_lat = np.max(self.Radar.fields[0].lat)
        else:
            min_lon, max_lon, min_lat, max_lat = extend
        XLON = np.arange(min_lon, max_lon, 0.01)
        YLAT = np.arange(min_lat, max_lat, 0.01)
        # ax.set_aspect("equal")
        self.Radar.add_product_CAPPI_lonlat(XLON, YLAT, level_height)
        radar_data = self.Radar.product["CAPPI_geo_%d" % level_height].values
        lon, lat = np.meshgrid(XLON, YLAT, indexing="ij")
        cmaps = plt.get_cmap(cmap)
        levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmaps.N, clip=True)
        pm = ax.pcolormesh(lon, lat, radar_data, transform=self.transform, cmap=cmap, norm=norm, zorder=4, **kwargs)
        # ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=self.transform)
        ax.add_feature(cfeature.OCEAN.with_scale('50m'), zorder=0)
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', \
                                                    edgecolor='none', facecolor="white"), zorder=1)
        ax.add_feature(cfeature.LAKES.with_scale('50m'), zorder=2)
        ax.add_feature(cfeature.RIVERS.with_scale('50m'), zorder=3)

        ax.add_feature(cfeature.ShapelyFeature(CN_shp_info.geometries(), self.transform, \
                                               edgecolor='k', facecolor='none'), linewidth=0.5, \
                       linestyle='-', zorder=5, alpha=0.8)
        parallels = np.arange(int(min_lat), np.ceil(max_lat) + 1, 1)
        meridians = np.arange(int(min_lon), np.ceil(max_lon) + 1, 1)
        ax.set_xticks(meridians, crs=self.transform)
        ax.set_yticks(parallels, crs=self.transform)
        lon_formatter = LongitudeFormatter()
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        if cbar:
            cb = plt.colorbar(mappable=pm, ax=ax, orientation=orientation)
            if cbar_ticks is None:
                ticks = levels
            else:
                ticks = cbar_ticks
            cb.set_ticks(ticks)
            if clabel is not None:
                cb.set_label(clabel)

        if cbar_ticklabels is not None:
            if orientation == "vertical":
                cb.ax.set_yticklabels(cbar_ticklabels)
            else:
                cb.ax.set_xticklabels(cbar_ticklabels)

        return ax

    def plot_crf_map(self, ax, extend=None, cmap=CINRAD_COLORMAP[CINRAD_field_mapping['dBZ']],
                    min_max=CINRAD_field_normvar[CINRAD_field_mapping['dBZ']], cmap_bins=CINRAD_field_bins[CINRAD_field_mapping['dBZ']],
                     cbar=True, orientation="vertical",cbar_ticks=None, cbar_ticklabels=None, clabel=None, **kwargs):
        """
        显示组合反射率因子
        :param ax: cartopy.mpl.geoaxes.GeoAxesSubplot, it should get from cartopy, eg:plt.axes(projection=ccrs.PlateCarree())
        :param extend: (min_lon, max_lon, min_lat, max_lat), Latitude and longitude range, units:degrees
        :param cmap: str or Colormap, optional, A Colormap instance or registered colormap name. to see cm.py!
        :param min_max: The colorbar range(vmin, vmax). If None, suitable min/max values are automatically chosen by min max of data!
        :param cmap_bins:bins of cmaps, int
        :param cbar: bool, if True, plot with colorbar,
        :param orientation: vertical or horizontal, it is vaild when cbar is True
        :param cbar_ticks: Set the locations of the tick marks from sequence ticks
        :param cbar_ticklabels: Set the text values of the tick labels.
        :param kwargs: kwargs: other arguments for pcolormesh!
        :return:
        """
        assert isinstance(ax, cartopy.mpl.geoaxes.GeoAxesSubplot), "axes is not cartopy axes!"
        vmin, vmax = min_max

        if extend is None:
            min_lon = np.min(self.Radar.fields[0].lon)
            max_lon = np.max(self.Radar.fields[0].lon)
            min_lat = np.min(self.Radar.fields[0].lat)
            max_lat = np.max(self.Radar.fields[0].lat)
        else:
            min_lon, max_lon, min_lat, max_lat = extend
        XLON = np.arange(min_lon, max_lon, 0.01)
        YLAT = np.arange(min_lat, max_lat, 0.01)
        #ax.set_aspect("equal")
        self.Radar.add_product_CR_lonlat(XLON, YLAT)
        radar_data = self.Radar.product["CR_geo"].values
        lon, lat = np.meshgrid(XLON, YLAT, indexing="ij")
        cmaps = plt.get_cmap(cmap)
        levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmaps.N, clip=True)
        pm = ax.pcolormesh(lon, lat, radar_data, transform=self.transform, cmap=cmap, norm=norm, zorder=4, **kwargs)
        #ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=self.transform)
        ax.add_feature(cfeature.OCEAN.with_scale('50m'), zorder=0)
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', \
                                                    edgecolor='none', facecolor="white"), zorder=1)
        ax.add_feature(cfeature.LAKES.with_scale('50m'), zorder=2)
        ax.add_feature(cfeature.RIVERS.with_scale('50m'), zorder=3)

        ax.add_feature(cfeature.ShapelyFeature(CN_shp_info.geometries(), self.transform, \
                                               edgecolor='k', facecolor='none'), linewidth=0.5, \
                       linestyle='-', zorder=5, alpha=0.8)
        parallels = np.arange(int(min_lat), np.ceil(max_lat) + 1, 1)
        meridians = np.arange(int(min_lon), np.ceil(max_lon) + 1, 1)
        ax.set_xticks(meridians, crs=self.transform)
        ax.set_yticks(parallels, crs=self.transform)
        lon_formatter = LongitudeFormatter()
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        if cbar:
            cb = plt.colorbar(mappable=pm, ax=ax, orientation=orientation)
            if cbar_ticks is None:
                ticks = levels
            else:
                ticks = cbar_ticks
            cb.set_ticks(ticks)
            if clabel is not None:
                cb.set_label(clabel)

        if cbar_ticklabels is not None:
            if orientation == "vertical":
                cb.ax.set_yticklabels(cbar_ticklabels)
            else:
                cb.ax.set_xticklabels(cbar_ticklabels)

        return ax

    def plot_vcs_map(self, ax, start_lonlat, end_lonlat, field_name, cmap=None, min_max=None,\
                 cmap_bins=None, cbar=True, orientation="vertical", cbar_ticks=None, cbar_ticklabels=None,\
                     clabel=None, **kwargs):
        """
        :param ax: axes.Axes object or array of Axes objects., eg: fig, ax = plt.subplots
        :param start_lonlat:(startlon, startlat),  VCS start position!
        :param end_lonlat:(endlon, endlat), VCS end position!
        :param field_name: field dict to select data, eg: "dBZ" "V"
        :param cmap:str or Colormap, optional, A Colormap instance or registered colormap name. to see cm.py!
        :param min_max:The colorbar range(vmin, vmax). If None, suitable min/max values are automatically chosen by min max of data!
        :param cmap_bins:bins of cmaps, int
        :param cbar:bool, if True, plot with colorbar,
        :param orientation:vertical or horizontal, it is vaild when cbar is True
        :param cbar_ticks: Set the locations of the tick marks from sequence ticks
        :param cbar_ticklabels: Set the text values of the tick labels.
        :param kwargs: other arguments for pcolormesh!
        :return:
        """
        assert isinstance(ax, matplotlib.axes._axes.Axes), "axes should be matplotlib axes not cartopy axes!"
        if field_name == "V":
            vmax = self.Radar.scan_info.nyquist_velocity[0].values
            vmin = -1 * vmax
        elif min_max is not None:
            vmin, vmax = min_max
        elif CINRAD_field_normvar[CINRAD_field_mapping[field_name]] == -1:
            vmax = np.nanmax(self.Radar.fields[0][field_name])
            vmin = np.nanmin(self.Radar.fields[0][field_name])
        else:
            vmin, vmax = CINRAD_field_normvar[CINRAD_field_mapping[field_name]]
        if cmap is None:
            cmap = CINRAD_COLORMAP[CINRAD_field_mapping[field_name]]
        if cmap_bins is None:
            cmap_bins = CINRAD_field_bins[CINRAD_field_mapping[field_name]]

        cmaps = plt.get_cmap(cmap)
        levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmaps.N, clip=True)
        start_x, start_y = geographic_to_cartesian_aeqd(lat=start_lonlat[1], lon=start_lonlat[0],
                                                        lat_0=self.Radar.scan_info.latitude.values,
                                                        lon_0=self.Radar.scan_info.longitude.values)
        end_x, end_y = geographic_to_cartesian_aeqd(lat=end_lonlat[1], lon=end_lonlat[0],
                                                    lat_0=self.Radar.scan_info.latitude.values,
                                                    lon_0=self.Radar.scan_info.longitude.values)
        mesh_xy, mesh_z, field_data = self.Radar.get_vcs_data((start_x[0], start_y[0]), (end_x[0], end_y[0]), field_name)
        for isweep, _ in enumerate(mesh_xy):
            gci = ax.pcolormesh(mesh_xy[isweep] / 1000., mesh_z[isweep] / 1000., field_data[isweep], cmap=cmaps,
                                norm=norm, **kwargs)

        xticks_data = ax.get_xticks()
        x_points_tk, y_points_tk = VerticalSection.get_points_from_ranges((start_x[0] / 1000., start_y[0] / 1000),
                                                                          (end_x[0] / 1000, end_y[0] / 1000),
                                                                          xticks_data)
        lon_point, lat_point = cartesian_to_geographic_aeqd(x_points_tk * 1000., y_points_tk * 1000.,
                                                            lat_0=self.Radar.scan_info.latitude.values,
                                                            lon_0=self.Radar.scan_info.longitude.values)  # to meters
        ax.set_xticklabels(["(%.2f, %.2f)" % (lon_point[i], lat_point[i]) \
                            for i, _ in enumerate(xticks_data)], rotation=15, fontsize=10)
        if cbar:
            cb = plt.colorbar(mappable=gci, ax=ax, orientation=orientation)
            if cbar_ticks is None:
                ticks = levels
            else:
                ticks = cbar_ticks
            cb.set_ticks(ticks)
            if clabel is not None:
                cb.set_label(clabel)

        if cbar_ticklabels is not None:
            if orientation == "vertical":
                cb.ax.set_yticklabels(cbar_ticklabels)
            else:
                cb.ax.set_xticklabels(cbar_ticklabels)

        return gci

    def add_lines_map(self, ax, start_lonlat, end_lonlat, color='red', marker='x', **kwargs):
        """
        :param ax: cartopy.mpl.geoaxes.GeoAxesSubplot, it should get from cartopy, eg:plt.axes(projection=ccrs.PlateCarree())
        :param start_lonlat: (startlon, startlat),  line start position!
        :param end_lonlat: (endlon, endlat),  line end position!
        :param color: str, color for line
        :param marker: str, style of marker
        :param markersize:float, size of marker
        :param kwargs:
        :return:
        """
        line_lon = [start_lonlat[0], end_lonlat[0]]
        line_lat = [start_lonlat[1], end_lonlat[1]]
        gci = ax.plot(line_lon, line_lat, color=color, marker=marker, transform=self.transform,zorder=20, **kwargs)
        return gci

def plot_xy(ax, x, y, data, cmap="CN_ref", bounds=np.arange(-5, 76, 5), cbar=True, orientation="vertical",
                cbar_ticks=None, cbar_ticklabels=None, clabel=None, **kwargs):
    """
    :param ax: axes.Axes object or array of Axes objects., eg: fig, ax = plt.subplots
    :param x: mesh grid x for data units: m
    :param y: y: mesh grid y for data units: m
    :param data: radar data ,dims like x,y
    :param cmap: str or Colormap, optional, A Colormap instance or registered colormap name. to see cm.py!
    :param orientation: vertical or horizontal, if cbar is True , this is vaild!, colorbar oriention!
    :param cbar: if True, plot with colorbar, else not!
    :param bounds: Monotonically increasing sequence of boundaries
    :param cbar_ticks: Set the locations of the tick marks from sequence ticks
    :param cbar_ticklabels: Set the text values of the tick labels.
    :param kwargs:
    :return:
    """
    assert isinstance(ax, matplotlib.axes._axes.Axes), "axes should be matplotlib axes not cartopy axes!"

    ax.set_aspect("equal")
    cmaps = plt.get_cmap(cmap)
    norm = BoundaryNorm(bounds, ncolors=cmaps.N, clip=True)
    gci = ax.pcolormesh(x / 1000., y / 1000., data, cmap=cmaps, \
                        zorder=0, norm=norm, **kwargs)
    if cbar:
        cb = plt.colorbar(mappable=gci, ax=ax, orientation=orientation)
        if cbar_ticks is not None and cbar_ticklabels is not None:
            cb.set_ticks(cbar_ticks)
            if orientation == "vertical":
                cb.ax.set_yticklabels(cbar_ticklabels)
            else:
                cb.ax.set_xticklabels(cbar_ticklabels)
        else:
            cb.set_ticks(bounds)

        if clabel is not None:
            cb.set_label(clabel)
    return gci

def add_rings(ax, rings, color="#5B5B5B", linestyle='-', linewidth=0.6, **kwargs):
    """
    :param ax: axes.Axes object or array of Axes objects., eg: fig, ax = plt.subplots
    :param rings: distance from radar (units:km)
    :param color: line color for rings
    :param linestyle: linestyle for rings, {'-', '--', '-.', ':', '', (offset, on-off-seq), ...}
    :param linewidth: linewidth for rings , float
    :param kwargs:  other arguments for ax.plot
    :return:
    """
    theta = np.linspace(0, 2 * np.pi, 200)
    for i in rings:
        x0 = i * np.cos(theta)
        y0 = i * np.sin(theta)
        gci = ax.plot(x0, y0, linestyle=linestyle, linewidth=linewidth, color=color, **kwargs)
    for rad in np.arange(0, np.pi, np.pi / 6):
        gci = ax.plot([-1 * rings[-1] * np.sin(rad), rings[-1] * np.sin(rad)], \
                [-1 * rings[-1] * np.cos(rad), rings[-1] * np.cos(rad)], \
                linestyle=linestyle, linewidth=linewidth, color=color, **kwargs)
    return gci


def plot_az_ranges(ax, _range, azimuth, elevation, data, cmap="CN_ref", bounds=np.arange(-5,76,5), cbar=True,
                   orientation="vertical", cbar_ticks=None, cbar_ticklabels=None, **kwargs):
    """
    :param ax:axes.Axes object or array of Axes objects., eg: fig, ax = plt.subplots
    :param _range: data second dim's range, 1d, units:meters, numpy.ndarray
    :param azimuth: data first dim's azimuth, 1d, units:degeree, numpy.ndarray
    :param elevation: data first dim's elevation, 1d, units:degeree, numpy.ndarray
    :param data: radar data ,dims like x,y
    :param cmap: str or Colormap, optional, A Colormap instance or registered colormap name. to see cm.py!
    :param min_max: The colorbar range(vmin, vmax). If None, suitable min/max values are automatically chosen by min max of data!
    :param cmap_bins: bins of colormaps
    :param cbar: if True, plot with colorbar, else not!
    :param orientation: vertical or horizontal, if cbar is True , this is vaild!, colorbar oriention!
    :param kwargs: other arguments for pcolormesh!
    :return:
    """
    assert isinstance(ax, matplotlib.axes._axes.Axes), "axes should be matplotlib axes not cartopy axes!"
    x, y, z = antenna_vectors_to_cartesian(_range, azimuth, elevation, edges=True)
    return plot_xy(ax, x, y, data, cmap=cmap, bounds=bounds, cbar=cbar, orientation=orientation, cbar_ticks=cbar_ticks,
                   cbar_ticklabels=cbar_ticklabels, **kwargs)

def plot_lonlat_map(ax, lon, lat, data, transform, extend=None, cmap="CN_ref", bounds=np.arange(-5,76,5),
                    cbar=True, orientation="vertical",cbar_ticks=None, cbar_ticklabels=None,  clabel=None,\
                    **kwargs):
    """
    :param ax:cartopy.mpl.geoaxes.GeoAxesSubplot, it should get from cartopy, eg:plt.axes(projection=ccrs.PlateCarree())
    :param lon: lon mesh grid for data units: degree
    :param lat: lat mesh grid for data units: degree
    :param data: radar data ,dims like lat, lon
    :param transform: The transform argument to plotting functions tells Cartopy what coordinate system your data are defined in.
    :param extend: (min_lon, max_lon, min_lat, max_lat), Latitude and longitude range, units:degrees
    :param cmap: str or Colormap, optional, A Colormap instance or registered colormap name. to see cm.py!
    :param min_max: The colorbar range(vmin, vmax). If None, suitable min/max values are automatically chosen by min max of data!
    :param cmap_bins:bins of cmaps, int
    :param cbar: bool, if True, plot with colorbar,
    :param orientation: vertical or horizontal, it is vaild when cbar is True
    :param kwargs: kwargs: other arguments for pcolormesh!
    :return:  pcolor result
    """
    assert isinstance(ax, cartopy.mpl.geoaxes.GeoAxesSubplot), "axes is not cartopy axes!"
    if extend is None:
        min_lon = np.min(lon)
        max_lon = np.max(lon)
        min_lat = np.min(lat)
        max_lat = np.max(lat)
    else:
        min_lon, max_lon, min_lat, max_lat = extend

    ax.set_aspect("equal")
    cmaps = plt.get_cmap(cmap)
    norm = BoundaryNorm(bounds, ncolors=cmaps.N, clip=True)
    pm = ax.pcolormesh(lon, lat, data, transform=transform, cmap=cmap, norm=norm, zorder=4, **kwargs)
    ax.add_feature(cfeature.OCEAN.with_scale('50m'), zorder=0)
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', \
                                                edgecolor='none', facecolor="white"), zorder=1)
    ax.add_feature(cfeature.LAKES.with_scale('50m'), zorder=2)
    ax.add_feature(cfeature.RIVERS.with_scale('50m'), zorder=3)

    ax.add_feature(cfeature.ShapelyFeature(CN_shp_info.geometries(), transform, \
                                           edgecolor='k', facecolor='none'), linewidth=0.5, \
                   linestyle='-', zorder=5, alpha=0.8)
    parallels = np.arange(int(min_lat), np.ceil(max_lat) + 1, 1)
    meridians = np.arange(int(min_lon), np.ceil(max_lon) + 1, 1)
    ax.set_xticks(meridians, crs=transform)
    ax.set_yticks(parallels, crs=transform)
    lon_formatter = LongitudeFormatter()
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    if cbar:
        cb = plt.colorbar(mappable=pm, ax=ax, orientation=orientation)
        if cbar_ticks is None:
            ticks = bounds
        else:
            ticks = cbar_ticks
        cb.set_ticks(ticks)
        if clabel is not None:
            cb.set_label(clabel)

    if cbar_ticklabels is not None:
        if orientation == "vertical":
            cb.ax.set_yticklabels(cbar_ticklabels)
        else:
            cb.ax.set_xticklabels(cbar_ticklabels)
    return pm

def plot_az_ranges_map(ax, _range, azimuth, elevation, data, main_point, transform, extend=None, cmap="CN_ref",
                       bounds=np.arange(-5,76,5), cbar=True, orientation="vertical",cbar_ticks=None, cbar_ticklabels=None,
                       **kwargs):
    """
    plot radar data with map, using range, azimuth, elevation
    :param ax:cartopy.mpl.geoaxes.GeoAxesSubplot, it should get from cartopy, eg:plt.axes(projection=ccrs.PlateCarree())
    :param _range:data second dim's range, 1d, units:meters, numpy.ndarray
    :param azimuth:data first dim's azimuth, 1d, units:degeree, numpy.ndarray
    :param elevation:data first dim's elevation, 1d, units:degeree, numpy.ndarray
    :param data:radar data ,dims like lat, lon
    :param main_point: list, (lon_0, lat_0) of radar station, units:degree
    :param transform: The transform argument to plotting functions tells Cartopy what coordinate system your data are defined in.
    :param extend: (min_lon, max_lon, min_lat, max_lat), Latitude and longitude range, units:degrees
    :param cmap: str or Colormap, optional, A Colormap instance or registered colormap name. to see cm.py!
    :param min_max: The colorbar range(vmin, vmax). If None, suitable min/max values are automatically chosen by min max of data!
    :param cmap_bins:bins of cmaps, int
    :param cbar: bool, if True, plot with colorbar,
    :param orientation: vertical or horizontal, it is vaild when cbar is True
    :param kwargs: kwargs: other arguments for pcolormesh!
    :return:  pcolor result
    """
    assert isinstance(ax, cartopy.mpl.geoaxes.GeoAxesSubplot), "axes is not cartopy axes!"
    main_lon, main_lat = main_point
    x, y, z = antenna_vectors_to_cartesian(_range, azimuth, elevation, edges=True)
    lon, lat = cartesian_to_geographic_aeqd(x, y, main_lon, main_lat)
    return plot_lonlat_map(ax, lon, lat, data,transform, extend, cmap, bounds, cbar, orientation,cbar_ticks, cbar_ticklabels **kwargs)