# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from ..configure.default_config import CINRAD_COLORMAP, CINRAD_field_bins,\
    CINRAD_field_normvar, CINRAD_field_mapping
import pandas as pd
from ..core.transforms import geographic_to_cartesian_aeqd, cartesian_to_geographic_aeqd

class VerticalSection(object):

    def __init__(self, NuistRadar):
        self.NRadar = NuistRadar

    def RHI(self, azimuth, field_name, height=(0,18), title=None, clabel=None, continuously=False):
        """
        :param azimuth:要剖面的方位角 degree
        :param field_name:要剖的数据场
        :param title: 剖面图的title
        :param clabel: colorbar的title
        :param continuously: 是否使用连续的colorbar
        :return:
        """
        mesh_xy, mesh_z, grid_field = self.NRadar.get_RHI_data(azimuth, field_name)
        assert self.NRadar.scan_info.scan_type.values == "ppi", "now only support ppi scan!"
        if field_name == "V":
            vmax = self.NRadar.scan_info.nyquist_velocity[0].values
            vmin = -1 * vmax
        elif CINRAD_field_normvar[CINRAD_field_mapping[field_name]] == -1:
            vmax = np.nanmax(self.NRadar.fields[0][field_name])
            vmin = np.nanmin(self.NRadar.fields[0][field_name])
        else:
            vmin, vmax = CINRAD_field_normvar[CINRAD_field_mapping[field_name]]

        cmap_bins = CINRAD_field_bins[CINRAD_field_mapping[field_name]]
        cmap = CINRAD_COLORMAP[CINRAD_field_mapping[field_name]]

        if title is None:
            title = pd.to_datetime(self.NRadar.fields[0].time[0].item()).strftime("UTC %Y-%m-%d %H:%M:%S")
        if clabel is None:
            clabel = CINRAD_field_mapping[field_name]

        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.3, 0.8, 0.6])
        cax = fig.add_axes([0.1, 0.1, 0.8, 0.06])

        return VerticalSection.SectionPlot_VCS(fig, ax, cax, mesh_xy, mesh_z, grid_field, height, title=title,
                                           normvar=(vmin, vmax), cmap=cmap, cmap_bins=cmap_bins, clabel=clabel,
                                           continuously=continuously)

    def section(self, start_point, end_point, field_name, height=(0, 18), title=None, clabel=None, continuously=False):
        """
        :param start_point: (start_x, start_y), units:km
        :param end_point: (end_x, end_y), units:km
        :param field_name: field names eg: dBZ
        :param height: height to show (min_height, max_height), units:km
        :param title:  the title to show top of graph
        :param clabel: the title of cbar
        :param continuously: if True, using continuously colormap
        :return:
        """
        assert len(start_point) == 2, "start pionts should be (start_x, start_y), units:km!"
        assert len(end_point) == 2, "end pionts should be (end_x, end_y), units:km!"

        start_point = (start_point[0]*1000., start_point[1]*1000)
        end_point = (end_point[0] * 1000., end_point[1] * 1000)
        mesh_xy, mesh_z, grid_field = self.NRadar.get_vcs_data(start_point, end_point, field_name)
        assert self.NRadar.scan_info.scan_type.values == "ppi", "now only support ppi scan!"
        if field_name == "V":
            vmax = self.NRadar.scan_info.nyquist_velocity[0].values
            vmin = -1 * vmax
        elif CINRAD_field_normvar[CINRAD_field_mapping[field_name]] == -1:
            vmax = np.nanmax(self.NRadar.fields[0][field_name])
            vmin = np.nanmin(self.NRadar.fields[0][field_name])
        else:
            vmin, vmax = CINRAD_field_normvar[CINRAD_field_mapping[field_name]]

        cmap_bins = CINRAD_field_bins[CINRAD_field_mapping[field_name]]
        cmap = CINRAD_COLORMAP[CINRAD_field_mapping[field_name]]

        if title is None:
            title = pd.to_datetime(self.NRadar.fields[0].time[0].item()).strftime("UTC %Y-%m-%d %H:%M:%S")
        if clabel is None:
            clabel = CINRAD_field_mapping[field_name]

        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.3, 0.8, 0.6])
        cax = fig.add_axes([0.1, 0.1, 0.8, 0.06])
        return VerticalSection.SectionPlot_VCS(fig, ax, cax, mesh_xy, mesh_z, grid_field, height=height, title=title,
                                           normvar=(vmin, vmax), cmap=cmap, cmap_bins=cmap_bins, clabel=clabel,
                                           continuously=continuously)

    def section_map(self, start_lonlat, end_lonlat, field_name, height=(0,18), title=None, \
                     orient="horizontal", label=None, clabel=None, continuously=False):
        """
        :param start_latlon: (start_lat, start_lon), units: degrees
        :param end_latlon: (end_lat, end_lon), units: degrees
        :param field_name: field names eg: dBZ
        :param height: height to show (min_height, max_height), units:km
        :param title: the title to show top of graph
        :param clabel: the title of cbar
        :param continuously:
        :return:
        """
        assert np.abs(start_lonlat[1]) <= 90, "check, (lon, lat), lon first, lat second!"

        if field_name == "V":
            vmax = self.NRadar.scan_info.nyquist_velocity[0].values
            vmin = -1 * vmax
        elif CINRAD_field_normvar[CINRAD_field_mapping[field_name]] == -1:
            vmax = np.nanmax(self.NRadar.fields[0][field_name])
            vmin = np.nanmin(self.NRadar.fields[0][field_name])
        else:
            vmin, vmax = CINRAD_field_normvar[CINRAD_field_mapping[field_name]]
        cmap_bins = CINRAD_field_bins[CINRAD_field_mapping[field_name]]
        cmap = CINRAD_COLORMAP[CINRAD_field_mapping[field_name]]
        if title is None:
            title = pd.to_datetime(self.NRadar.fields[0].time[0].item()).strftime("UTC %Y-%m-%d %H:%M:%S")
        if clabel is None:
            clabel = CINRAD_field_mapping[field_name]
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.35, 0.8, 0.55])
        cax = fig.add_axes([0.1, 0.1, 0.8, 0.06])
        return VerticalSection.SectionPlot_VCS_map(fig, ax, cax, start_lonlat, end_lonlat, field_name=field_name,
                                                   NRadar=self.NRadar, height=height, title=title,normvar=(vmin, vmax),
                                                   cmap=cmap, cmap_bins=cmap_bins, clabel=clabel, continuously=continuously)


    @staticmethod
    def GUI_section(fig, ax, cax, NRadar, start_point, end_point, field_name, title=None, clabel=None, continuously=False):
        """
        :param start_point:剖面的起点(x, y)
        :param end_point:剖面的结束点(x, y)
        :param field_name:要剖的数据场
        :param title: 剖面图的title
        :param clabel: colorbar的title
        :param continuously: 是否使用连续的colorbar
        :return:
        """
        mesh_xy, mesh_z, grid_field = NRadar.get_vcs_data(start_point, end_point, field_name)
        assert NRadar.scan_info.scan_type.values == "ppi", "now only support ppi scan!"
        if field_name == "V":
            vmax = NRadar.scan_info.nyquist_velocity[0].values
            vmin = -1 * vmax
        elif CINRAD_field_normvar[CINRAD_field_mapping[field_name]] == -1:
            vmax = np.nanmax(NRadar.fields[0][field_name])
            vmin = np.nanmin(NRadar.fields[0][field_name])
        else:
            vmin, vmax = CINRAD_field_normvar[CINRAD_field_mapping[field_name]]

        cmap_bins = CINRAD_field_bins[CINRAD_field_mapping[field_name]]
        cmap = CINRAD_COLORMAP[CINRAD_field_mapping[field_name]]

        if title is None:
            title = pd.to_datetime(NRadar.fields[0].time[0].item()).strftime("UTC %Y-%m-%d %H:%M:%S")
        if clabel is None:
            clabel = CINRAD_field_mapping[field_name]
        return VerticalSection.SectionPlot_VCS(fig, ax, cax, mesh_xy, mesh_z, grid_field, title=title,
                                           normvar=(vmin, vmax), cmap=cmap, cmap_bins=cmap_bins, clabel=clabel,
                                           continuously=continuously)

    @staticmethod
    def GUI_section_map(fig, ax, cax, NRadar, start_lonlat, end_lonlat, field_name, title=None, clabel=None,
                    continuously=False):
        """
        :param start_point:剖面的起点(x, y)
        :param end_point:剖面的结束点(x, y)
        :param field_name:要剖的数据场
        :param title: 剖面图的title
        :param clabel: colorbar的title
        :param continuously: 是否使用连续的colorbar
        :return:
        """
        if field_name == "V":
            vmax = NRadar.scan_info.nyquist_velocity[0].values
            vmin = -1 * vmax
        elif CINRAD_field_normvar[CINRAD_field_mapping[field_name]] == -1:
            vmax = np.nanmax(NRadar.fields[0][field_name])
            vmin = np.nanmin(NRadar.fields[0][field_name])
        else:
            vmin, vmax = CINRAD_field_normvar[CINRAD_field_mapping[field_name]]

        cmap_bins = CINRAD_field_bins[CINRAD_field_mapping[field_name]]
        cmap = CINRAD_COLORMAP[CINRAD_field_mapping[field_name]]

        if title is None:
            title = pd.to_datetime(NRadar.fields[0].time[0].item()).strftime("UTC %Y-%m-%d %H:%M:%S")
        if clabel is None:
            clabel = CINRAD_field_mapping[field_name]
        return VerticalSection.SectionPlot_VCS_map(fig, ax, cax, start_lonlat, end_lonlat, field_name=field_name, NRadar=NRadar,
                                               title=title, normvar=(vmin, vmax), cmap=cmap, cmap_bins=cmap_bins, clabel=clabel,
                                               continuously=continuously)

    @staticmethod
    def SectionPlot_VCS(fig, ax, cx, mesh_xy, mesh_z, field_data, height=(0, 18), title=None, normvar=None, cmap=None, \
                    cmap_bins=16, orient="horizontal", label=None, clabel=None, continuously=False):
        """
        画剖面图像
        :param fig:matplotlib的fig对象
        :param ax:axes对象
        :param cx:cax对象
        :param mesh_xy:剖面的线
        :param mesh_z:垂直的线
        :param title:标题
        :param normvar:最小最大值
        :param cmap: str
        :param cmap_bins: cmap要分的色段
        :param orient:map的方向, 默认是水平
        :param label:x和y轴的label
        :param clabel:cbar的label
        :param continuously:连续色标
        :return:
        """
        if label is None:
            xlabel = "Distance From Section Start (Uints:km)"
            ylabel = "Height (Uints:km)"
        else:
            xlabel, ylable = label
        if normvar is None:
            vmin = -5
            vmax = 75
        else:
            vmin, vmax = normvar
        if title is None:
            title = ""
        if clabel is None:
            clabel = ""
        if continuously:
            cmap_bins = 256
        min_h, max_h = height
        cmaps = plt.get_cmap(cmap)
        levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmaps.N, clip=True)
        for isweep, _ in enumerate(mesh_xy):
            gci = ax.pcolormesh(mesh_xy[isweep] / 1000., mesh_z[isweep] / 1000., field_data[isweep], cmap=cmaps, norm=norm)
        ax.set_xlabel(xlabel, fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)
        ax.set_title(title, fontsize=16)
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.set_ylim([min_h, max_h])
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.tick_params(axis="y", which="both", direction='in')
        ax.tick_params(axis="x", which="both", direction='in')
        if continuously:
            cbar = fig.colorbar(mappable=gci, cax=cx, orientation=orient)
        else:
            cbar = fig.colorbar(mappable=gci, cax=cx, orientation=orient, ticks=levels)
            cbar.set_ticklabels(VerticalSection._FixTicks(levels))
        cbar.set_label(clabel, fontsize=14)
        ax.grid(True, zorder=15)
        plt.style.use('default')

    @staticmethod
    def SectionPlot_VCS_map(fig, ax, cx, start_lonlat, end_lonlat, field_name, NRadar, height=(0, 18),
                            title=None, normvar=None, cmap=None, cmap_bins=16, orient="horizontal",
                            label=None, clabel=None, continuously=False):
        """

        :param fig: matplotlib figure
        :param ax: matplotlib axes
        :param cx: matplotlib cbar axes
        :param start_lonlat: (start_lon, start_lat), units:degree
        :param end_lonlat: (end_lon, end_lat), units:degree
        :param field_name: for exmaple "dBZ", "V"
        :param NRadar:
        :param height:(MIN_H, MAX_H)， units:km
        :param title:
        :param normvar:
        :param cmap:
        :param cmap_bins:
        :param orient:
        :param label:
        :param clabel:
        :param continuously:
        :return:
        """
        if label is None:
            xlabel = "Latitude, Longitude"
            ylabel = "Height (km) "
        else:
            xlabel, ylable = label
        if title is None:
            title = ""
        if clabel is None:
            clabel = ""
        if continuously:
            cmap_bins = 256
        if normvar is None:
            vmin = -5
            vmax = 75
        else:
            vmin, vmax = normvar
        min_h, max_h = height
        cmaps = plt.get_cmap(cmap)
        lat_0 = NRadar.scan_info.latitude.values
        lon_0 = NRadar.scan_info.longitude.values
        levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmaps.N, clip=True)
        start_x, start_y = geographic_to_cartesian_aeqd(lat=start_lonlat[1], lon=start_lonlat[0], lat_0=lat_0, lon_0=lon_0)
        end_x, end_y = geographic_to_cartesian_aeqd(lat=end_lonlat[1], lon=end_lonlat[0], lat_0=lat_0, lon_0=lon_0)
        mesh_xy, mesh_z, field_data = NRadar.get_vcs_data((start_x[0], start_y[0]), (end_x[0], end_y[0]), field_name)
        for isweep, _ in enumerate(mesh_xy):
            gci = ax.pcolormesh(mesh_xy[isweep] / 1000., mesh_z[isweep] / 1000., field_data[isweep], cmap=cmaps, norm=norm)
        ax.set_xlabel(xlabel, fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)
        ax.set_title(title, fontsize=16)
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.tick_params(axis="y", which="both", direction='in')
        ax.tick_params(axis="x", which="both", direction='in')
        ax.set_ylim([min_h, max_h])
        ax.tick_params(axis='both', which='major', labelsize=12)
        xticks_data = ax.get_xticks()
        x_points_tk, y_points_tk = VerticalSection.get_points_from_ranges((start_x[0]/1000., start_y[0]/1000),\
                                                                          (end_x[0]/1000, end_y[0]/1000), xticks_data)
        lon_point, lat_point = cartesian_to_geographic_aeqd(x_points_tk*1000., y_points_tk*1000., \
                                                            lon_0=lon_0, lat_0=lat_0) #to meters

        ax.set_xticklabels(["(%.2f, %.2f)"%(lon_point[i], lat_point[i])\
                            for i, _ in enumerate(xticks_data)], rotation=15, fontsize=10)
        if continuously:
            cbar = fig.colorbar(mappable=gci, cax=cx, orientation=orient)
        else:
            cbar = fig.colorbar(mappable=gci, cax=cx, orientation=orient, ticks=levels)
            cbar.set_ticklabels(VerticalSection._FixTicks(levels))
        cbar.set_label(clabel, fontsize=14)
        ax.grid(True, zorder=15)
        plt.style.use('default')

    @staticmethod
    def _FixTicks(ticks):
        """修改ticks的显示label"""
        if (ticks % 1).sum() == 0:
            temp = ["%2.f" % i for i in ticks]
        else:
            temp = ["%.2f" % i for i in ticks]
        return temp

    @staticmethod
    def get_points_from_ranges(start_points, end_points, distance_from_start):
        """
        根据距离起始点的距离，计算该点的坐标
        :param start_points: units:km
        :param end_points:  units:km
        :param distance_from_start: units:km
        :return:
        """
        start_x, start_y = start_points
        end_x, end_y = end_points

        rgs = np.sqrt((start_x-end_x)**2 + (start_y-end_y)**2)
        return distance_from_start/rgs*(end_x - start_x) + start_x, distance_from_start/rgs*(end_y - start_y) + start_y
