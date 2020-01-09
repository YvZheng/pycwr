# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 22:19:11 2019

this script is created to plot Stereotypical radar  images

@author: zy
"""
import cartopy.feature as cfeature
from ..configure.location_config import CN_shp_info
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import math
import cartopy.crs as ccrs
import numpy as np
import pandas as pd
from ..core.transforms import cartesian_to_geographic_aeqd, antenna_vectors_to_cartesian_cwr
from ..configure.default_config import CINRAD_COLORMAP, CINRAD_field_bins, \
    CINRAD_field_normvar, CINRAD_field_mapping, DEFAULT_METADATA

class RadarGraphMap(object):

    def __init__(self, NuistRadar = None):
        self.NRadar = NuistRadar

    def plot(self, sweep, field_name, normvar=None, title=None, clabel=None, continuously=False):
        """
        带地图画图模块
        :param sweep:扫描的sweep序号，从0开始
        :param field_name:扫描的数据场
        :param normvar:最小最大值
        :param title:画图的title
        :param clabel:colorbar的title
        :param continuously: 是否使用连续的cmaps
        :return:
        """
        assert self.NRadar.scan_info.scan_type.values == "ppi", "now only support ppi scan!"
        if field_name == "V":
            vmax = self.NRadar.scan_info.nyquist_velocity[sweep].values
            vmin = -1 * vmax
        elif CINRAD_field_normvar[CINRAD_field_mapping[field_name]] == -1:
            vmax = np.nanmax(self.NRadar.fields[sweep][field_name])
            vmin = np.nanmin(self.NRadar.fields[sweep][field_name])
        else:
            vmin, vmax = CINRAD_field_normvar[CINRAD_field_mapping[field_name]]

        if normvar is not None:
            vmin, vmax = normvar
        cmap_bins = CINRAD_field_bins[CINRAD_field_mapping[field_name]]
        cmap = CINRAD_COLORMAP[CINRAD_field_mapping[field_name]]

        if title is None:
            title = pd.to_datetime(self.NRadar.fields[0].time[0].item()).strftime("UTC %Y-%m-%d %H:%M:%S")
            title = title + " Elevation : %.1f" % self.NRadar.scan_info.fixed_angle[sweep].values
        if clabel is None:
            clabel = CINRAD_field_mapping[field_name]  +\
                     " (%s)"%DEFAULT_METADATA[CINRAD_field_mapping[field_name]]['units']
        field = self.NRadar.fields[sweep][field_name]

        RadarGraphMap.simple_plot_ppi_map(radar_data=field, title=title, normvar=(vmin, vmax), cmap=cmap, cmap_bins=cmap_bins,\
                                          clabel=clabel, continuously=continuously)

    @staticmethod
    def GUI_plot(NRadar, fig, ax, cax, sweep, field_name,\
                 normvar=None, title=None, clabel=None, continuously=False):
        """
        带地图画图模块
        :param sweep:扫描的sweep序号，从0开始
        :param field_name:扫描的数据场
        :param normvar:最小最大值
        :param title:画图的title
        :param clabel:colorbar的title
        :param continuously: 是否使用连续的cmaps
        :return:
        """
        assert NRadar.scan_info.scan_type.values == "ppi", "now only support ppi scan!"
        if field_name == "V":
            vmax = NRadar.scan_info.nyquist_velocity[sweep].values
            vmin = -1 * vmax
        elif CINRAD_field_normvar[CINRAD_field_mapping[field_name]] == -1:
            vmax = np.nanmax(NRadar.fields[sweep][field_name])
            vmin = np.nanmin(NRadar.fields[sweep][field_name])
        else:
            vmin, vmax = CINRAD_field_normvar[CINRAD_field_mapping[field_name]]

        if normvar is not None:
            vmin, vmax = normvar

        cmap_bins = CINRAD_field_bins[CINRAD_field_mapping[field_name]]
        cmap = CINRAD_COLORMAP[CINRAD_field_mapping[field_name]]
        if title is None:
            title = pd.to_datetime(NRadar.fields[0].time[0].item()).strftime("UTC %Y-%m-%d %H:%M:%S")
            title = title + " Elevation : %.1f" % NRadar.scan_info.fixed_angle[sweep].values
        if clabel is None:
            clabel = CINRAD_field_mapping[field_name] + " (%s)"%DEFAULT_METADATA[CINRAD_field_mapping[field_name]]['units']

        field = NRadar.fields[sweep][field_name]
        assert hasattr(field, 'lat') and hasattr(field, 'lon'), "check NRadar's lat lon coords!"
        lon, lat = field.lon, field.lat
        RadarGraphMap.plot_ppi_map(fig, ax, cax, lon, lat, field, title=title,\
                                    normvar=(vmin, vmax), cmap = cmap, cmap_bins = cmap_bins,\
                                   clabel=clabel, continuously=continuously)

    @staticmethod
    def simple_plot_ppi_map(radar_data, _range=None, azimuth=None, elevation=None, main_piont=None,\
                            title=None, normvar=None, cmap=None, cmap_bins=16, extend=None, \
                            projection=ccrs.PlateCarree(), orient="vertical", \
                            clabel=None, continuously=False):
        """
        使用方位角 距离库作ppi图像
        :param _range: 距离库 units:m
        :param azimuth: 方位角 units:degree
        :param elevation: 仰角 units:degree
        :param radar_data: 数据产品 dBZ, V, W
        :param main_piont: 雷达中心站点的经纬度 lon, lat
        :param title: 画图的title
        :param normvar: cmap的最大最小值
        :param cmap: cmap, str
        :param cmap_bins: cmap分多少段
        :param extend: 画的图像的范围
        :param projection: 投影方式, 目前支持ccrs.PlateCarree()
        :param orient: colorbar的朝向，目前支持vertical
        :param clabel: colorbar的title
        :param continuously: 是否使用连续的cmaps
        :return:
        """
        if hasattr(elevation, "values"):
            elevation = elevation.values
        if hasattr(azimuth, "values"):
            azimuth = azimuth.values
        if hasattr(_range, "values"):
            _range = _range.values
        assert radar_data is not None, "radar_data should not be None!"
        if hasattr(radar_data, 'lat') and hasattr(radar_data, 'lon'):
            lon, lat = radar_data.lon, radar_data.lat
        else:
            assert main_piont is not None, "should input (station_lon, station_lat) as main_point!"
            main_lon, main_lat = main_piont
            x, y, z = antenna_vectors_to_cartesian_cwr(_range, azimuth, elevation, h=0)
            lon, lat = cartesian_to_geographic_aeqd(x, y, main_lon, main_lat)
        fig = plt.figure()
        ax = fig.add_axes([0.04, 0.1, 0.82, 0.82], projection=projection)
        cax = fig.add_axes([0.85, 0.1, 0.028, 0.82])
        ax.tick_params(axis="y", which="both", direction='in')
        ax.tick_params(axis="x", which="both", direction='in')
        return RadarGraphMap.plot_ppi_map(fig, ax, cax, lon, lat, radar_data, title, normvar, cmap, \
                 cmap_bins, extend, projection, orient, clabel, continuously)

    @staticmethod
    def simple_plot_ppi_xy_map(x, y, radar_data, main_piont, title=None, normvar=None, cmap=None, \
                 cmap_bins=16, extend=None, projection=ccrs.PlateCarree(), orient="vertical", \
                clabel=None, continuously=False):

        main_lon, main_lat = main_piont
        lon, lat = cartesian_to_geographic_aeqd(x, y, main_lon, main_lat)
        fig = plt.figure()
        ax = fig.add_axes([0.04, 0.1, 0.82, 0.82], projection=projection)
        cax = fig.add_axes([0.85, 0.1, 0.028, 0.82])
        ax.tick_params(axis="y", which="both", direction='in')
        ax.tick_params(axis="x", which="both", direction='in')
        return RadarGraphMap.plot_ppi_map(fig, ax, cax, lon, lat, radar_data, title, normvar, cmap, \
                                          cmap_bins, extend, projection, orient, clabel, continuously)

    @staticmethod
    def plot_ppi_map(fig, ax, cx, lon, lat, radar_data, title=None, normvar=None, cmap=None, \
                 cmap_bins=16, extend=None, projection=ccrs.PlateCarree(), orient="vertical", \
                clabel=None, continuously=False):
        """
        显示带底图的雷达图像
        :param fig:matplotlib的fig对象
        :param ax: axes对象
        :param cx:colorbar axes对象
        :param lon: 经度网格
        :param lat: 纬度网格
        :param radar_data: 雷达数据产品
        :param title: 画图的title
        :param normvar: cmap的最小最大值
        :param cmap:cmap的类型
        :param cmap_bins: cmap分多少段
        :param extend: 图像的lat lon的范围（minlon, maxlon, minlat, maxlat）
        :param projection: 地图的投影方式
        :param orient: colorbar的朝向
        :param clabel: colorbar的title
        :param continuously:使用采用连续的colorbar
        :return:
        """
        if normvar is None:
            vmin = -5
            vmax = 75
        else:
            vmin, vmax = normvar
        if extend is None:
            min_lon = np.min(lon)
            max_lon = np.max(lon)
            min_lat = np.min(lat)
            max_lat = np.max(lat)
        else:
            min_lon, max_lon, min_lat, max_lat = extend
        if title is None:
            title = ""
        if clabel is None:
            clabel = ""
        if continuously:
            cmap_bins = 256

        cmaps = plt.get_cmap(cmap)
        levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmaps.N, clip=True)
        pm = ax.pcolormesh(lon, lat, radar_data, transform=projection, cmap=cmap, norm = norm, zorder=4)

        ax.add_feature(cfeature.OCEAN.with_scale('50m'), zorder=0)
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m',\
                                                    edgecolor='none', facecolor="white"), zorder=1)
        ax.add_feature(cfeature.LAKES.with_scale('50m'), zorder=2)
        ax.add_feature(cfeature.RIVERS.with_scale('50m'), zorder=3)
        ax.set_extent([min_lon, max_lon, min_lat, max_lat], projection)
        ax.add_feature(cfeature.ShapelyFeature(CN_shp_info.geometries(), projection,\
                                               edgecolor='k', facecolor='none'), linewidth=0.5, \
                                                linestyle='-' , zorder=5, alpha=0.8)
        if continuously:
            cbar = fig.colorbar(mappable=pm, cax=cx, ax=ax, orientation=orient)
        else:
            cbar = fig.colorbar(mappable=pm, cax=cx, ax=ax, orientation=orient, ticks=levels)
            cbar.set_ticklabels(RadarGraphMap._FixTicks(levels))
        cbar.set_label(clabel)
        #ax.set_aspect("equal")
        ax.set_title(title, fontsize=15)
        parallels = np.arange(int(min_lat), math.ceil(max_lat) + 1, 1)
        meridians = np.arange(int(min_lon), math.ceil(max_lon) + 1, 1)
        ax.set_xticks(meridians, crs=projection)
        ax.set_yticks(parallels, crs=projection)
        lon_formatter = LongitudeFormatter()
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)

    @staticmethod
    def _FixTicks(ticks):
        """修改ticks的显示label"""
        if (ticks % 1).sum() == 0:
            temp = ["%2.f" % i for i in ticks]
        else:
            temp = ["%.2f" % i for i in ticks]
        return temp




