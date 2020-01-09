# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 22:19:11 2019

this script is created to plot Stereotypical radar  images

@author: zy
"""
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from ..configure.default_config import CINRAD_COLORMAP, CINRAD_field_bins, \
    CINRAD_field_normvar, CINRAD_field_mapping, DEFAULT_METADATA
from ..core.transforms import antenna_vectors_to_cartesian_cwr

class RadarGraph(object):
    """雷达绘图显示部分"""

    def __init__(self, NuistRadar=None):
        self.NRadar = NuistRadar

    def plot(self, sweep, field_name, normvar=None, title=None, clabel=None, dark=False, continuously=False):
        """
        绘图
        :param sweep: sweep从0开始
        :param field_name: 产品场的名称
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
            title = title + " Elevation : %.1f"%self.NRadar.scan_info.fixed_angle[sweep].values
        if clabel is None:
            clabel = CINRAD_field_mapping[field_name] + " (%s)"%DEFAULT_METADATA[CINRAD_field_mapping[field_name]]['units']
        field = self.NRadar.fields[sweep][field_name]
        RadarGraph.simple_plot_ppi(radar_data=field, normvar=(vmin, vmax),
                                   title=title, cmap=cmap, cmap_bins=cmap_bins,
                                   clabel=clabel, dark=dark, continuously=continuously)
    @staticmethod
    def GUI_plot(NRadar, fig, ax, cx, sweep, field_name, normvar=None, title=None, clabel=None, continuously=False):
        """
        绘图
        :param sweep: sweep从0开始
        :param field_name: 产品场的名称
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
            title = title + " Elevation : %.1f"%NRadar.scan_info.fixed_angle[sweep].values
        if clabel is None:
            clabel = CINRAD_field_mapping[field_name]  + " (%s)"%DEFAULT_METADATA[CINRAD_field_mapping[field_name]]['units']
        field = NRadar.fields[sweep][field_name]
        assert hasattr(field, "x") and hasattr(field, "y"), "check Nradar coords!"
        x, y = field.x, field.y
        RadarGraph.plot_ppi(fig, ax, cx, x, y, field, normvar=(vmin, vmax),
                                title=title, cmap=cmap, cmap_bins=cmap_bins,
                                clabel=clabel, continuously=continuously)

    @staticmethod
    def simple_plot_ppi_xy(x, y, radar_data, normvar=None, cmap=None,
                           max_range=None, title=None, cmap_bins=16, orient="vertical",
                           label=None, clabel=None, dark=False, continuously=False):
        """
        根据方位角 距离画图
        :param x: array x轴方向的距离(units:m)
        :param y: array y轴方向的距离(units:m)
        :param radar_data: array 雷达数据的产品
        :param normvar:[vmin, vmax]要画的最大值最小值
        :param cmap: str要使用的cmaps
        :param max_range: 要画的最大半径
        :param title: 图的title
        :param cmap_bins: cmap要分的色段
        :param orient: cmap的方向, 默认是垂直
        :param label: x和y轴的label
        :param clabel: cbar的label
        :param dark: 是否使用黑色背景
        :return:
        """
        if dark:
            plt.style.use('dark_background')
        fig = plt.figure()
        ax = fig.add_axes([0.08, 0.1, 0.82, 0.82])
        cax = fig.add_axes([0.85, 0.1, 0.028, 0.82])
        ax.tick_params(axis="y", which="both", direction='in')
        ax.tick_params(axis="x", which="both", direction='in')
        return RadarGraph.plot_ppi(fig, ax, cax, x, y, radar_data, max_range, title, normvar, cmap, \
                                   cmap_bins, orient, label, clabel, continuously)

    @staticmethod
    def simple_plot_ppi(radar_data, _range=None, azimuth=None, elevation=None, normvar=None, cmap=None,
                        max_range=None, title=None, cmap_bins=16, orient="vertical",
                        label=None, clabel=None, dark=False, continuously=False):
        """
        根据方位角 距离画图
        :param _range: array 每个库的距离(units:m)
        :param azimuth: array 方位角(units:degree)
        :param elevation: array 每个径向的仰角(units:degree)
        :param radar_data: array 雷达数据的产品
        :param normvar:[vmin, vmax]要画的最大值最小值
        :param cmap: str要使用的cmaps
        :param max_range: 要画的最大半径
        :param title: 图的title
        :param cmap_bins: cmap要分的色段
        :param orient: cmap的方向, 默认是垂直
        :param label: x和y轴的label
        :param clabel: cbar的label
        :param dark: 是否使用黑色背景
        :return:
        """
        if hasattr(elevation, "values"):
            elevation = elevation.values
        if hasattr(azimuth, "values"):
            azimuth = azimuth.values
        if hasattr(_range, "values"):
            _range = _range.values
        assert radar_data is not None, "radar_data should not be None!"
        if hasattr(radar_data, "x") and hasattr(radar_data, "y"):
            x, y = radar_data.x, radar_data.y
        else:
            x, y, z = antenna_vectors_to_cartesian_cwr(_range, azimuth, elevation, h=0)
        if dark:
            plt.style.use('dark_background')
        fig = plt.figure()
        ax = fig.add_axes([0.08, 0.1, 0.82, 0.82])
        cax = fig.add_axes([0.85, 0.1, 0.028, 0.82])
        ax.tick_params(axis="y", which="both", direction='in')
        ax.tick_params(axis="x", which="both", direction='in')
        return RadarGraph.plot_ppi(fig, ax, cax, x, y, radar_data, max_range, title, normvar, cmap, \
                                   cmap_bins, orient, label, clabel, continuously)

    @staticmethod
    def plot_ppi(fig, ax, cx, x, y, radar_data, max_range=None, title=None, normvar=None, cmap=None, \
                 cmap_bins=16, orient="vertical", label=None, clabel=None, continuously=False):
        """
        绘制ppi图像
        :param fig matplotlib的fig对象
        :param ax axes的对象
        :param cx colorbar的axes的对象
        :param x: x是指数据对应的x轴坐标, units:m
        :param y: y是指数据对应的y轴坐标, units:m
        :param radar_data: 雷达数据 //dBZ, V, W, KDP, PhiDP, ZDR, RHV
        :param max_range: 最大显示的距离,(min_x, max_x, min_y, max_y), units:m
        :param title: 图像显示的title
        :param normvar: list or tuple, (vmin, vmax), 作图的最大值与最小值
        :param cmap: 选择画图的colorbar
        :param orient:代表colormap的方向
        :param label:label=[xlabel,ylabel]代表x轴y轴要显示的字符串
        :param clabel:代表cbar上的label
        :param dark: 是否打开黑色北京
        :param continuously: bool,采用连续的colorbar,cmap_bins失效
        :return:
        """
        if label is None:
            xlabel = "Distance From Radar In East (Uints:km)"
            ylabel = "Distance From Radar In North (Uints:km)"
        else:
            xlabel, ylable = label
        if normvar is None:
            vmin = -5
            vmax = 75
        else:
            vmin, vmax = normvar
        if max_range is None:
            range_cycle = np.max(x)
            min_x, max_x, min_y, max_y = -range_cycle / 1000., range_cycle / 1000.,\
                                         -range_cycle/1000., range_cycle/1000.
        else:
            min_x, max_x, min_y, max_y = max_range
            range_cycle = max(min_x, max_x, min_y, max_y)
        if title is None:
            title = ""
        if clabel is None:
            clabel = ""
        if continuously:
            cmap_bins = 256
        cmaps = plt.get_cmap(cmap)
        levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmaps.N, clip=True)
        gci = ax.pcolormesh(x / 1000., y / 1000., radar_data, cmap=cmaps, \
                            zorder=0, norm=norm)
        RadarGraph._SetGrids(ax, range_cycle / 1000.)
        ax.set_aspect("equal")
        ax.set_xlim([min_x, max_x])
        ax.set_ylim([min_y, max_y])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if continuously:
            cbar = fig.colorbar(mappable=gci, cax=cx, orientation=orient)
        else:
            cbar = fig.colorbar(mappable=gci, cax=cx, orientation=orient, ticks=levels)
            cbar.set_ticklabels(RadarGraph._FixTicks(levels))
        cbar.set_label(clabel)
        RadarGraph._SetAxis(ax, 50, 5, 8)
        plt.style.use('default')

    @staticmethod
    def _SetAxis(ax, major, minor, fontsize):
        """设置主次刻度, major代表主刻度每隔多少一个,minor代表次刻度每隔多少一个,以及字体大小"""
        ax.xaxis.set_major_locator(MultipleLocator(major))
        ax.xaxis.set_minor_locator(MultipleLocator(minor))
        ax.yaxis.set_major_locator(MultipleLocator(major))
        ax.yaxis.set_minor_locator(MultipleLocator(minor))
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)

    @staticmethod
    def _SetGrids(ax, max_range, deltaR=50):
        """绘制背景格点, ranges代表数据距离中心点的最大值, deltaR代表每隔多少画一个圆"""
        theta = np.linspace(0, 2 * np.pi, 200)
        rings = np.arange(deltaR, max_range + 1, deltaR)
        for i in rings:
            x0 = i * np.cos(theta)
            y0 = i * np.sin(theta)
            ax.plot(x0, y0, linestyle='-', linewidth=0.6, color='#5B5B5B')

        for rad in np.arange(0, np.pi, np.pi / 6):
            ax.plot([-1 * rings[-1] * np.sin(rad), rings[-1] * np.sin(rad)], \
                    [-1 * rings[-1] * np.cos(rad), rings[-1] * np.cos(rad)], \
                    linestyle='-', linewidth=0.6, color='#5B5B5B')

    @staticmethod
    def _FixTicks(ticks):
        """修改ticks的显示label"""
        if (ticks % 1).sum() == 0:
            temp = ["%2.f" % i for i in ticks]
        else:
            temp = ["%.2f" % i for i in ticks]
        return temp
