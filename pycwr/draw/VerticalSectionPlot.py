# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from ..configure.default_config import CINRAD_COLORMAP, CINRAD_field_bins,\
    CINRAD_field_normvar, CINRAD_field_mapping
import pandas as pd

class VerticalSection(object):

    def __init__(self, NuistRadar):
        self.NRadar = NuistRadar

    def section(self, start_point, end_point, field_name, title=None, clabel=None, continuously=False):
        """
        :param start_point:剖面的起点(x, y)
        :param end_point:剖面的结束点(x, y)
        :param field_name:要剖的数据场
        :param title: 剖面图的title
        :param clabel: colorbar的title
        :param continuously: 是否使用连续的colorbar
        :return:
        """
        self.mesh_xy, self.mesh_z, self.grid_field = self.NRadar.get_vertical_section(start_point,
                                                                            end_point, field_name)
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
        return VerticalSection.SectionPlot(fig, ax, cax, self.mesh_xy, self.mesh_z, self.grid_field, title=title,
                                           normvar=(vmin, vmax), cmap=cmap, cmap_bins=cmap_bins, clabel=clabel,
                                           continuously=continuously)

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
        mesh_xy, mesh_z, grid_field = NRadar.get_vertical_section(start_point, end_point, field_name)
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
        return VerticalSection.SectionPlot(fig, ax, cax, mesh_xy, mesh_z, grid_field, title=title,
                                           normvar=(vmin, vmax), cmap=cmap, cmap_bins=cmap_bins, clabel=clabel,
                                           continuously=continuously)

    @staticmethod
    def SectionPlot(fig, ax, cx, mesh_xy, mesh_z, field_data, title=None, normvar=None, cmap=None, \
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
            ylabel = "Distance Above on Radar (Uints:km)"
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
        cmaps = plt.get_cmap(cmap)
        levels = MaxNLocator(nbins=cmap_bins).tick_values(vmin, vmax)
        norm = BoundaryNorm(levels, ncolors=cmaps.N, clip=True)
        gci = ax.pcolormesh(mesh_xy / 1000., mesh_z / 1000., field_data, cmap=cmaps, \
                            zorder=10, norm=norm)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if continuously:
            cbar = fig.colorbar(mappable=gci, cax=cx, orientation=orient)
        else:
            cbar = fig.colorbar(mappable=gci, cax=cx, orientation=orient, ticks=levels)
            cbar.set_ticklabels(VerticalSection._FixTicks(levels))
        cbar.set_label(clabel)
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
