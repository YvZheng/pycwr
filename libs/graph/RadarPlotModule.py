# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 22:04:11 2018

this script is created to plot Stereotypical radar  images

@author: zy
"""

#import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as col
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from colormap import *
import matplotlib.pyplot as plt

class graph(object):
    '''plot radar graph'''                            
            
    @staticmethod
    def plot_ppi(fig, ax, cx, x, y, dat, label, title, normvar ,orient ,cmap,
                 ticks, clabel, nws = False):
        """绘制ppi图像,x是指数据对应的x轴坐标,单位:千米,y是指数据对应y轴坐标,单位:km, dat是指数据
        label=[xlabel,ylabel]代表x轴y轴要显示的字符串, title是指图像的标题, 
        normvar=[vmin,vmax]画图时要画的范围,orient代表colormap的方向,clabel代表cbar上的label"""
#        global ax
#        plt.style.use('dark_background')
#        plt.style.use('default')
        
        xlabel, ylabel = label
        vmin, vmax = normvar
        ranges = graph._MaxR(x.max())
        xmin = ymin = -1*ranges
        xmax = ymax = ranges
        if nws:
            dat = np.ma.masked_where((dat == -999.) | (dat == 999.), dat)
        else:
            dat = np.ma.masked_where((dat == -999.) | (dat == 999.), dat)
                
        graph._SetGrids(ax,ranges)
        gci = ax.pcolormesh(x, y, dat, cmap=cmap, zorder=100,norm=col.Normalize(vmin,vmax))
        ax.set_aspect("equal")
        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        #cax = make_axes_locatable(ax).append_axes("right", size="5%", pad="5%")
        cbar = fig.colorbar(mappable=gci,cax=cx,orientation=orient , ticks=ticks) 
        #cbar = colorbar(mappable=gci,cax=cx,orientation=orient, ticks=ticks,)
        cbar.set_label(clabel)
        cbar.set_ticklabels(graph._FixTicks(ticks,nws))
        graph._SetAxis(ax,50,5,8)
#        plt.style.use('default')
        
    @staticmethod    
    def _SetAxis(ax,major,minor,fontsize):
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
    def _SetGrids(ax, ranges, deltaR = 50):
        """绘制背景格点, ranges代表数据距离中心点的最大值, deltaR代表每隔多少画一个圆"""        
        theta = np.linspace(0, 2 * np.pi, 200)
        for i in np.arange(deltaR, ranges+1 , deltaR):
            x0 = i*np.cos(theta)
            y0 = i*np.sin(theta)
            ax.plot(x0,y0,linestyle = '-',linewidth=0.6,color='#6F6F6F')
                    
        for rad in np.arange(0,np.pi,np.pi/6):
            ax.plot([-1*ranges*np.sin(rad),ranges*np.sin(rad)],\
                     [-1*ranges*np.cos(rad),ranges*np.cos(rad)],\
                     linestyle = '-',linewidth=0.6,color='#6F6F6F')
    @staticmethod
    def _MaxR(ranges,deltaR = 50):
        """计算图像应该显示的范围"""
        return ranges-ranges%deltaR + deltaR
    
    @staticmethod
    def _FixTicks(ticks,nws):
        """修改ticks的显示label，让最后一个值显示为RF,nws不存在RF"""
        if (ticks%1).sum() == 0:
            temp = ["%2.f"%i for i in ticks]
        else:    
            temp = ["%.1f"%i for i in ticks]
        #if nws == False:
        #    temp[-1] = "FV"
        return temp
    
    @staticmethod
    def plot_ppi_map(fig, ax, cx, lon, lat, dat, title, normvar ,clocation ,cmap, ticks,
                     clabel, extend, main_point, projection, resolution, alpha,
                     parallels, meridians, shapefile = None):
        """绘制带有地图底图的PPI图像, lon数据的经度坐标,lat数据的纬度坐标,title图像显示的标题,
        normvar=[vmin,vmax]画图时要画的范围控制pcolor,clocation是cbar的位置
        (bottom,left,right,up)可选,cmap是画图所选取的cmap,ticks是cbar上显示的刻度,clabel
        是cbar上的label, extend=(min_lon, min_lat, max_lon, max_lat)确定图像显示的范围
        main_point = (lon_0, lat_0)雷达站点的位置 ,projection是basemap的投影方式,lcc
        是basemap图像的分辨率(l,c,h),parallels是绘制的纬度网格,meridians绘制的经度网格,
        均为一维矩阵,shapefile保存的省界数据"""
        
        min_lon, min_lat, max_lon, max_lat = extend
        lon_0, lat_0 = main_point
        vmin, vmax = normvar
        
        dat = np.ma.masked_where((dat == -999) | (dat == 999), dat)
        
        m = Basemap(llcrnrlon=min_lon, llcrnrlat=min_lat,urcrnrlon=max_lon,
                    urcrnrlat=max_lat,lat_0=lat_0, lon_0=lon_0, projection=\
                    projection,resolution=resolution, ax=ax,)
               
        m.fillcontinents(color = "#D3D3D3",ax = ax, zorder=0)  
        if shapefile is not None:
            m.readshapefile(shapefile, 'states', drawbounds=True,
                            zorder=1, ax = ax)
        else:
            m.drawcoastlines()
            m.drawcountries()
            m.drawrivers()
            
        pm = m.pcolormesh(lon, lat, dat, cmap = cmap, norm=col.Normalize(vmin,vmax),
                          alpha = alpha,latlon = True,zorder=2)
        
        m.plot(lon_0, lat_0, marker = 'o', color = '#FFFFFF', latlon = True,
               markersize = 3)
        #m.plot(119.9,32.516,latlon=True,marker="*",color="r",markersize=10)
        #cax = make_axes_locatable(ax).append_axes("right", size="5%", pad="5%")
        #cbar = colorbar(mappable=pm, cax=cx, orientation="vertical", ticks=ticks,)
        cbar = fig.colorbar(mappable=pm,cax=cx,orientation="vertical" , ticks=ticks)        
        cbar.set_label(clabel)
        cbar.set_ticklabels(graph._FixTicks(ticks,True))
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=16)

        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10,color="#696969",linewidth = 0.4) # 绘制纬线
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10,color="#696969",linewidth = 0.4) # 绘制经线
        #plt.savefig("taizhou.pdf",format="pdf",ppi=350)