# -*- coding: utf-8 -*-
from pycwr.io.auto_io import radar_io
from pycwr.draw.VerticalSectionPlot import VerticalSection
from pycwr.draw.RadarPlot import Graph, GraphMap, plot_az_ranges, plot_az_ranges_map
from pycwr.draw.SingleRadarPlotMap import RadarGraphMap
from pycwr.draw.SingleRadarPlot import RadarGraph
from pycwr.qc.attenuation import correct_attenuation_HB, correct_attenuation
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

#file = r"E:\RadarBaseData\CINRAD-SA\温州\2015080816.59A"
#file = r"E:\RadarBaseData\CINRAD-CC\2016070817.03V"
file = r"C:\Users\zy\Desktop\Z9577.20190809.100149.AR2.bz2"
#file = r"E:\RadarBaseData\WSR98D\hangzhou\Z9040.20190515.040618.AR2.bz2"
#radar = SCBaseData(file)
Radar = radar_io(file)
NRadar = Radar.ToPRD()


fig, ax = plt.subplots()
graph = Graph(NRadar)
graph.plot_ppi(ax, 0, "dBZ", cmap="pyart_NWSRef")
graph.add_rings(ax, [0, 50, 100, 150, 200, 250, 300])
ax.set_title("example of PPI", fontsize=16)
ax.set_xlabel("Distance From Radar In East (km)", fontsize=14)
ax.set_ylabel("Distance From Radar In North (km)", fontsize=14)
#ax.set_ylim([0,15])
#ax.set_ylabel("Height (km)", fontsize=14)
#ax.set_xlabel("Latitude, Longitude", fontsize=14)
#ax.set_title("VCS exmaple", fontsize=16)
#ax.set_title("example of PPI with map", fontsize=16)

#ax2 = plt.axes(projection=ccrs.PlateCarree())
#graph.plot_ppi_map(ax2, 0, "dBZ", cmap="pyart_NWSRef")
#graph.add_lines_map(ax2, (120.8, 27.8), (122.9, 26.8))
#plt.show()

#vcs.section((0,0), (150, 0), "dBZ", (0, 10))
#ax.set_extent([117, 124, 25, 31], ccrs.PlateCarree())
#graph.add_lines_map(ax, (120.8, 27.8), (122.9, 26.8))
#ax.plot([120.8, 122.9], [27.8, 26.8], color="k", marker="x", transform=ccrs.PlateCarree(),zorder=20)
#gci2 = graph_obj.plot_ppi(ax2, 0, "V", cmap="pyart_NWSVel", cmap_bins=256)
#graph_obj.plot_vcs(ax1, (0,0), (150,0), "dBZ")
#ax1.set_ylim([0,10])
#graph_obj.add_rings(ax1, [0, 50, 100, 150])
#graph_obj.add_lines(ax, [0, 150], [0, 0], linewidth=2, color="k", markersize=12)
#vcs.section_map((120.8, 27.8), (122.9, 26.8), "dBZ", (0, 10))
#dbz = NRadar.fields[0]["dBZ"]
#y = correct_attenuation_HB(dbz)
#graph = RadarGraph(NRadar)
#z = correct_attenuation(dbz, "X")
#graph.plot(0, "dBZ", continuously=False)
#graph.simple_plot_ppi_map(dbz.range, dbz.azimuth, dbz.elevation, dbz,(0,0), cmap="CN_ref")
plt.savefig("../examples/graph.png")
plt.show()
#ax = plt.axes()
#point = (NRadar.scan_info.longitude.values, NRadar.scan_info.latitude.values)
#plot_az_ranges_map(ax, NRadar.fields[0].range.values, NRadar.fields[0].azimuth.values,\
#               NRadar.fields[0].elevation.values, NRadar.fields[0].dBZ, main_point=point ,\
#                   transform=ccrs.PlateCarree(), min_max=(-5,75))
