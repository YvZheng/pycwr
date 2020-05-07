from pycwr.io import read_auto
import matplotlib.pyplot as plt
import numpy as np
from pycwr.draw.RadarPlot import Graph

file = r"C:\Users\zy\Desktop\HID\NUIST.20160707.001054.AR2"

NRadar = read_auto(file)
num = 3
NRadar.fields[num]['dBZ'][:] = np.where(NRadar.fields[num].CC>0.9, NRadar.fields[num].dBZ, np.nan)
NRadar.fields[num]['KDP'][:] = np.where(NRadar.fields[num].CC>0.9, NRadar.fields[num].KDP, np.nan)
NRadar.fields[num]['ZDR'][:] = np.where(NRadar.fields[num].CC>0.9, NRadar.fields[num].ZDR, np.nan)
NRadar.fields[num]['CC'][:] = np.where(NRadar.fields[num].CC>0.9, NRadar.fields[num].CC, np.nan)

fig, ax = plt.subplots()
graph = Graph(NRadar)
str_="CC"
graph.plot_ppi(ax, 3, str_, min_max=[0.85,1])
graph.add_rings(ax, [0, 50, 100, 150])
ax.set_title("(d) Cross correlation ratio, El : 3.4", fontsize=14, loc="left")
ax.set_xlabel("Distance From Radar In East (km)", fontsize=12)
ax.set_ylabel("Distance From Radar In North (km)", fontsize=12)
ax.set_xlim([-150, 150])
ax.set_ylim([-150, 150])
plt.savefig(r"C:\Users\zy\Desktop\HID\201607070010\%s.png"%str_, dpi=600)
plt.show()