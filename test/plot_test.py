from pycwr.io import read_PA, read_auto, PAFile
from pycwr.draw.RadarPlot import Graph
import matplotlib.pyplot as plt
import numpy as np

file = "/Users/zhengyu/Downloads/Z_RADR_I_ZGZ01_20200820220246_O_DOR_DXK_CAR.bin.bz2"
PRD = read_PA(file)
x = PAFile.PABaseData(file)


fig, ax = plt.subplots()
graph = Graph(PRD)
graph.plot_ppi(ax, 0, "dBZ")
plt.show()
