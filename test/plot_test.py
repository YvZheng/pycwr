from pycwr.io import read_PA, read_auto
from pycwr.draw.SingleRadarPlotMap import RadarGraphMap
import matplotlib.pyplot as plt

file = "/Users/zhengyu/Downloads/Z_RADR_I_ZGZ01_20200820220246_O_DOR_DXK_CAR.bin.bz2"
x = read_auto(file)

# PRD = read_auto(filename)
# graph = RadarGraphMap(PRD)
# graph.plot(0, "dBZ")
# plt.show()