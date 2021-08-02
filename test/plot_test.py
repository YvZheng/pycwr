from pycwr.io import read_PA, read_auto, read_WSR98D, WSR98DFile
from pycwr.draw.SingleRadarPlotMap import RadarGraphMap
import matplotlib.pyplot as plt

file = "/Users/zhengyu/Downloads/2021051121050000.can"
x = WSR98DFile.WSR98DBaseData(file)

y = WSR98DFile.WSR98D2NRadar(x)

z = read_WSR98D(file)
# PRD = read_auto(filename)
# graph = RadarGraphMap(PRD)
# graph.plot(0, "dBZ")
# plt.show()