from pycwr.io import read_PA, read_auto, PAFile
from pycwr.draw.SingleRadarPlotMap import RadarGraphMap
import matplotlib.pyplot as plt

file = "/Users/zhengyu/Downloads/Z_RADR_I_Z0001_20200521191950_O_DOR_DXK_CAR.bin"
x = read_PA(file)
y = PAFile.PABaseData(file)

# PRD = read_auto(filename)
# graph = RadarGraphMap(PRD)
# graph.plot(0, "dBZ")
# plt.show()