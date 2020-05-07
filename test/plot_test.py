from pycwr.io import read_auto
from pycwr.draw.SingleRadarPlotMap import RadarGraphMap
import matplotlib.pyplot as plt

filename = r"E:\RadarBaseData\WSR98D\Z9898\Z_RADR_I_Z9898_20190828192401_O_DOR_SAD_CAP_FMT.bin.bz2"

PRD = read_auto(filename)
graph = RadarGraphMap(PRD)
graph.plot(0, "dBZ")
plt.show()