import sys
sys.path.append("../")
from libs.io.auto_io import radar_io

radar_obj = radar_io(r"E:\RadarBaseData\StandardFormat\厦门\Z9592.20160728.111443.AR2.bz2")
NRadar = radar_obj.ToNuistRadar()
PyartRadar = radar_obj.ToPyartRadar()