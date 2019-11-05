# -*- coding: utf-8 -*-
from NuistRadar.io.auto_io import radar_io
from NuistRadar.io.SCFile import SCBaseData
from NuistRadar.draw.VerticalSectionPlot import VerticalSection
import numpy as np
import matplotlib.pyplot as plt

#file = r"E:\RadarBaseData\CINRAD-SA\温州\2015080816.59A"
file = r"C:\Users\zy\Desktop\Z_RADR_I_Z9831_20190427200200_O_DOR_SC_CAP.bin"
#radar = SCBaseData(file)
Radar = radar_io(file)
NRadar = Radar.ToNuistRadar()
PyartRadar = Radar.ToPyartRadar()
