# -*- coding: utf-8 -*-
from NuistRadar.io.auto_io import radar_io
from NuistRadar.draw.VerticalSectionPlot import VerticalSection
import numpy as np
import matplotlib.pyplot as plt

file = r"E:\RadarBaseData\CINRAD-SA\温州\2015080816.59A"
NRadar = radar_io(file).ToNuistRadar()

sec = VerticalSection(NRadar)
dat = sec.section((0, 0), (-80000, 80000), "dBZ")
plt.show()

