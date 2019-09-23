# -*- coding: utf-8 -*-
import sys
sys.path.append("../")
import pandas as pd
from cfgs import Radar_info_Path

radar_info = pd.read_json(Radar_info_Path)