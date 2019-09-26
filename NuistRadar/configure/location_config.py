# -*- coding: utf-8 -*-
import cartopy.io.shapereader as shpreader
import pandas as pd
import os

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Radar_info_Path = os.path.join(ROOT_DIR, "data", "radar_info.json")
CN_shp_path = os.path.join(ROOT_DIR, "data", "CHN_pronvices.shp")

radar_info = pd.read_json(Radar_info_Path)
CN_shp_info = shpreader.Reader(CN_shp_path)