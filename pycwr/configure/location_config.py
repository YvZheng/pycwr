# -*- coding: utf-8 -*-
import os

import pandas as pd

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Radar_info_Path = os.path.join(ROOT_DIR, "data", "radar_info.json")
CN_shp_path = os.path.join(ROOT_DIR, "data", "CHN_pronvices.shp")
last_open_dir = os.path.join(ROOT_DIR, "data", "default_opendir.json")

radar_info = pd.read_json(Radar_info_Path)
_CN_SHP_INFO = None


def get_cn_shp_info():
    """Return the packaged province shapefile reader when cartopy is installed."""
    global _CN_SHP_INFO
    if _CN_SHP_INFO is None:
        try:
            import cartopy.io.shapereader as shpreader
        except ImportError as exc:
            raise ImportError(
                "Map plotting requires cartopy. Install the full optional stack with "
                "`pip install \"pycwr[full]\"`."
            ) from exc
        _CN_SHP_INFO = shpreader.Reader(CN_shp_path)
    return _CN_SHP_INFO


def __getattr__(name):
    if name == "CN_shp_info":
        return get_cn_shp_info()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
