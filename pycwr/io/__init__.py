from . import SCFile, WSR98DFile, SABFile, CCFile
from .util import radar_format

__all__ = ["read_auto", "CCFile", "SCFile", "WSR98DFile", "SABFile"]

def read_auto(filename, station_lon=None, station_lat=None, station_alt=None):
    """
    :param filename:  radar basedata filename
    :param station_lon:  radar station longitude //units: degree east
    :param station_lat:  radar station latitude //units:degree north
    :param station_alt:  radar station altitude //units: meters
    """
    radar_type = radar_format(filename)
    if radar_type == "WSR98D":
        return WSR98DFile.WSR98D2NRadar(WSR98DFile.WSR98DBaseData(filename, station_lon, station_lat, station_alt)).ToPRD()
    elif radar_type == "SAB":
        return SABFile.SAB2NRadar(SABFile.SABBaseData(filename, station_lon, station_lat, station_alt)).ToPRD()
    elif radar_type == "CC":
        return CCFile.CC2NRadar(CCFile.CCBaseData(filename, station_lon, station_lat, station_alt)).ToPRD()
    elif radar_type == "SC":
        return SCFile.SC2NRadar(SCFile.SCBaseData(filename, station_lon, station_lat, station_alt)).ToPRD()
    else:
        raise TypeError("unsupported radar type!")
