# -*- coding: utf-8 -*-
from .util import radar_format
from . import CCFile, SABFile, SCFile, WSR98DFile

class radar_io(object):
    def __init__(self, filename, station_lon=None, station_lat=None, station_alt=None):
        """
        :param filename:  radar basedata filename
        :param station_lon:  radar station longitude //units: degree east
        :param station_lat:  radar station latitude //units:degree north
        :param station_alt:  radar station altitude //units: meters
        """
        radar_type = radar_format(filename)
        if radar_type == "WSR98D":
            self.radar_obj = WSR98DFile.WSR98D2NRadar(WSR98DFile.WSR98DBaseData(filename, station_lon, station_lat, station_alt))
        elif radar_type == "SAB":
            self.radar_obj = SABFile.SAB2NRadar(SABFile.SABBaseData(filename, station_lon, station_lat, station_alt))
        elif radar_type == "CC":
            self.radar_obj = CCFile.CC2NRadar(CCFile.CCBaseData(filename, station_lon, station_lat, station_alt))
        elif radar_type == "SC":
            self.radar_obj = SCFile.SC2NRadar(SCFile.SCBaseData(filename, station_lon, station_lat, station_alt))
        else:
            raise TypeError("unsupported radar type!")

    def ToPRD(self, withlatlon=True):
        return self.radar_obj.ToPRD(withlatlon=withlatlon)

    def ToPyartRadar(self):
        return self.radar_obj.ToPyartRadar()


if __name__ == "__main__":
    import glob
    files = glob.glob("E:\RadarBaseData\**\*.bz2", recursive=True)
    for ifile in files:
        try:
            res = radar_io(ifile)
            c = res.ToPRD()
            d = res.ToPyartRadar()
        except Exception:
            print(ifile)
