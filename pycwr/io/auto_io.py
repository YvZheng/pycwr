# -*- coding: utf-8 -*-
from .util import radar_format
from . import CCFile, SABFile, SCFile, WSR98DFile

class radar_io(object):
    def __init__(self, filename):
        radar_type = radar_format(filename)
        if radar_type == "WSR98D":
            self.radar_obj = WSR98DFile.WSR98D2NRadar(WSR98DFile.WSR98DBaseData(filename))
        elif radar_type == "SAB":
            self.radar_obj = SABFile.SAB2NRadar(SABFile.SABBaseData(filename))
        elif radar_type == "CC":
            self.radar_obj = CCFile.CC2NRadar(CCFile.CCBaseData(filename))
        elif radar_type == "SC":
            self.radar_obj = SCFile.SC2NRadar(SCFile.SCBaseData(filename))
        else:
            raise TypeError("unsupported radar type!")

    def ToNuistRadar(self):
        return self.radar_obj.ToNuistRadar()

    def ToPyartRadar(self):
        return self.radar_obj.ToPyartRadar()


if __name__ == "__main__":
    import glob
    files = glob.glob("E:\RadarBaseData\**\*.bz2", recursive=True)
    for ifile in files:
        try:
            res = radar_io(ifile)
            c = res.ToNuistRadar()
            d = res.ToPyartRadar()
        except Exception:
            print(ifile)
