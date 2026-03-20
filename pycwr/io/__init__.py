from importlib import import_module
from .util import radar_format

_READER_MODULES = {"CCFile", "SCFile", "WSR98DFile", "SABFile", "PAFile", "NEXRADLevel2File"}

__all__ = [
    "read_auto",
    "read_CC",
    "read_SC",
    "read_WSR98D",
    "read_SAB",
    "read_PA",
    "write_wsr98d",
    "write_nexrad_level2_msg31",
    "write_nexrad_level2_msg1",
    *sorted(_READER_MODULES),
]

def __getattr__(name):
    if name in _READER_MODULES:
        module = import_module(f".{name}", __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def read_auto(filename, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
    """
    :param filename:  radar basedata filename
    :param station_lon:  radar station longitude //units: degree east
    :param station_lat:  radar station latitude //units:degree north
    :param station_alt:  radar station altitude //units: meters
    """
    radar_type = radar_format(filename)
    if radar_type == "WSR98D":
        WSR98DFile = __getattr__("WSR98DFile")
        return WSR98DFile.WSR98D2NRadar(
            WSR98DFile.WSR98DBaseData(filename, station_lon, station_lat, station_alt)
        ).ToPRD(effective_earth_radius=effective_earth_radius)
    elif radar_type == "SAB":
        SABFile = __getattr__("SABFile")
        return SABFile.SAB2NRadar(
            SABFile.SABBaseData(filename, station_lon, station_lat, station_alt)
        ).ToPRD(effective_earth_radius=effective_earth_radius)
    elif radar_type == "CC":
        CCFile = __getattr__("CCFile")
        return CCFile.CC2NRadar(
            CCFile.CCBaseData(filename, station_lon, station_lat, station_alt)
        ).ToPRD(effective_earth_radius=effective_earth_radius)
    elif radar_type == "SC":
        SCFile = __getattr__("SCFile")
        return SCFile.SC2NRadar(
            SCFile.SCBaseData(filename, station_lon, station_lat, station_alt)
        ).ToPRD(effective_earth_radius=effective_earth_radius)
    elif radar_type == "PA":
        PAFile = __getattr__("PAFile")
        return PAFile.PA2NRadar(
            PAFile.PABaseData(filename, station_lon, station_lat, station_alt)
        ).ToPRD(effective_earth_radius=effective_earth_radius)
    else:
        raise TypeError("unsupported radar type!")

def read_SAB(filename, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
    """
    :param filename:  radar basedata filename
    :param station_lon:  radar station longitude //units: degree east
    :param station_lat:  radar station latitude //units:degree north
    :param station_alt:  radar station altitude //units: meters
    """
    SABFile = __getattr__("SABFile")
    return SABFile.SAB2NRadar(
        SABFile.SABBaseData(filename, station_lon, station_lat, station_alt)
    ).ToPRD(effective_earth_radius=effective_earth_radius)

def read_CC(filename, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
    """
    :param filename:  radar basedata filename
    :param station_lon:  radar station longitude //units: degree east
    :param station_lat:  radar station latitude //units:degree north
    :param station_alt:  radar station altitude //units: meters
    """
    CCFile = __getattr__("CCFile")
    return CCFile.CC2NRadar(
        CCFile.CCBaseData(filename, station_lon, station_lat, station_alt)
    ).ToPRD(effective_earth_radius=effective_earth_radius)

def read_SC(filename, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
    """
    :param filename:  radar basedata filename
    :param station_lon:  radar station longitude //units: degree east
    :param station_lat:  radar station latitude //units:degree north
    :param station_alt:  radar station altitude //units: meters
    """
    SCFile = __getattr__("SCFile")
    return SCFile.SC2NRadar(
        SCFile.SCBaseData(filename, station_lon, station_lat, station_alt)
    ).ToPRD(effective_earth_radius=effective_earth_radius)

def read_WSR98D(filename, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
    """
    :param filename:  radar basedata filename
    :param station_lon:  radar station longitude //units: degree east
    :param station_lat:  radar station latitude //units:degree north
    :param station_alt:  radar station altitude //units: meters
    """
    WSR98DFile = __getattr__("WSR98DFile")
    return WSR98DFile.WSR98D2NRadar(
        WSR98DFile.WSR98DBaseData(filename, station_lon, station_lat, station_alt)
    ).ToPRD(effective_earth_radius=effective_earth_radius)

def read_PA(filename, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
    """
    :param filename:  radar basedata filename
    :param station_lon:  radar station longitude //units: degree east
    :param station_lat:  radar station latitude //units:degree north
    :param station_alt:  radar station altitude //units: meters
    """
    PAFile = __getattr__("PAFile")
    return PAFile.PA2NRadar(
        PAFile.PABaseData(filename, station_lon, station_lat, station_alt)
    ).ToPRD(effective_earth_radius=effective_earth_radius)


def write_wsr98d(prd, filename, **kwargs):
    WSR98DFile = __getattr__("WSR98DFile")
    return WSR98DFile.write_wsr98d(prd, filename, **kwargs)


def write_nexrad_level2_msg31(prd, filename, **kwargs):
    NEXRADLevel2File = __getattr__("NEXRADLevel2File")
    return NEXRADLevel2File.write_nexrad_level2_msg31(prd, filename, **kwargs)


def write_nexrad_level2_msg1(prd, filename, **kwargs):
    NEXRADLevel2File = __getattr__("NEXRADLevel2File")
    return NEXRADLevel2File.write_nexrad_level2_msg1(prd, filename, **kwargs)
