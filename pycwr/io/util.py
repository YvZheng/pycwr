import sys
sys.path.append("../../")
import struct
import bz2
import gzip
import datetime
import os
from ..configure.location_config import radar_info

def _structure_size(structure):
    """计算structure的字节大小"""
    return struct.calcsize('<' + ''.join([i[1] for i in structure]))

def _unpack_from_buf(buf, pos, structure):
    '''unpack a structure from buf'''
    size = _structure_size(structure)
    return _unpack_structure(buf[pos:pos + size], structure), size

def _unpack_structure(string, structure):
    '''unpack a structure from a string'''
    '''unpack a structure from a string'''
    fmt = '<' + ''.join([i[1] for i in structure])
    lst = struct.unpack(fmt, string)
    return dict(zip([i[0] for i in structure], lst))

def _prepare_for_read(filename):
    """
    Return a file like object read for reading.
    Open a file for reading in binary mode with transparent decompression of
    Gzip and BZip2 files.  The resulting file-like object should be closed.
    Parameters
    ----------
    filename : str or file-like object
        Filename or file-like object which will be opened.  File-like objects
        will not be examined for compressed data.
    Returns
    -------
    file_like : file-like object
        File like object from which data can be read.
    """
    # if a file-like object was provided, return
    if hasattr(filename, 'read'):  # file-like object
        return filename
    # look for compressed data by examining the first few bytes
    fh = open(filename, 'rb')
    magic = fh.read(3)
    fh.close()
    if magic.startswith(b'\x1f\x8b'):
        f = gzip.GzipFile(filename, 'rb')
    elif magic.startswith(b'BZh'):
        f = bz2.BZ2File(filename, 'rb')
    else:
        f = open(filename, 'rb')
    return f

def julian2date(JulianDate, Msec):
    """
    faster than num2date in netcdf4
    :param JulianDate: Julian Date
    :param Msec: msec from 00:00
    :return:
    """
    deltday = datetime.timedelta(days=JulianDate)
    deltsec = datetime.timedelta(milliseconds=Msec)
    scantime = datetime.datetime(1969, 12, 31) + deltday + deltsec
    return scantime

def julian2date_SEC(Sec, Msec):
    """
    faster than num2date in netcdf4
    :param Sec: seconds
    :param Msec: microseconds
    :return:
    """
    deltSec = datetime.timedelta(seconds=Sec)
    deltMSec = datetime.timedelta(microseconds=Msec)
    scantime = datetime.datetime(1970, 1, 1) + deltSec + deltMSec
    return scantime

def get_radar_info(filename):
    """
    根据雷达名称找雷达的经纬度信息
    :param filename:
    :return:(lat(deg), lon(deg), elev(m), frequency(GHZ))
    """
    name = os.path.basename(filename)
    try:
        station_id = [int(name[idx:idx+4]) for idx in range(len(name)-4) if name[idx:idx+4].isdigit()][0]
    except Exception:
        station_id = 9250
    if station_id not in radar_info.index:
        station_id = 9250 ###找不到站点信息返回南京雷达
    return radar_info.loc[station_id, "Latitude"], radar_info.loc[station_id, "Longitude"],\
           radar_info.loc[station_id, "Elevation"], radar_info.loc[station_id, "Frequency"]

def get_radar_sitename(filename):
    name = os.path.basename(filename)
    station_id = [int(name[idx:idx + 4]) for idx in range(len(name) - 4) if name[idx:idx + 4].isdigit()][0]
    if station_id not in radar_info.index:
        station_id = 9250  ###找不到站点信息返回南京雷达
    return radar_info.loc[station_id, "Name"]

def _get_radar_type(filename):
    """
    根据雷达名称找雷达类型
    :param filename:
    :return:
    """
    name = os.path.basename(filename)
    station_id = [int(name[idx:idx + 4]) for idx in range(len(name) - 4) if name[idx:idx + 4].isdigit()][0]
    if station_id not in radar_info.index:
        return None
    Datatype = radar_info.loc[station_id, "Datatype"]
    if Datatype in ["SA", "SB", "CB", "SC", "CD"]:
        return "SAB"
    elif Datatype in ["CC", "CCJ"]:
        return "CC"
    else:
        return None

def radar_format(filename):
    if hasattr(filename, 'read'):
        return filename
    fh = _prepare_for_read(filename)
    flag = fh.read(28)
    size = len(fh.read()) + 28
    fh.seek(100, 0)
    sc_flag = fh.read(9)
    fh.seek(116, 0)
    cc_flag = fh.read(9)
    fh.close()
    if flag[:4] == b'RSTM':
        return "WSR98D"
    elif flag[14:16] == b'\x01\x00':
        return "SAB"
    elif (size-1024)%3000 == 0 and cc_flag == b"CINRAD/CC":
        return "CC"
    elif (size-1024)%4000 == 0 and (sc_flag == b"CINRAD/SC" or sc_flag == b"CINRAD/CD"):
        return "SC"
    elif flag[8:12] == b'\x10\x00\x00\x00':
        return "PA"
    else:
        return _get_radar_type(filename)

def make_time_unit_str(dtobj):
    """ Return a time unit string from a datetime object. """
    return "seconds since " + dtobj.strftime("%Y-%m-%dT%H:%M:%SZ")

