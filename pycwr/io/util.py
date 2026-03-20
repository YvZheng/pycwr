import struct
import bz2
import gzip
import datetime
import os
import numpy as np
from ..configure.location_config import radar_info


DEFAULT_MAX_READ_BYTES = int(os.environ.get("PYCWR_MAX_READ_BYTES", str(256 * 1024 * 1024)))
DEFAULT_CHUNK_SIZE = 1024 * 1024


def _validate_max_read_bytes(max_bytes=None):
    limit = DEFAULT_MAX_READ_BYTES if max_bytes is None else int(max_bytes)
    if limit <= 0:
        raise ValueError("Maximum read size must be positive.")
    return limit


def _raise_truncated(context, expected=None, actual=None):
    if expected is None or actual is None:
        raise ValueError("%s is truncated or malformed." % context)
    raise ValueError("%s is truncated or malformed: expected %s bytes, got %s." % (context, expected, actual))

def _structure_size(structure):
    """Return the byte size of a parsed binary structure."""
    return struct.calcsize('<' + ''.join([i[1] for i in structure]))

def _unpack_from_buf(buf, pos, structure):
    '''unpack a structure from buf'''
    size = _structure_size(structure)
    if pos < 0:
        raise ValueError("Buffer offset must be non-negative.")
    if len(buf) < pos + size:
        _raise_truncated("Binary structure", expected=size, actual=max(len(buf) - pos, 0))
    return _unpack_structure(buf[pos:pos + size], structure), size

def _unpack_structure(string, structure):
    '''unpack a structure from a string'''
    '''unpack a structure from a string'''
    fmt = '<' + ''.join([i[1] for i in structure])
    expected = struct.calcsize(fmt)
    if len(string) != expected:
        _raise_truncated("Binary structure", expected=expected, actual=len(string))
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


def _read_exact(file_obj, size, context):
    size = int(size)
    if size < 0:
        raise ValueError("%s requested a negative read size." % context)
    data = file_obj.read(size)
    if len(data) != size:
        _raise_truncated(context, expected=size, actual=len(data))
    return data


def _read_all(file_obj, context, max_bytes=None):
    max_bytes = _validate_max_read_bytes(max_bytes)
    chunks = []
    total = 0
    while True:
        chunk = file_obj.read(DEFAULT_CHUNK_SIZE)
        if not chunk:
            break
        total += len(chunk)
        if total > max_bytes:
            raise ValueError(
                "%s exceeds the maximum allowed decoded size of %s bytes." % (context, max_bytes)
            )
        chunks.append(chunk)
    return b"".join(chunks)


def _validate_count(name, value, minimum=0, maximum=None):
    value = int(value)
    if value < minimum:
        raise ValueError("%s must be >= %s." % (name, minimum))
    if maximum is not None and value > maximum:
        raise ValueError("%s exceeds the maximum allowed value of %s." % (name, maximum))
    return value


def _validate_offset_length(buffer_size, offset, length, context):
    offset = int(offset)
    length = int(length)
    if offset < 0 or length < 0:
        raise ValueError("%s uses a negative offset or length." % context)
    if offset > buffer_size or length > (buffer_size - offset):
        raise ValueError("%s points outside the available buffer." % context)
    return offset, length

def julian2date(JulianDate, Msec):
    """
    faster than num2date in netcdf4
    :param JulianDate: Julian Date
    :param Msec: msec from 00:00
    :return:
    """
    deltday = datetime.timedelta(days=int(JulianDate))
    deltsec = datetime.timedelta(milliseconds=int(Msec))
    scantime = datetime.datetime(1969, 12, 31) + deltday + deltsec
    return scantime

def julian2date_SEC(Sec, Msec):
    """
    faster than num2date in netcdf4
    :param Sec: seconds
    :param Msec: microseconds
    :return:
    """
    deltSec = datetime.timedelta(seconds=int(Sec))
    deltMSec = datetime.timedelta(microseconds=int(Msec))
    scantime = datetime.datetime(1970, 1, 1) + deltSec + deltMSec
    return scantime

def get_radar_info(filename):
    """
    Look up radar site metadata from the station identifier embedded in the filename.
    :param filename:
    :return:(lat(deg), lon(deg), elev(m), frequency(GHZ))
    """
    name = os.path.basename(filename)
    try:
        station_id = [int(name[idx:idx+4]) for idx in range(len(name)-4) if name[idx:idx+4].isdigit()][0]
    except Exception:
        station_id = 9250
    if station_id not in radar_info.index:
        station_id = 9250
    return radar_info.loc[station_id, "Latitude"], radar_info.loc[station_id, "Longitude"],\
           radar_info.loc[station_id, "Elevation"], radar_info.loc[station_id, "Frequency"]

def get_radar_sitename(filename):
    name = os.path.basename(filename)
    station_ids = [int(name[idx:idx + 4]) for idx in range(len(name) - 4) if name[idx:idx + 4].isdigit()]
    if not station_ids:
        return "Unknown"
    station_id = station_ids[0]
    if station_id not in radar_info.index:
        station_id = 9250
    return radar_info.loc[station_id, "Name"]

def _get_radar_type(filename):
    """
    Infer the radar family from the station identifier embedded in the filename.
    :param filename:
    :return:
    """
    name = os.path.basename(filename)
    station_ids = [int(name[idx:idx + 4]) for idx in range(len(name) - 4) if name[idx:idx + 4].isdigit()]
    if not station_ids:
        return None
    station_id = station_ids[0]
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
    try:
        flag = fh.read(28)
        if flag[:4] == b'RSTM':
            return "WSR98D"
        if flag[14:16] == b'\x01\x00':
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            if size >= 2432 and ((size % 2432 == 0) or (size % 4132 == 0) or (size % 3132 == 0)):
                return "SAB"
            fh.seek(28, 0)
        if flag[8:12] == b'\x10\x00\x00\x00':
            return "PA"

        fh.seek(100, 0)
        sc_flag = fh.read(9)
        fh.seek(116, 0)
        cc_flag = fh.read(9)

        if cc_flag == b"CINRAD/CC" or sc_flag in (b"CINRAD/SC", b"CINRAD/CD"):
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            if (size - 1024) % 3000 == 0 and cc_flag == b"CINRAD/CC":
                return "CC"
            if (size - 1024) % 4000 == 0 and sc_flag in (b"CINRAD/SC", b"CINRAD/CD"):
                return "SC"
        return _get_radar_type(filename)
    finally:
        fh.close()

def make_time_unit_str(dtobj):
    """ Return a time unit string from a datetime object. """
    return "seconds since " + dtobj.strftime("%Y-%m-%dT%H:%M:%SZ")

def date2num(dates, units):
    """Convert datetimes to numeric seconds for the unit format used by pycwr."""
    prefix = "seconds since "
    assert units.startswith(prefix), "only second-based time units are supported"
    base_time = datetime.datetime.strptime(units[len(prefix):], "%Y-%m-%dT%H:%M:%SZ")
    if np.isscalar(dates):
        return (dates - base_time).total_seconds()
    return np.array([(date - base_time).total_seconds() for date in dates], dtype=np.float64)
