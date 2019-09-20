import struct
import bz2
import gzip
import datetime

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