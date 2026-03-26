# -*- coding: utf-8 -*-
import datetime
import os
import struct
from collections import OrderedDict

import numpy as np
from .BaseDataProtocol.WSR98DProtocol import dtype_98D
from .util import (
    _prepare_for_read,
    _read_all,
    _read_exact,
    _unpack_from_buf,
    _validate_count,
    julian2date_SEC,
    make_time_unit_str,
    date2num,
)
from ..core.NRadar import PRD
from ..configure.pyart_config import get_metadata, get_fillvalue
from ..configure.default_config import CINRAD_field_mapping
from ..core.PyartRadar import Radar


MAX_WSR98D_SWEEPS = 512
MAX_WSR98D_MOMENTS_PER_RADIAL = 64
MAX_WSR98D_MOMENT_DATA_BYTES = 32 * 1024 * 1024
WSR98D_MAX_DECODED_PAYLOAD_BYTES = 1024 * 1024 * 1024

WSR98D_WRITE_FIELD_SPECS = OrderedDict(
    [
        ("dBT", {"data_type": 1, "scale": 2, "offset": 66, "bin_length": 1}),
        ("dBZ", {"data_type": 2, "scale": 2, "offset": 66, "bin_length": 1}),
        ("V", {"data_type": 3, "scale": 2, "offset": 129, "bin_length": 1}),
        ("W", {"data_type": 4, "scale": 2, "offset": 129, "bin_length": 1}),
        ("SQI", {"data_type": 5, "scale": 200, "offset": 0, "bin_length": 1}),
        ("CPA", {"data_type": 6, "scale": 2, "offset": 129, "bin_length": 1}),
        ("ZDR", {"data_type": 7, "scale": 16, "offset": 130, "bin_length": 1}),
        ("LDR", {"data_type": 8, "scale": 16, "offset": 130, "bin_length": 1}),
        ("CC", {"data_type": 9, "scale": 200, "offset": 5, "bin_length": 1}),
        ("PhiDP", {"data_type": 10, "scale": 100, "offset": 50, "bin_length": 2}),
        ("KDP", {"data_type": 11, "scale": 10, "offset": 50, "bin_length": 1}),
        ("CP", {"data_type": 12, "scale": 2, "offset": 129, "bin_length": 1}),
        ("HCL", {"data_type": 14, "scale": 1, "offset": 0, "bin_length": 1}),
        ("CF", {"data_type": 15, "scale": 2, "offset": 129, "bin_length": 1}),
        ("SNRH", {"data_type": 16, "scale": 2, "offset": 20, "bin_length": 1}),
        ("SNRV", {"data_type": 17, "scale": 2, "offset": 20, "bin_length": 1}),
        ("Zc", {"data_type": 32, "scale": 2, "offset": 66, "bin_length": 1}),
        ("Vc", {"data_type": 33, "scale": 2, "offset": 129, "bin_length": 1}),
        ("Wc", {"data_type": 34, "scale": 2, "offset": 129, "bin_length": 1}),
        ("ZDRc", {"data_type": 35, "scale": 16, "offset": 130, "bin_length": 1}),
    ]
)

_WSR98D_RADIAL_HEADER = dtype_98D.RadialHeader()
_WSR98D_RADIAL_HEADER_KEYS = tuple(field_name for field_name, _ in _WSR98D_RADIAL_HEADER)
_WSR98D_RADIAL_HEADER_STRUCT = struct.Struct("<" + "".join(field_type for _, field_type in _WSR98D_RADIAL_HEADER))
_WSR98D_MOMENT_HEADER = dtype_98D.RadialData()
_WSR98D_MOMENT_HEADER_KEYS = tuple(field_name for field_name, _ in _WSR98D_MOMENT_HEADER)
_WSR98D_MOMENT_HEADER_STRUCT = struct.Struct("<" + "".join(field_type for _, field_type in _WSR98D_MOMENT_HEADER))


def _decode_wsr98d_resolution(raw_value):
    """Decode meter resolution fields, including the vendor 62.5 m fixed-point encoding."""
    raw_value = int(raw_value)
    if raw_value >= 32768:
        return (raw_value - 32768) / 100.0
    return float(raw_value)


def _encode_wsr98d_resolution(resolution_m):
    """Encode a meter resolution back into the on-disk WSR98D integer representation."""
    resolution_m = float(resolution_m)
    rounded_meter = round(resolution_m)
    if abs(resolution_m - rounded_meter) < 1.0e-6:
        return int(rounded_meter)
    return int(round(resolution_m * 100.0 + 32768.0))


def _ensure_ppi_prd_volume(prd):
    scan_type = str(np.asarray(prd.scan_info["scan_type"].values).item())
    if scan_type != "ppi":
        raise ValueError("WSR98D export only supports ppi PRD volumes.")


def _ensure_output_path(filename, overwrite=False):
    path = os.path.abspath(filename)
    if os.path.exists(path) and not overwrite:
        raise FileExistsError("Output file already exists: %s" % path)
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    return path


def _python_datetime(value):
    array = np.asarray(value)
    if np.issubdtype(array.dtype, np.datetime64):
        return array.astype("datetime64[us]").tolist()
    if isinstance(value, datetime.datetime):
        return value
    raise TypeError("Unable to convert %r to datetime." % (value,))


def _epoch_parts(dt):
    epoch = datetime.datetime(1970, 1, 1, tzinfo=dt.tzinfo)
    delta = dt - epoch
    total_seconds = int(delta.total_seconds())
    microseconds = int(delta.microseconds)
    return total_seconds, microseconds


def _normalize_ascii_bytes(value, length, fallback):
    text = (str(value).strip() if value is not None else "") or fallback
    return text[:length].ljust(length).encode("ascii", "ignore")


def _available_wsr98d_fields(prd, field_names=None, strict=True):
    available = set(prd.available_fields(range_mode=None))
    if field_names is None:
        selected = [name for name in WSR98D_WRITE_FIELD_SPECS if name in available]
    else:
        selected = []
        for name in field_names:
            if name not in WSR98D_WRITE_FIELD_SPECS:
                if strict:
                    raise ValueError("Field %s is not supported by WSR98D export." % name)
                continue
            if name not in available:
                if strict:
                    raise KeyError(name)
                continue
            selected.append(name)
    if not selected:
        raise ValueError("No exportable WSR98D fields are available.")
    return selected


def _range_geometry(range_values):
    range_values = np.asarray(range_values, dtype=np.float64)
    if range_values.size == 0:
        raise ValueError("Export range vector cannot be empty.")
    if range_values.size > 1:
        spacing = float(range_values[1] - range_values[0])
        if spacing <= 0.0:
            raise ValueError("Range spacing must be positive.")
    else:
        spacing = float(range_values[0]) if float(range_values[0]) > 0.0 else 1.0
    return int(round(float(range_values[0]))), int(round(spacing))


def _select_resolution_range(prd, sweep, candidates):
    for field_name in candidates:
        source_name = field_name
        if hasattr(prd, "resolve_field_name"):
            source_name = prd.resolve_field_name(field_name, sweep=sweep, range_mode=None, required=False)
        if source_name is not None and source_name in prd.available_fields(sweep=sweep, range_mode=None):
            field = prd.get_sweep_field(sweep, source_name, range_mode=None)
            return np.asarray(field["range"].values, dtype=np.float64)
    return np.asarray(prd.fields[sweep]["range"].values, dtype=np.float64)


def _encode_quantized(values, scale, offset, bin_length, missing_code=3):
    max_code = 255 if bin_length == 1 else 65535
    codes = np.round(np.asarray(values, dtype=np.float64) * float(scale) + float(offset))
    codes[~np.isfinite(values)] = missing_code
    codes = np.clip(codes.astype(np.int64), 0, max_code)
    if bin_length == 1:
        return codes.astype(np.uint8).tobytes()
    if bin_length == 2:
        return codes.astype("<u2").tobytes()
    raise ValueError("WSR98D bin length must be 1 or 2.")


class WSR98DBaseData(object):
    """Decode standard WSR-98D dual-polarization base data."""

    def __init__(self, filename, station_lon=None, station_lat=None, station_alt=None):
        """
        :param filename:  radar basedata filename
        :param station_lon:  radar station longitude //units: degree east
        :param station_lat:  radar station latitude //units:degree north
        :param station_alt:  radar station altitude //units: meters
        """
        super(WSR98DBaseData, self).__init__()
        self.filename = filename
        self.station_lon = station_lon
        self.station_lat = station_lat
        self.station_alt = station_alt
        self.fid = _prepare_for_read(self.filename)
        self._check_standard_basedata()
        self.header = self._parse_BaseDataHeader()
        self.radial, self._status, self._azimuth, self._elevation, self._seconds, self._microseconds = \
            self._parse_radial()
        self._scan_time = None
        self._nyquist_velocity = None
        self._unambiguous_range = None
        self.nrays = len(self.radial)
        self.sweep_start_ray_index = np.where((self._status == 0) | (self._status == 3))[0]
        self.sweep_end_ray_index = np.where((self._status == 2) | (self._status == 4))[0]
        self.nsweeps = len(self.sweep_start_ray_index)
        self.fid.close()

    def _check_standard_basedata(self):
        """
        :param fid: file fid
        :return:
        """
        if _read_exact(self.fid, 4, "WSR98D file signature") != b'RSTM':
            raise ValueError("File is not a standard WSR-98D file.")
        self.fid.seek(0, 0)
        return

    def _parse_BaseDataHeader(self):
        BaseDataHeader = {}
        fixed_buf = _read_exact(self.fid, dtype_98D.CutConfigurationBlockPos, "WSR98D fixed header")

        BaseDataHeader['GenericHeader'], _ = _unpack_from_buf(fixed_buf, \
                                                              dtype_98D.GenericHeaderBlockPos,
                                                              dtype_98D.BaseDataHeader['GenericHeaderBlock'])
        BaseDataHeader['SiteConfig'], _ = _unpack_from_buf(fixed_buf, \
                                                           dtype_98D.SiteConfigurationBlockPos,
                                                           dtype_98D.BaseDataHeader['SiteConfigurationBlock'])
        BaseDataHeader['TaskConfig'], _ = _unpack_from_buf(fixed_buf,
                                                           dtype_98D.TaskConfigurationBlockPos,
                                                           dtype_98D.BaseDataHeader['TaskConfigurationBlock'])
        cut_count = _validate_count(
            "WSR98D CutNumber",
            BaseDataHeader['TaskConfig']['CutNumber'],
            minimum=1,
            maximum=MAX_WSR98D_SWEEPS,
        )
        cut_buf = _read_exact(
            self.fid,
            dtype_98D.CutConfigurationBlockSize * cut_count,
            "WSR98D cut configuration block",
        )
        BaseDataHeader['CutConfig'] = np.frombuffer(cut_buf, dtype_98D.BaseDataHeader['CutConfigurationBlock'])
        return BaseDataHeader

    def _parse_radial(self):
        radial = []
        azimuth = []
        elevation = []
        status = []
        seconds = []
        microseconds = []
        buf = _read_all(
            self.fid,
            "WSR98D radial payload",
            max_bytes=WSR98D_MAX_DECODED_PAYLOAD_BYTES,
        )
        buffer_view = memoryview(buf)
        pos = 0
        size = len(buf)
        header_size = _WSR98D_RADIAL_HEADER_STRUCT.size
        radial_keys = _WSR98D_RADIAL_HEADER_KEYS
        radial_struct = _WSR98D_RADIAL_HEADER_STRUCT
        while pos + header_size <= size:
            radial_values = radial_struct.unpack_from(buf, pos)
            pos += header_size
            RadialDict = dict(zip(radial_keys, radial_values))
            moment_num = _validate_count(
                "WSR98D MomentNumber",
                RadialDict['MomentNumber'],
                minimum=0,
                maximum=MAX_WSR98D_MOMENTS_PER_RADIAL,
            )
            RadialDict['fields'], pos = self._parse_radial_single(buffer_view, pos, moment_num)
            radial.append(RadialDict)
            azimuth.append(RadialDict['Azimuth'])
            elevation.append(RadialDict['Elevation'])
            status.append(RadialDict['RadialState'])
            seconds.append(RadialDict['Seconds'])
            microseconds.append(RadialDict['MicroSeconds'])
        if pos != size:
            raise ValueError("WSR98D radial payload contains trailing truncated data.")
        return radial, np.asarray(status), np.asarray(azimuth), np.asarray(elevation), \
            np.asarray(seconds), np.asarray(microseconds)

    def _parse_radial_single(self, buf, pos, moment_num):
        radial_var = {}
        moment_header_size = _WSR98D_MOMENT_HEADER_STRUCT.size
        moment_header_struct = _WSR98D_MOMENT_HEADER_STRUCT
        moment_header_keys = _WSR98D_MOMENT_HEADER_KEYS
        flag_to_product = dtype_98D.flag2Product
        buffer_size = len(buf)
        for _ in range(moment_num):
            if pos + moment_header_size > buffer_size:
                raise ValueError("WSR98D moment header is truncated or malformed.")
            moment_header_values = moment_header_struct.unpack_from(buf, pos)
            Momheader = dict(zip(moment_header_keys, moment_header_values))
            pos += moment_header_size
            data_len = _validate_count(
                "WSR98D moment length",
                Momheader['Length'],
                minimum=0,
                maximum=MAX_WSR98D_MOMENT_DATA_BYTES,
            )
            if pos + data_len > buffer_size:
                raise ValueError("WSR98D moment payload extends beyond the available buffer.")
            Data_buf = buf[pos:pos + data_len]
            pos += data_len
            if Momheader['BinLength'] not in (1, 2):
                raise ValueError("WSR98D moment bin length must be 1 or 2 bytes.")
            if data_len % int(Momheader['BinLength']) != 0:
                raise ValueError("WSR98D moment payload length does not match the bin length.")
            if int(Momheader['Scale']) == 0:
                raise ValueError("WSR98D moment scale cannot be zero.")
            if Momheader['BinLength'] == 1:
                dat_tmp = np.frombuffer(Data_buf, dtype=np.uint8, count=data_len)
            else:
                dat_tmp = np.frombuffer(Data_buf, dtype=np.uint16, count=data_len // 2)
            if Momheader['DataType'] <= 35:
                dat = np.full(dat_tmp.shape, np.nan, dtype=np.float32)
                valid = dat_tmp >= 5
                if np.any(valid):
                    dat[valid] = (dat_tmp[valid].astype(np.float32) - Momheader['Offset']) / Momheader['Scale']
                field_name = flag_to_product.get(int(Momheader['DataType']))
                if field_name is not None:
                    radial_var[field_name] = dat
        return radial_var, pos

    def get_nyquist_velocity(self):
        """Return the per-ray Nyquist velocity."""
        if self._nyquist_velocity is None:
            self._nyquist_velocity = np.repeat(self.header['CutConfig']['NyquistSpeed'], self.get_rays_per_sweep())
        return self._nyquist_velocity

    def get_unambiguous_range(self):
        """Return the per-ray unambiguous range."""
        if self._unambiguous_range is None:
            self._unambiguous_range = np.repeat(self.header['CutConfig']['MaximumRange'], self.get_rays_per_sweep())
        return self._unambiguous_range

    def get_scan_time(self):
        """Return the acquisition time for each ray."""
        if self._scan_time is None:
            self._scan_time = np.array(
                [julian2date_SEC(sec, usec) for sec, usec in zip(self._seconds, self._microseconds)]
            )
        return self._scan_time

    def get_sweep_end_ray_index(self):
        """Return the inclusive end index of each sweep."""
        return self.sweep_end_ray_index

    def get_sweep_start_ray_index(self):
        """Return the start index of each sweep."""
        return self.sweep_start_ray_index

    def get_rays_per_sweep(self):
        """Return the number of rays in each sweep."""
        return self.sweep_end_ray_index - self.sweep_start_ray_index + 1

    def get_azimuth(self):
        """Return the azimuth angle for each ray."""
        return self._azimuth

    def get_elevation(self):
        """Return the elevation angle for each ray."""
        return self._elevation

    def get_latitude_longitude_altitude_frequency(self):
        """Return latitude, longitude, altitude, and radar frequency."""
        lat, lon, alt, frequency = self.header['SiteConfig']['Latitude'], self.header['SiteConfig']['Longitude'], \
               self.header['SiteConfig']['Height'], self.header['SiteConfig']['Frequency'] / 1000.
        if self.station_lon is not None:
            lon = self.station_lon
        if self.station_lat is not None:
            lat = self.station_lat
        if self.station_alt is not None:
            alt = self.station_alt
        return lat, lon, alt, frequency

    def get_scan_type(self):
        if self.header['TaskConfig']['ScanType'] in [0, 1]:
            return "ppi"
        elif self.header['TaskConfig']['ScanType'] in [2, 5]:
            return "rhi"
        elif self.header['TaskConfig']['ScanType'] in [3, 4]:
            return "sector"
        else:
            return "other"
    def get_sitename(self):
        return (self.header['SiteConfig']['SiteName']).decode('UTF-8', 'ignore').strip().strip(b'\x00'.decode())


class WSR98D2NRadar(object):
    """Bridge from raw WSR98D data to an NRadar object."""

    _REFLECTIVITY_LIKE_FIELDS = ("dBZ", "Zc", "dBT")
    _VELOCITY_LIKE_FIELDS = ("V", "Vc", "W", "Wc")
    _RANGE_REFERENCE_FIELDS = ("V", "Vc", "W", "Wc", "dBZ", "Zc", "dBT")

    def __init__(self, WSR98D):
        self.WSR98D = WSR98D
        self.flag_match = np.all(self.WSR98D.header['CutConfig']['LogResolution'] == \
                       self.WSR98D.header['CutConfig']['DopplerResolution'])
        self._dbz_index_cache = {}
        self._sweep_descriptors = self._describe_raw_sweeps()
        self._paired_sweeps = self._pair_split_sweeps()
        self.dBZ_index_alone = np.array([source for source, _ in self._paired_sweeps], dtype=np.int32)
        self.v_index_alone = np.array([target for _, target in self._paired_sweeps], dtype=np.int32)

        for index_with_dbz, index_with_v in self._paired_sweeps:
            self.interp_dBZ(index_with_dbz, index_with_v)

        keep_mask = np.ones(self.WSR98D.nrays, dtype=bool)
        keep_mask[self.get_remove_radial_indices()] = False
        self.radial = [iray for iray, keep in zip(self.WSR98D.radial, keep_mask) if keep]
        self._azimuth = self.WSR98D.get_azimuth()[keep_mask]
        self._elevation = self.WSR98D.get_elevation()[keep_mask]
        self._seconds = self.WSR98D._seconds[keep_mask]
        self._microseconds = self.WSR98D._microseconds[keep_mask]
        self._scan_time = None
        self._nyquist_velocity = None
        self._unambiguous_range = None

        status = self.WSR98D._status[keep_mask]
        self.sweep_start_ray_index = np.where((status == 0) | (status == 3))[0]
        self.sweep_end_ray_index = np.where((status == 2) | (status == 4))[0]
        self.nsweeps = len(self.sweep_start_ray_index)
        self.nrays = len(self.radial)
        self.scan_type = self.WSR98D.get_scan_type()
        self.latitude, self.longitude, self.altitude, self.frequency = \
            self.WSR98D.get_latitude_longitude_altitude_frequency()
        self.header = self.WSR98D.header
        self.header['CutConfig'] = np.delete(self.header['CutConfig'], self.dBZ_index_alone)
        self.bins_per_sweep = self.get_nbins_per_sweep()
        self.max_bins = self.bins_per_sweep.max()
        self.range = self.get_range_per_radial(self.max_bins)
        self.azimuth = self.get_azimuth()
        self.elevation = self.get_elevation()
        self.extended_fields = self._build_extended_fields()
        self.fields = self._get_fields()
        self.sitename = self.WSR98D.get_sitename()

    def _describe_raw_sweeps(self):
        """Classify raw sweeps as reflectivity-only, Doppler-only, or combined."""
        descriptors = []
        for sweep_index, start_idx in enumerate(self.WSR98D.sweep_start_ray_index):
            fields = tuple(self.WSR98D.radial[start_idx]["fields"].keys())
            field_names = set(fields)
            has_reflectivity = any(name in field_names for name in self._REFLECTIVITY_LIKE_FIELDS)
            has_doppler = any(name in field_names for name in self._VELOCITY_LIKE_FIELDS)
            descriptors.append(
                {
                    "index": sweep_index,
                    "fields": fields,
                    "elevation": float(self.WSR98D.header["CutConfig"]["Elevation"][sweep_index]),
                    "reflectivity_only": has_reflectivity and not has_doppler,
                    "doppler_only": has_doppler and not has_reflectivity,
                    "combined": has_reflectivity and has_doppler,
                }
            )
        return descriptors

    def _pair_split_sweeps(self, elevation_tolerance=0.5):
        """Pair reflectivity-only sweeps with the closest compatible Doppler-only sweeps."""
        pending_reflectivity = []
        paired = []
        for descriptor in self._sweep_descriptors:
            sweep_index = descriptor["index"]
            if descriptor["reflectivity_only"]:
                pending_reflectivity.append(sweep_index)
                continue
            if not descriptor["doppler_only"]:
                continue
            match = None
            for candidate in reversed(pending_reflectivity):
                candidate_descriptor = self._sweep_descriptors[candidate]
                if abs(candidate_descriptor["elevation"] - descriptor["elevation"]) <= float(elevation_tolerance):
                    match = candidate
                    break
            if match is None:
                continue
            pending_reflectivity.remove(match)
            paired.append((match, sweep_index))
        return paired

    def get_remove_radial_indices(self):
        """Return the radial indices removed after sweep alignment."""
        dBZ_alone = self.dBZ_index_alone
        index_romove = []
        for isweep in dBZ_alone:
            index_romove.extend(range(self.WSR98D.sweep_start_ray_index[isweep], \
                                      self.WSR98D.sweep_end_ray_index[isweep] + 1))
        return index_romove

    def get_reomve_radial_num(self):
        """Backward-compatible alias for ``get_remove_radial_indices``."""
        return self.get_remove_radial_indices()

    def get_v_idx(self):
        """Return sweep indices that contain Doppler moments without reflectivity."""
        return np.array(
            [descriptor["index"] for descriptor in self._sweep_descriptors if descriptor["doppler_only"]],
            dtype=np.int32,
        )

    def get_dbz_idx(self):
        """Return sweep indices that contain reflectivity without Doppler moments."""
        return np.array(
            [descriptor["index"] for descriptor in self._sweep_descriptors if descriptor["reflectivity_only"]],
            dtype=np.int32,
        )

    def interp_VCP26(self, dBZ_sweep_index, V_sweep_index):
        """
        Handle the VCP26/VCP27 case where reflectivity and velocity sweeps do not align.
        :param dBZ_sweep_index: reflectivity-only sweep indices.
        :param V_sweep_index: Doppler-only sweep indices.
        """
        same_sweeps = min(len(dBZ_sweep_index), len(V_sweep_index))
        for isweep in range(same_sweeps):
            self.interp_dBZ(dBZ_sweep_index[isweep], V_sweep_index[isweep])

    def interp_dBZ(self, field_with_dBZ_num, field_without_dBZ_num):
        """
        Copy reflectivity rays onto the matching Doppler-only sweep by azimuth.
        :param field_with_dBZ_num: source sweep index, 0-based.
        :param field_without_dBZ_num: target sweep index, 0-based.
        """
        azimuth = self.WSR98D.get_azimuth()
        dbz_az = azimuth[self.WSR98D.sweep_start_ray_index[field_with_dBZ_num]: \
                         self.WSR98D.sweep_end_ray_index[field_with_dBZ_num] + 1]
        v_az = azimuth[self.WSR98D.sweep_start_ray_index[field_without_dBZ_num]: \
                       self.WSR98D.sweep_end_ray_index[field_without_dBZ_num] + 1]
        dbz_idx = np.argmin(np.abs(dbz_az.reshape(-1, 1) - v_az.reshape(1, -1)), axis=0) + \
                  self.WSR98D.sweep_start_ray_index[field_with_dBZ_num]

        v_idx = np.arange(self.WSR98D.sweep_start_ray_index[field_without_dBZ_num], \
                          self.WSR98D.sweep_end_ray_index[field_without_dBZ_num] + 1)
        keys = self.WSR98D.radial[self.WSR98D.sweep_start_ray_index[field_with_dBZ_num]]['fields'].keys()
        for ind_dbz, ind_v in zip(dbz_idx, v_idx):
            for ikey in keys:
                self.WSR98D.radial[ind_v]["fields"][ikey] = self.WSR98D.radial[ind_dbz]["fields"][ikey]

    def get_azimuth(self):
        """Return the azimuth angle for each retained ray."""
        return self._azimuth

    def get_elevation(self):
        """Return the elevation angle for each retained ray."""
        return np.where(self._elevation > 180, self._elevation - 360, self._elevation)

    def get_rays_per_sweep(self):
        """Return the number of rays in each retained sweep."""
        return self.sweep_end_ray_index - self.sweep_start_ray_index + 1

    def get_scan_time(self):
        """Return the acquisition time for each retained ray."""
        if self._scan_time is None:
            self._scan_time = np.array(
                [julian2date_SEC(sec, usec) for sec, usec in zip(self._seconds, self._microseconds)]
            )
        return self._scan_time

    def get_nyquist_velocity(self):
        """Return the per-ray Nyquist velocity."""
        if self._nyquist_velocity is None:
            self._nyquist_velocity = np.repeat(self.header['CutConfig']['NyquistSpeed'], self.get_rays_per_sweep())
        return self._nyquist_velocity

    def get_unambiguous_range(self):
        """Return the per-ray unambiguous range."""
        if self._unambiguous_range is None:
            self._unambiguous_range = np.repeat(self.header['CutConfig']['MaximumRange'], self.get_rays_per_sweep())
        return self._unambiguous_range

    def get_sweep_end_ray_index(self):
        """Return the inclusive end index of each retained sweep."""
        return self.sweep_end_ray_index

    def get_sweep_start_ray_index(self):
        """Return the start index of each retained sweep."""
        return self.sweep_start_ray_index

    def get_nbins_per_sweep(self):
        """Return the Doppler gate count for each retained sweep."""
        bins = []
        for idx in self.sweep_start_ray_index:
            radial_fields = self.radial[idx]["fields"]
            selected = None
            for field_name in self._RANGE_REFERENCE_FIELDS:
                selected = radial_fields.get(field_name)
                if selected is not None and selected.size:
                    bins.append(int(selected.size))
                    break
            else:
                sizes = [int(values.size) for values in radial_fields.values() if getattr(values, "size", 0)]
                bins.append(max(sizes) if sizes else 0)
        return np.asarray(bins, dtype=np.int32)

    def get_range_per_radial(self, length, sweep=0):
        """Return Doppler gate-center ranges for a radial of ``length`` bins."""
        Resolution = _decode_wsr98d_resolution(self.header['CutConfig']['DopplerResolution'][int(sweep)])
        return np.linspace(Resolution, Resolution * length, length)

    def get_dbz_range_per_radial(self, length, sweep=0):
        """Return reflectivity gate-center ranges for a radial of ``length`` bins."""
        Resolution = _decode_wsr98d_resolution(self.header['CutConfig']['LogResolution'][int(sweep)])
        return np.linspace(Resolution, Resolution * length, length)

    def _get_fields(self):
        """Assemble the retained fields into dense 2-D arrays."""
        field_keys = []
        seen = set()
        for radial in self.radial:
            for key in radial["fields"].keys():
                if key not in seen:
                    seen.add(key)
                    field_keys.append(key)
        field_keys = tuple(field_keys)
        fields = {ikey: np.full((self.nrays, self.max_bins), np.nan, dtype=np.float64) for ikey in field_keys}
        for sweep_idx, (start, end) in enumerate(zip(self.sweep_start_ray_index, self.sweep_end_ray_index)):
            for iray in range(start, end + 1):
                radial_fields = self.radial[iray]['fields']
                for ikey in field_keys:
                    fields[ikey][iray, :] = self._add_or_del_field(
                        radial_fields.get(ikey),
                        ikey,
                        self.flag_match,
                        sweep_idx=sweep_idx,
                    )
        return fields

    def _build_extended_fields(self):
        """Collect native-range reflectivity for sweeps longer than the aligned grid."""
        extended_sweeps = {}
        scan_time = self.get_scan_time()
        for sweep_idx, (start, end, aligned_bins) in enumerate(
            zip(self.sweep_start_ray_index, self.sweep_end_ray_index, self.bins_per_sweep)
        ):
            native_bins = max(
                (
                    self.radial[iray]["fields"].get("dBZ", np.empty((0,), dtype=np.float32)).size
                    for iray in range(start, end + 1)
                ),
                default=0,
            )
            if native_bins <= int(aligned_bins):
                continue

            data = np.full((end - start + 1, native_bins), np.nan, dtype=np.float32)
            for out_idx, iray in enumerate(range(start, end + 1)):
                dat_ray = self.radial[iray]["fields"].get("dBZ")
                if dat_ray is None or dat_ray.size == 0:
                    continue
                copy_len = min(dat_ray.size, native_bins)
                data[out_idx, :copy_len] = dat_ray[:copy_len].astype(np.float32, copy=False)

            extended_sweeps[sweep_idx] = {
                "data": data,
                "range": self.get_dbz_range_per_radial(native_bins, sweep=sweep_idx).astype(np.float32, copy=False),
                "time": np.asarray(scan_time[start:end + 1]),
                "azimuth": np.asarray(self._azimuth[start:end + 1], dtype=np.float32),
                "elevation": np.asarray(self._elevation[start:end + 1], dtype=np.float32),
                "aligned_bins": int(aligned_bins),
            }

        return {"dBZ": extended_sweeps} if extended_sweeps else {}

    def _get_dbz_resample_index(self, source_length, sweep_idx=0):
        cache_key = (int(sweep_idx), int(source_length))
        cached = self._dbz_index_cache.get(cache_key)
        if cached is not None:
            return cached

        source = self.get_dbz_range_per_radial(source_length, sweep=sweep_idx)
        target = self.range
        right = np.searchsorted(source, target, side="left")
        left = np.clip(right - 1, 0, source_length - 1)
        right = np.clip(right, 0, source_length - 1)
        choose_right = np.abs(source[right] - target) < np.abs(target - source[left])
        nearest = np.where(choose_right, right, left)
        valid = (target >= source[0]) & (target <= source[-1])
        self._dbz_index_cache[cache_key] = (valid, nearest)
        return valid, nearest

    def _add_or_del_field(self, dat_ray, key, flag_match=True, sweep_idx=0):
        """
        Normalize a radial field to the common Doppler-aligned range grid.
        :param dat_ray: radial field data.
        :param key: field name.
        :param flag_match: whether Doppler and reflectivity resolutions match.
        """
        length = self.max_bins
        if dat_ray is None:
            return np.full((length,), np.nan)

        if not flag_match and key == "dBZ":
            valid, nearest = self._get_dbz_resample_index(dat_ray.size, sweep_idx=sweep_idx)
            out = np.full((length,), np.nan)
            out[valid] = dat_ray[nearest[valid]]
            return out

        assert dat_ray.ndim == 1, "check dat_ray"
        if dat_ray.size >= length:
            return dat_ray[:length]
        else:
            out = np.full((length,), np.nan)
            out[:dat_ray.size] = dat_ray
            return out

    def get_nradar_nyquist_speed(self):
        """array shape (nsweeps)"""
        return self.header['CutConfig']['NyquistSpeed']

    def get_NRadar_nyquist_speed(self):
        """Backward-compatible alias for ``get_nradar_nyquist_speed``."""
        return self.get_nradar_nyquist_speed()

    def get_nradar_unambiguous_range(self):
        """array shape (nsweeps)"""
        return self.header['CutConfig']['MaximumRange']

    def get_NRadar_unambiguous_range(self):
        """Backward-compatible alias for ``get_nradar_unambiguous_range``."""
        return self.get_nradar_unambiguous_range()

    def get_fixed_angle(self):
        if self.scan_type == "rhi":
            return self.header['CutConfig']['Azimuth']
        else:
            return self.header['CutConfig']['Elevation']

    def to_prd(self, effective_earth_radius=None):
        """Build an ``NRadar.PRD`` object from the decoded WSR-98D volume."""
        return PRD(fields=self.fields, scan_type=self.scan_type, time=self.get_scan_time(), \
                          range=self.range, azimuth=self.azimuth, elevation=self.elevation, latitude=self.latitude, \
                          longitude=self.longitude, altitude=self.altitude,
                          sweep_start_ray_index=self.sweep_start_ray_index, \
                          sweep_end_ray_index=self.sweep_end_ray_index, fixed_angle=self.get_fixed_angle(), \
                          bins_per_sweep=self.bins_per_sweep, nyquist_velocity=self.get_nradar_nyquist_speed(), \
                          frequency=self.frequency, unambiguous_range=self.get_nradar_unambiguous_range(), \
                          nrays=self.nrays, nsweeps=self.nsweeps, sitename = self.sitename,
                          pyart_radar=None, effective_earth_radius=effective_earth_radius,
                          extended_fields=self.extended_fields,
                          metadata={"original_container": "WSR98D", "radar_name": "WSR98D"})

    def ToPRD(self, effective_earth_radius=None):
        """Backward-compatible alias for ``to_prd``."""
        return self.to_prd(effective_earth_radius=effective_earth_radius)

    def to_pyart_radar(self, effective_earth_radius=None, **kwargs):
        """Export the decoded WSR-98D volume through the PRD Py-ART adapter."""
        return self.to_prd(effective_earth_radius=effective_earth_radius).to_pyart_radar(**kwargs)

    def ToPyartRadar(self, effective_earth_radius=None, **kwargs):
        """Backward-compatible alias for ``to_pyart_radar``."""
        return self.to_pyart_radar(effective_earth_radius=effective_earth_radius, **kwargs)

    def _get_instrument_parameters(self):
        """ Return a dictionary containing instrument parameters. """

        # pulse width
        pulse_width = get_metadata('pulse_width')
        pulse_width['data'] = np.array([self.header["TaskConfig"]['PulseWidth']/10**9,], dtype='float32')  # nanosec->sec
        # assume that the parameters in the first ray represent the beam widths,
        # bandwidth and frequency in the entire volume
        wavelength_hz = self.frequency * 10 ** 9
        # radar_beam_width_h
        radar_beam_width_h = get_metadata('radar_beam_width_h')
        radar_beam_width_h['data'] = np.array([self.header['SiteConfig']['BeamWidthHori'], ], dtype='float32')
        # radar_beam_width_v
        radar_beam_width_v = get_metadata('radar_beam_width_v')
        radar_beam_width_v['data'] = np.array([self.header['SiteConfig']['BeamWidthVert'], ], dtype='float32')
        # frequency
        frequency = get_metadata('frequency')
        frequency['data'] = np.array([wavelength_hz], dtype='float32')
        instrument_parameters = {
            'pulse_width': pulse_width,
            'radar_beam_width_h': radar_beam_width_h,
            'radar_beam_width_v': radar_beam_width_v,
            'frequency': frequency, }

        # nyquist velocity if defined
        nyquist_velocity = get_metadata('nyquist_velocity')
        nyquist_velocity['data'] = self.get_nyquist_velocity()
        instrument_parameters['nyquist_velocity'] = nyquist_velocity
        return instrument_parameters


def write_wsr98d(
    prd,
    filename,
    field_names=None,
    strict=True,
    overwrite=False,
    site_code=None,
    site_name=None,
    task_name=None,
    frequency=None,
    rda_version=722951,
    radar_type=1,
):
    """Write a PRD volume to standard WSR-98D base-data format."""
    _ensure_ppi_prd_volume(prd)
    path = _ensure_output_path(filename, overwrite=overwrite)
    selected_fields = _available_wsr98d_fields(prd, field_names=field_names, strict=strict)

    latitude = float(prd.scan_info["latitude"].values)
    longitude = float(prd.scan_info["longitude"].values)
    altitude = float(prd.scan_info["altitude"].values)
    frequency_ghz = float(frequency) if frequency is not None else float(prd.scan_info["frequency"].values)
    beam_width = np.asarray(prd.scan_info["beam_width"].values, dtype=np.float64)
    fixed_angle = np.asarray(prd.scan_info["fixed_angle"].values, dtype=np.float64)
    rays_per_sweep = np.asarray(prd.scan_info["rays_per_sweep"].values, dtype=np.int32)
    nyquist_velocity = np.asarray(prd.scan_info["nyquist_velocity"].values, dtype=np.float64)
    unambiguous_range = np.asarray(prd.scan_info["unambiguous_range"].values, dtype=np.float64)
    start_time = _python_datetime(prd.scan_info["start_time"].values)

    site_code_bytes = _normalize_ascii_bytes(site_code, 8, "Z0000")
    site_name_bytes = _normalize_ascii_bytes(site_name, 32, getattr(prd, "sitename", "pycwr"))
    task_name_bytes = _normalize_ascii_bytes(task_name, 32, "PYCWR")

    fmt_generic_header = "<" + "".join(item[1] for item in dtype_98D.BaseDataHeader["GenericHeaderBlock"])
    fmt_site_configuration = "<" + "".join(item[1] for item in dtype_98D.BaseDataHeader["SiteConfigurationBlock"])
    fmt_task_configuration = "<" + "".join(item[1] for item in dtype_98D.BaseDataHeader["TaskConfigurationBlock"])
    fmt_radial_header = "<" + "".join(item[1] for item in dtype_98D.RadialHeader())
    fmt_moment_header = "<" + "".join(item[1] for item in dtype_98D.RadialData())

    generic_header = struct.pack(fmt_generic_header, 1297371986, -1, -1, 1, -1, b"")
    site_configuration = struct.pack(
        fmt_site_configuration,
        site_code_bytes,
        site_name_bytes,
        latitude,
        longitude,
        int(round(altitude)),
        int(round(altitude)),
        frequency_ghz * 1000.0,
        float(np.nanmean(beam_width)) if beam_width.size else 1.0,
        float(np.nanmean(beam_width)) if beam_width.size else 1.0,
        int(rda_version),
        int(radar_type),
        b"",
    )
    task_configuration = struct.pack(
        fmt_task_configuration,
        task_name_bytes,
        b"",
        -1,
        0,
        -1,
        int((_python_datetime(start_time) - datetime.datetime(1970, 1, 1, tzinfo=start_time.tzinfo)).total_seconds()),
        int(prd.nsweeps),
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        b"",
    )

    cut_config = np.zeros(int(prd.nsweeps), dtype=dtype_98D.BaseDataHeader["CutConfigurationBlock"])
    for sweep in range(int(prd.nsweeps)):
        log_range = _select_resolution_range(prd, sweep, ("dBZ", "Zc", "dBT"))
        dop_range = _select_resolution_range(prd, sweep, ("V", "W", "ZDR", "CC", "PhiDP", "KDP", "SNRH", "SNRV"))
        _, log_resolution = _range_geometry(log_range)
        _, dop_resolution = _range_geometry(dop_range)
        cut_config["ProcessMode"][sweep] = 1
        cut_config["WaveForm"][sweep] = 0
        cut_config["PRF_1"][sweep] = 0.0
        cut_config["PRF_2"][sweep] = 0.0
        cut_config["UnfoldMode"][sweep] = 0
        cut_config["Azimuth"][sweep] = 0.0
        cut_config["Elevation"][sweep] = fixed_angle[sweep]
        cut_config["StartAngle"][sweep] = 0.0
        cut_config["EndAngle"][sweep] = 360.0
        cut_config["AngleResolution"][sweep] = float(beam_width[sweep]) if beam_width.size > sweep else 1.0
        cut_config["ScanSpeed"][sweep] = 0.0
        cut_config["LogResolution"][sweep] = _encode_wsr98d_resolution(log_resolution)
        cut_config["DopplerResolution"][sweep] = _encode_wsr98d_resolution(dop_resolution)
        cut_config["MaximumRange"][sweep] = int(round(unambiguous_range[sweep])) if unambiguous_range.size > sweep else int(round(log_range[-1]))
        cut_config["MaximumRange2"][sweep] = 0
        cut_config["StartRange"][sweep] = int(round(min(log_range[0], dop_range[0])))
        cut_config["Sample_1"][sweep] = int(rays_per_sweep[sweep])
        cut_config["Sample_2"][sweep] = 0
        cut_config["PhaseMode"][sweep] = 0
        cut_config["AtmosphericLoss"][sweep] = 0.0
        cut_config["NyquistSpeed"][sweep] = float(nyquist_velocity[sweep]) if nyquist_velocity.size > sweep else 0.0
        cut_config["MomentsMask"][sweep] = 0
        cut_config["MomentsSizeMask"][sweep] = 0

    with open(path, "wb") as handle:
        handle.write(generic_header)
        handle.write(site_configuration)
        handle.write(task_configuration)
        handle.write(cut_config.tobytes())

        sequence_number = 1
        for sweep in range(int(prd.nsweeps)):
            sweep_dataset = prd.fields[sweep]
            rays = int(rays_per_sweep[sweep])
            for iray in range(rays):
                azimuth = float(sweep_dataset["azimuth"].values[iray])
                elevation = float(sweep_dataset["elevation"].values[iray])
                dt = _python_datetime(sweep_dataset["time"].values[iray])
                seconds, microseconds = _epoch_parts(dt)
                radial_state = 3 if (sweep == 0 and iray == 0) else 4 if (sweep == int(prd.nsweeps) - 1 and iray == rays - 1) else 0 if iray == 0 else 2 if iray == rays - 1 else 1

                moment_headers = []
                moment_payloads = []
                for field_name in selected_fields:
                    if field_name not in prd.available_fields(sweep=sweep, range_mode=None):
                        continue
                    spec = WSR98D_WRITE_FIELD_SPECS[field_name]
                    field = prd.get_sweep_field(sweep, field_name, range_mode=None)
                    values = np.asarray(field.values[iray], dtype=np.float32)
                    data_bytes = _encode_quantized(
                        values,
                        spec["scale"],
                        spec["offset"],
                        spec["bin_length"],
                    )
                    moment_headers.append(
                        struct.pack(
                            fmt_moment_header,
                            spec["data_type"],
                            spec["scale"],
                            spec["offset"],
                            spec["bin_length"],
                            0,
                            len(data_bytes),
                            b"",
                        )
                    )
                    moment_payloads.append(data_bytes)

                length_of_data = sum(len(header) + len(payload) for header, payload in zip(moment_headers, moment_payloads))
                radial_header = struct.pack(
                    fmt_radial_header,
                    radial_state,
                    0,
                    sequence_number,
                    iray + 1,
                    sweep + 1,
                    azimuth,
                    elevation,
                    seconds,
                    microseconds,
                    length_of_data,
                    len(moment_headers),
                    b"",
                )
                handle.write(radial_header)
                for header, payload in zip(moment_headers, moment_payloads):
                    handle.write(header)
                    handle.write(payload)
                sequence_number += 1
    return path
