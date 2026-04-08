# -*- coding: utf-8 -*-
import numpy as np
from .BaseDataProtocol.PAProtocol import dtype_PA
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


MAX_PA_BEAMS = 4096
MAX_PA_SWEEPS = 512
MAX_PA_MOMENTS_PER_RADIAL = 64
MAX_PA_MOMENT_DATA_BYTES = 32 * 1024 * 1024


class PABaseData(object):
    """
    解码新一代双偏振的数据格式
    """

    def __init__(self, filename, station_lon=None, station_lat=None, station_alt=None):
        """
        :param filename:  radar basedata filename
        :param station_lon:  radar station longitude //units: degree east
        :param station_lat:  radar station latitude //units:degree north
        :param station_alt:  radar station altitude //units: meters
        """
        super(PABaseData, self).__init__()
        self.filename = filename
        self.station_lon = station_lon
        self.station_lat = station_lat
        self.station_alt = station_alt
        self.fid = _prepare_for_read(self.filename)  ##对压缩的文件进行解码
        self._check_standard_basedata()  ##确定文件是standard文件
        self.header = self._parse_BaseDataHeader()
        self.radial = self._parse_radial()
        self.nrays = len(self.radial)
        self.nsweeps = _validate_count("CutNumber", self.header['TaskConfig']['CutNumber'], minimum=1, maximum=MAX_PA_SWEEPS)
        if self.nrays < self.nsweeps:
            raise ValueError("PA radial count is smaller than the declared sweep count.")
        self.sweep_start_ray_index, self.sweep_end_ray_index = self._build_sweep_indices()
        self.fid.close()

    def _check_standard_basedata(self):
        """
        :param fid: file fid
        :return:
        """
        if _read_exact(self.fid, 12, "PA file signature")[8:] != b'\x10\x00\x00\x00':
            raise ValueError("File is not a standard phase-array base file.")
        self.fid.seek(0, 0)
        return

    def _parse_BaseDataHeader(self):
        BaseDataHeader = {}
        fixed_buf = _read_exact(self.fid, dtype_PA.BeamConfigurationBlockPos, "PA fixed header")  ##读取前面固定头的信息

        BaseDataHeader['GenericHeader'], _ = _unpack_from_buf(fixed_buf, \
                                                              dtype_PA.GenericHeaderBlockPos,
                                                              dtype_PA.BaseDataHeader['GenericHeaderBlock'])
        BaseDataHeader['SiteConfig'], _ = _unpack_from_buf(fixed_buf, \
                                                           dtype_PA.SiteConfigurationBlockPos,
                                                           dtype_PA.BaseDataHeader['SiteConfigurationBlock'])
        BaseDataHeader['TaskConfig'], _ = _unpack_from_buf(fixed_buf,
                                                           dtype_PA.TaskConfigurationBlockPos,
                                                           dtype_PA.BaseDataHeader['TaskConfigurationBlock'])

        beam_count = _validate_count("BeamNumber", BaseDataHeader['TaskConfig']['BeamNumber'], minimum=0, maximum=MAX_PA_BEAMS)
        cut_count = _validate_count("CutNumber", BaseDataHeader['TaskConfig']['CutNumber'], minimum=1, maximum=MAX_PA_SWEEPS)
        beam_buf = _read_exact(
            self.fid,
            dtype_PA.BeamConfigurationBlockSize * beam_count,
            "PA beam configuration block",
        )
        cut_buf = _read_exact(
            self.fid,
            dtype_PA.CutConfigurationBlockSize * cut_count,
            "PA cut configuration block",
        )
        BaseDataHeader["BeamConfig"] = np.frombuffer(beam_buf, dtype_PA.BaseDataHeader['BeamConfigurationBlock'])
        BaseDataHeader['CutConfig'] = np.frombuffer(cut_buf, dtype_PA.BaseDataHeader['CutConfigurationBlock'])
        return BaseDataHeader

    def _parse_radial(self):
        radial = []
        buf = self.fid.read(dtype_PA.RadialHeaderBlockSize)
        while len(buf) == dtype_PA.RadialHeaderBlockSize:  ##read until EOF
            RadialDict, _ = _unpack_from_buf(buf, 0, dtype_PA.RadialHeader())
            self.MomentNum = _validate_count(
                "PA MomentNumber",
                RadialDict['MomentNumber'],
                minimum=0,
                maximum=MAX_PA_MOMENTS_PER_RADIAL,
            )
            self.LengthOfData = RadialDict['LengthOfData']
            RadialDict['fields'] = self._parse_radial_single()
            if RadialDict["fields"]:
                radial.append(RadialDict)
            buf = self.fid.read(dtype_PA.RadialHeaderBlockSize)
        if buf:
            raise ValueError("PA radial header is truncated.")
        return radial

    def _build_sweep_indices(self):
        """Build sweep start/end indices from explicit PA radial metadata."""
        if not self.radial:
            return np.array([], dtype=np.int32), np.array([], dtype=np.int32)

        elevation_starts = [0]
        state_starts = [0]
        for idx in range(1, self.nrays):
            prev = self.radial[idx - 1]
            curr = self.radial[idx]
            if int(curr["ElevationNumber"]) != int(prev["ElevationNumber"]):
                elevation_starts.append(idx)
            if int(curr["RadialNumber"]) == 1 or int(curr["RadialState"]) in (0, 3):
                state_starts.append(idx)

        elevation_starts = np.asarray(sorted(set(elevation_starts)), dtype=np.int32)
        state_starts = np.asarray(sorted(set(state_starts)), dtype=np.int32)

        if state_starts.size == self.nsweeps:
            starts = state_starts
        elif elevation_starts.size == self.nsweeps:
            starts = elevation_starts
        else:
            raise ValueError(
                "PA sweep segmentation does not match the declared CutNumber: "
                "declared %d, elevation-based %d, state-based %d."
                % (self.nsweeps, elevation_starts.size, state_starts.size)
            )

        ends = np.concatenate((starts[1:] - 1, np.array([self.nrays - 1], dtype=np.int32)))
        self._validate_sweep_layout(starts, ends)
        return starts, ends

    def _validate_sweep_layout(self, starts, ends):
        if starts.size != ends.size:
            raise ValueError("PA sweep start/end index arrays are inconsistent.")
        if starts.size != self.nsweeps:
            raise ValueError("PA sweep count does not match CutNumber.")
        if int(starts[0]) != 0 or int(ends[-1]) != self.nrays - 1:
            raise ValueError("PA sweep boundaries do not cover the full radial sequence.")

        for start, end in zip(starts, ends):
            if int(end) < int(start):
                raise ValueError("PA sweep contains a negative-length radial segment.")
            segment = self.radial[int(start): int(end) + 1]
            elevation_numbers = {int(radial["ElevationNumber"]) for radial in segment}
            if len(elevation_numbers) != 1:
                raise ValueError("PA sweep contains multiple elevation numbers.")
            radial_numbers = np.array([int(radial["RadialNumber"]) for radial in segment], dtype=np.int32)
            if radial_numbers[0] != 1:
                raise ValueError("PA sweep does not start at radial number 1.")
            if np.any(np.diff(radial_numbers) <= 0):
                raise ValueError("PA radial numbers are not strictly increasing within a sweep.")

    def _parse_radial_single(self):
        radial_var = {}
        for _ in range(self.MomentNum):
            Mom_buf = _read_exact(self.fid, dtype_PA.MomentHeaderBlockSize, "PA moment header")
            Momheader, _ = _unpack_from_buf(Mom_buf, 0, dtype_PA.RadialData())
            data_len = _validate_count(
                "PA moment length",
                Momheader['Length'],
                minimum=0,
                maximum=MAX_PA_MOMENT_DATA_BYTES,
            )
            Data_buf = _read_exact(self.fid, data_len, "PA moment payload")
            if Momheader['BinLength'] not in (1, 2):
                raise ValueError("PA moment bin length must be 1 or 2 bytes.")
            if data_len % int(Momheader['BinLength']) != 0:
                raise ValueError("PA moment payload length does not match the bin length.")
            if int(Momheader['Scale']) == 0:
                raise ValueError("PA moment scale cannot be zero.")
            if Momheader['BinLength'] == 1:
                dat_tmp = (np.frombuffer(Data_buf, dtype="u1", offset=0)).astype(int)
            else:
                dat_tmp = (np.frombuffer(Data_buf, dtype="u2", offset=0)).astype(int)
            field_name = dtype_PA.flag2Product.get(int(Momheader['DataType']))
            if field_name is None:
                continue
            radial_var[field_name] = np.where(
                dat_tmp >= 5,
                (dat_tmp - Momheader['Offset']) / Momheader['Scale'],
                np.nan,
            ).astype(np.float32)
        return radial_var

    def get_nyquist_velocity(self):
        """get nyquist vel per ray
        获取每根径向的不模糊速度
        :return:(nRays)
        """
        return np.concatenate([[nyquist, ] * ray for nyquist, ray in zip(self.header['CutConfig']['NyquistSpeed'], \
                                                                         self.get_rays_per_sweep())], axis=0)

    def get_unambiguous_range(self):
        """
        获取每根径向的不模糊距离
        :return:(nRays)
        """
        return np.concatenate([[nyquist, ] * ray for nyquist, ray in zip(self.header['CutConfig']['MaximumRange'], \
                                                                         self.get_rays_per_sweep())], axis=0)

    def get_scan_time(self):
        """
        获取每根径向的扫描时间
        :return:(nRays)
        """
        return np.array(
            [julian2date_SEC(int(radial["Seconds"]), int(radial["MicroSeconds"])) for radial in self.radial]
        )

    def get_sweep_end_ray_index(self):
        """
        获取每个sweep的结束的index，包含在内
        :return:(nsweep)
        """
        return self.sweep_end_ray_index

    def get_sweep_start_ray_index(self):
        """
        获取每个sweep的开始的index
        :return:(nsweep)
        """
        return self.sweep_start_ray_index

    def get_rays_per_sweep(self):
        """
        获取每个sweep的径向数
        :return:(nsweep)
        """
        return self.sweep_end_ray_index - self.sweep_start_ray_index + 1

    def get_azimuth(self):
        """
        获取每根径向的方位角
        :return:(nRays)
        """
        return np.array([self.radial[iray]['Azimuth'] for iray in range(self.nrays)])

    def get_elevation(self):
        """
        获取每根径向的仰角
        :return: (nRays)
        """
        return np.array([self.radial[iray]['Elevation'] for iray in range(self.nrays)])

    def get_latitude_longitude_altitude_frequency(self):
        """
        获取经纬度高度，雷达频率
        :return:lat, lon, alt, frequency
        """
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


class PA2NRadar(object):
    """到NRadar object 的桥梁"""

    def __init__(self, pa_reader):
        super(PA2NRadar, self).__init__()
        self.pa_reader = pa_reader
        self.radial = list(self.pa_reader.radial)
        self.nrays = len(self.radial)
        self.nsweeps = int(self.pa_reader.header['TaskConfig']['CutNumber'])
        self.scan_type = self.pa_reader.get_scan_type()
        self.latitude, self.longitude, self.altitude, self.frequency = \
            self.pa_reader.get_latitude_longitude_altitude_frequency()
        self.header = self.pa_reader.header
        self.sweep_start_ray_index = np.asarray(self.pa_reader.sweep_start_ray_index, dtype=np.int32)
        self.sweep_end_ray_index = np.asarray(self.pa_reader.sweep_end_ray_index, dtype=np.int32)
        self.rays_per_sweep = self.get_rays_per_sweep()
        self.bins_per_sweep = self.get_nbins_per_sweep()
        self.range = self.get_range_per_radial(self.bins_per_sweep.max())
        self.azimuth = self.get_azimuth()
        self.elevation = self.get_elevation()
        self.extended_fields = self._build_extended_fields()
        self.fields = self._get_fields()
        self.sitename = self.pa_reader.get_sitename()


    def get_nbins_per_sweep(self):
        """
        确定每个 sweep 的 aligned 距离库长度。
        :return:
        """
        return np.array(
            [
                self._aligned_bins_for_sweep(int(start), int(end))
                for start, end in zip(self.sweep_start_ray_index, self.sweep_end_ray_index)
            ],
            dtype=np.int32,
        )

    def _aligned_bins_for_sweep(self, start, end):
        preferred_lengths = []
        reflectivity_lengths = []
        for iray in range(start, end + 1):
            for field_name, values in self.radial[iray]["fields"].items():
                if field_name in ("dBT", "dBZ", "Zc"):
                    reflectivity_lengths.append(values.size)
                else:
                    preferred_lengths.append(values.size)
        if preferred_lengths:
            return int(max(preferred_lengths))
        if reflectivity_lengths:
            return int(max(reflectivity_lengths))
        return 0

    def get_rays_per_sweep(self):
        return self.sweep_end_ray_index - self.sweep_start_ray_index + 1

    def get_azimuth(self):
        """
        获取每根径向的方位角
        :return:(nRays)
        """
        return np.array([self.radial[iray]['Azimuth'] for iray in range(self.nrays)])

    def get_elevation(self):
        """
        获取每根径向的仰角
        :return: (nRays)
        """
        elevation = np.array([self.radial[iray]['Elevation'] for iray in range(self.nrays)])
        return np.where(elevation>180, elevation-360, elevation)

    def get_scan_time(self):
        """
        获取每根径向的扫描时间
        :return:(nRays)
        """
        return np.array(
            [julian2date_SEC(int(radial["Seconds"]), int(radial["MicroSeconds"])) for radial in self.radial]
        )

    def get_nyquist_velocity(self):
        """get nyquist vel per ray
        获取每根径向的不模糊速度
        :return:(nRays)
        """
        return np.concatenate(
            [[nyquist] * int(ray_count) for nyquist, ray_count in zip(self.header['CutConfig']['NyquistSpeed'], self.rays_per_sweep)],
            axis=0,
        )

    def get_unambiguous_range(self):
        """
        获取每根径向的不模糊距离
        :return:(nRays)
        """
        return np.concatenate(
            [[nyquist] * int(ray_count) for nyquist, ray_count in zip(self.header['CutConfig']['MaximumRange'], self.rays_per_sweep)],
            axis=0,
        )

    def get_range_per_radial(self, length):
        """
        确定径向每个库的距离 range变量
        :param length:
        :return:
        """
        Resolution = self.header['CutConfig']['DopplerResolution'][0]
        if Resolution == 0:
            return np.linspace(30, 30 * length, length)
        else:
            return np.linspace(Resolution, Resolution * length, length)

    def get_dbz_range_per_radial(self, length):
        """
        确定径向每个库的距离
        :param length:
        :return:
        """
        Resolution = self.header['CutConfig']['LogResolution'][0]
        if Resolution == 0:
            return np.linspace(30, 30 * length, length)
        else:
            return np.linspace(Resolution, Resolution * length, length)

    def _get_fields(self):
        """将所有的field的数据提取出来"""
        fields = {}
        field_keys = sorted({key for radial in self.radial for key in radial['fields'].keys()})
        for ikey in field_keys:
            fields[ikey] = np.array([self._add_or_del_field(iray['fields'], ikey) for iray in self.radial])
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
                "range": self.get_dbz_range_per_radial(native_bins).astype(np.float32, copy=False),
                "time": np.asarray(scan_time[start:end + 1]),
                "azimuth": np.asarray(self.azimuth[start:end + 1], dtype=np.float32),
                "elevation": np.asarray(self.elevation[start:end + 1], dtype=np.float32),
                "aligned_bins": int(aligned_bins),
            }

        return {"dBZ": extended_sweeps} if extended_sweeps else {}

    def _add_or_del_field(self, dat_fields, key):
        """
        根据fields的key提取数据
        :param dat_fields: fields的数据
        :param key: key words
        :param flag_match: dop和dbz分辨率是否匹配, 匹配则为True，不匹配为False
        :return:
        """
        length = self.bins_per_sweep.max()
        if key not in dat_fields.keys():
            return np.full((length,), np.nan)

        dat_ray = dat_fields[key]
        assert dat_ray.ndim == 1, "check dat_ray"
        if dat_ray.size >= length:
            return dat_ray[:length]
        else:
            out = np.full((length,), np.nan)
            out[:dat_ray.size] = dat_ray
            return out

    def get_NRadar_nyquist_speed(self):
        """array shape (nsweeps)"""
        return self.header['CutConfig']['NyquistSpeed']

    def get_NRadar_unambiguous_range(self):
        """array shape (nsweeps)"""
        return self.header['CutConfig']['MaximumRange']

    def get_fixed_angle(self):
        if self.scan_type == "rhi":
            return self.header['CutConfig']['Azimuth']
        else:
            return self.header['CutConfig']['Elevation']

    def ToPRD(self, effective_earth_radius=None):
        """将WSR98D数据转为PRD的数据格式"""
        return PRD(fields=self.fields, scan_type=self.scan_type, time=self.get_scan_time(), \
                          range=self.range, azimuth=self.azimuth, elevation=self.elevation, latitude=self.latitude, \
                          longitude=self.longitude, altitude=self.altitude,
                          sweep_start_ray_index=self.sweep_start_ray_index, \
                          sweep_end_ray_index=self.sweep_end_ray_index, fixed_angle=self.get_fixed_angle(), \
                          bins_per_sweep=self.bins_per_sweep, nyquist_velocity=self.get_NRadar_nyquist_speed(), \
                          frequency=self.frequency, unambiguous_range=self.get_NRadar_unambiguous_range(), \
                          nrays=self.nrays, nsweeps=self.nsweeps, sitename = self.sitename,
                          pyart_radar=None, effective_earth_radius=effective_earth_radius,
                          extended_fields=self.extended_fields,
                          metadata={"original_container": "PA", "radar_name": "PA"})

    def to_pyart_radar(self, effective_earth_radius=None, **kwargs):
        """Export the decoded PA volume through the PRD Py-ART adapter."""
        return self.ToPRD(effective_earth_radius=effective_earth_radius).to_pyart_radar(**kwargs)

    def ToPyartRadar(self, effective_earth_radius=None, **kwargs):
        """Backward-compatible alias for ``to_pyart_radar``."""
        return self.to_pyart_radar(effective_earth_radius=effective_earth_radius, **kwargs)

    def _get_instrument_parameters(self):
        """ Return a dictionary containing instrument parameters. """

        # pulse width
        pulse_width = get_metadata('pulse_width')
        pulse_width['data'] = np.array([self.header["BeamConfig"]['SubPulseBandWidth']/10**9,], dtype='float32')  # nanosec->sec
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
