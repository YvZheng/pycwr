# -*- coding: utf-8 -*-
import numpy as np
from .BaseDataProtocol.WSR98DProtocol import dtype_98D
from .util import _prepare_for_read, _unpack_from_buf, julian2date_SEC, make_time_unit_str
from ..core.NRadar import PRD
from ..configure.pyart_config import get_metadata, get_fillvalue
from ..configure.default_config import CINRAD_field_mapping
from ..core.PyartRadar import Radar
from netCDF4 import date2num
from scipy import interpolate

class WSR98DBaseData(object):
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
        super(WSR98DBaseData, self).__init__()
        self.filename = filename
        self.station_lon = station_lon
        self.station_lat = station_lat
        self.station_alt = station_alt
        self.fid = _prepare_for_read(self.filename)  ##对压缩的文件进行解码
        self._check_standard_basedata()  ##确定文件是standard文件
        self.header = self._parse_BaseDataHeader()
        self.radial = self._parse_radial()
        self.nrays = len(self.radial)
        status = np.array([istatus['RadialState'] for istatus in self.radial[:]])
        self.sweep_start_ray_index = np.where((status == 0) | (status == 3))[0]
        self.sweep_end_ray_index = np.where((status == 2) | (status == 4))[0]
        self.nsweeps = len(self.sweep_start_ray_index)
        self.fid.close()

    def _check_standard_basedata(self):
        """
        :param fid: file fid
        :return:
        """
        assert self.fid.read(4) == b'RSTM', 'file in not a stardand WSR-98D file!'
        self.fid.seek(0, 0)
        return

    def _parse_BaseDataHeader(self):
        BaseDataHeader = {}
        fixed_buf = self.fid.read(dtype_98D.CutConfigurationBlockPos)  ##读取前面固定头的信息

        BaseDataHeader['GenericHeader'], _ = _unpack_from_buf(fixed_buf, \
                                                              dtype_98D.GenericHeaderBlockPos,
                                                              dtype_98D.BaseDataHeader['GenericHeaderBlock'])
        BaseDataHeader['SiteConfig'], _ = _unpack_from_buf(fixed_buf, \
                                                           dtype_98D.SiteConfigurationBlockPos,
                                                           dtype_98D.BaseDataHeader['SiteConfigurationBlock'])
        BaseDataHeader['TaskConfig'], _ = _unpack_from_buf(fixed_buf,
                                                           dtype_98D.TaskConfigurationBlockPos,
                                                           dtype_98D.BaseDataHeader['TaskConfigurationBlock'])
        cut_buf = self.fid.read(dtype_98D.CutConfigurationBlockSize * \
                                BaseDataHeader['TaskConfig']['CutNumber'])
        BaseDataHeader['CutConfig'] = np.frombuffer(cut_buf, dtype_98D.BaseDataHeader['CutConfigurationBlock'])
        return BaseDataHeader

    def _parse_radial(self):
        radial = []
        buf = self.fid.read(dtype_98D.RadialHeaderBlockSize)
        while len(buf) == dtype_98D.RadialHeaderBlockSize:  ##read until EOF
            RadialDict, _ = _unpack_from_buf(buf, 0, dtype_98D.RadialHeader())
            self.MomentNum = RadialDict['MomentNumber']
            self.LengthOfData = RadialDict['LengthOfData']
            RadialDict['fields'] = self._parse_radial_single()
            radial.append(RadialDict)
            buf = self.fid.read(dtype_98D.RadialHeaderBlockSize)
        return radial

    def _parse_radial_single(self):
        radial_var = {}
        for _ in range(self.MomentNum):
            Mom_buf = self.fid.read(dtype_98D.MomentHeaderBlockSize)
            Momheader, _ = _unpack_from_buf(Mom_buf, 0, dtype_98D.RadialData())
            Data_buf = self.fid.read(Momheader['Length'])
            assert (Momheader['BinLength'] == 1) | (Momheader['BinLength'] == 2), "Bin Length has problem!"
            if Momheader['BinLength'] == 1:
                dat_tmp = (np.frombuffer(Data_buf, dtype="u1", offset=0)).astype(int)
            else:
                dat_tmp = (np.frombuffer(Data_buf, dtype="u2", offset=0)).astype(int)
            radial_var[dtype_98D.flag2Product[Momheader['DataType']]] = np.where(dat_tmp >= 5, \
                                                                                 (dat_tmp - Momheader['Offset']) /
                                                                                 Momheader['Scale'],
                                                                                 np.nan).astype(np.float32)
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
        return np.array([julian2date_SEC(iray['Seconds'], iray['MicroSeconds']) for iray in self.radial])

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


class WSR98D2NRadar(object):
    """到NusitRadar object 的桥梁"""

    def __init__(self, WSR98D):
        self.WSR98D = WSR98D
        self.flag_match = np.all(self.WSR98D.header['CutConfig']['LogResolution'] == \
                       self.WSR98D.header['CutConfig']['DopplerResolution'])
        self.v_index_alone = self.get_v_idx()
        self.dBZ_index_alone = self.get_dbz_idx()

        if self.WSR98D.header["TaskConfig"]["TaskName"][:6] in [b'VCP26D', b'VCP27D']: #HZ Radar problem!
            self.interp_VCP26(self.dBZ_index_alone, self.v_index_alone) ##process V and dBZ are not match
            self.dBZ_index_alone = self.dBZ_index_alone[:min(len(self.dBZ_index_alone), len(self.v_index_alone))] ##要去除的sweep
        else:
            assert np.all(self.v_index_alone == (self.dBZ_index_alone + 1)), """v and dBZ not equal!"""
            for index_with_dbz, index_with_v in zip(self.dBZ_index_alone, self.v_index_alone):
                assert (self.WSR98D.header["CutConfig"]["Elevation"][index_with_v] - \
                        self.WSR98D.header["CutConfig"]["Elevation"][index_with_dbz]) < 0.5, \
                    "warning! maybe it is a problem."
                self.interp_dBZ(index_with_dbz, index_with_v)

        ind_remove = self.get_reomve_radial_num()
        self.radial = [iray for ind, iray in enumerate(self.WSR98D.radial) if ind not in ind_remove]

        status = np.array([istatus['RadialState'] for istatus in self.radial[:]])
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
        self.fields = self._get_fields()
        self.sitename = self.WSR98D.get_sitename()

    def get_reomve_radial_num(self):
        """获得需要remove的radial的index"""
        dBZ_alone = self.dBZ_index_alone
        index_romove = []
        for isweep in dBZ_alone:
            index_romove.extend(range(self.WSR98D.sweep_start_ray_index[isweep], \
                                      self.WSR98D.sweep_end_ray_index[isweep] + 1))
        return index_romove

    def get_v_idx(self):
        """获取需要插值的sweep, 插值到有径向速度仰角"""
        flag = np.array([(("V" in self.WSR98D.radial[idx]['fields'].keys()) and (
                    "dBZ" not in self.WSR98D.radial[idx]['fields'].keys())) \
                         for idx in self.WSR98D.sweep_start_ray_index])
        return np.where(flag == 1)[0]

    def get_dbz_idx(self):
        """获取含有dbz的sweep"""
        flag = np.array([(("dBZ" in self.WSR98D.radial[idx]['fields'].keys()) and (
                    "V" not in self.WSR98D.radial[idx]['fields'].keys())) \
                         for idx in self.WSR98D.sweep_start_ray_index])
        return np.where(flag == 1)[0]

    def interp_VCP26(self, dBZ_sweep_index, V_sweep_index):
        """
        处理VCP26不匹配问题
        :param dBZ_sweep_index: array, dbz单独扫描的仰角
        :param V_sweep_index: array, dBZ单独扫描的仰角
        :return:
        """
        add_keys = ["V", "W"]
        same_sweeps = min(len(dBZ_sweep_index), len(V_sweep_index))
        for isweep in range(same_sweeps):
            self.interp_dBZ(dBZ_sweep_index[isweep], V_sweep_index[isweep])
        for dbz_dense in dBZ_sweep_index[same_sweeps:]:
            for index in range(self.WSR98D.sweep_start_ray_index[dbz_dense], self.WSR98D.sweep_end_ray_index[dbz_dense]+1):
                for ikey in add_keys:
                    self.WSR98D.radial[index]["fields"][ikey] = np.full_like(self.WSR98D.radial[index]["fields"]["dBZ"],\
                                                                             np.nan, dtype=np.float32)

    def interp_dBZ(self, field_with_dBZ_num, field_without_dBZ_num):
        """
        将dBZ插值到不含dBZ的仰角
        :param field_with_dBZ_num: 要插值的sweep num, （从0开始）
        :param field_without_dBZ_num: 要插值到的sweep num, (从0开始)  which to evaluate the interpolated values
        :return:
        """
        azimuth = self.WSR98D.get_azimuth() ##获取98D的方位角
        assert (field_with_dBZ_num + 1) == field_without_dBZ_num, "check interp sweep!"
        dbz_az = azimuth[self.WSR98D.sweep_start_ray_index[field_with_dBZ_num]: \
                         self.WSR98D.sweep_end_ray_index[field_with_dBZ_num] + 1]
        v_az = azimuth[self.WSR98D.sweep_start_ray_index[field_without_dBZ_num]: \
                       self.WSR98D.sweep_end_ray_index[field_without_dBZ_num] + 1]
        dbz_idx = np.argmin(np.abs(dbz_az.reshape(-1, 1) - v_az.reshape(1, -1)), axis=0) + \
                  self.WSR98D.sweep_start_ray_index[field_with_dBZ_num] ##最邻近插值

        v_idx = np.arange(self.WSR98D.sweep_start_ray_index[field_without_dBZ_num], \
                          self.WSR98D.sweep_end_ray_index[field_without_dBZ_num] + 1)
        keys = self.WSR98D.radial[self.WSR98D.sweep_start_ray_index[field_with_dBZ_num]]['fields'].keys()
        for ind_dbz, ind_v in zip(dbz_idx, v_idx):
            for ikey in keys:
                self.WSR98D.radial[ind_v]["fields"][ikey] = self.WSR98D.radial[ind_dbz]["fields"][ikey]

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

    def get_rays_per_sweep(self):
        """
        获取每个sweep的径向数
        :return:(nsweep)
        """
        return self.sweep_end_ray_index - self.sweep_start_ray_index + 1

    def get_scan_time(self):
        """
        获取每根径向的扫描时间
        :return:(nRays)
        """
        return np.array([julian2date_SEC(iray['Seconds'], iray['MicroSeconds']) for iray in self.radial])

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

    def get_nbins_per_sweep(self):
        """
        确定每个sweep V探测的库数
        :return:
        """
        return np.array([self.radial[idx]['fields']['V'].size for idx in self.sweep_start_ray_index])

    def get_range_per_radial(self, length):
        """
        确定径向每个库的距离 range变量
        :param length:
        :return:
        """
        Resolution = self.header['CutConfig']['DopplerResolution'][0]
        return np.linspace(Resolution, Resolution * length, length)

    def get_dbz_range_per_radial(self, length):
        """
        确定径向每个库的距离
        :param length:
        :return:
        """
        Resolution = self.header['CutConfig']['LogResolution'][0]
        return np.linspace(Resolution, Resolution * length, length)

    def _get_fields(self):
        """将所有的field的数据提取出来"""
        fields = {}
        field_keys = self.radial[0]['fields'].keys()
        for ikey in field_keys:
            fields[ikey] = np.array([self._add_or_del_field(iray['fields'], ikey, self.flag_match) for iray in self.radial])
        return fields

    def _add_or_del_field(self, dat_fields, key, flag_match=True):
        """
        根据fields的key提取数据
        :param dat_fields: fields的数据
        :param key: key words
        :param flag_match: dop和dbz分辨率是否匹配, 匹配则为True，不匹配为False
        :return:
        """
        length = self.max_bins
        if key not in dat_fields.keys():
            return np.full((length,), np.nan)

        if flag_match == False:
            if key == "dBZ":
                dbz_range = self.get_dbz_range_per_radial(dat_fields[key].size)
                dop_range = self.get_range_per_radial(length)
                match_data = interpolate.interp1d(dbz_range, dat_fields[key], kind="nearest",
                                                  bounds_error=False, fill_value=np.nan)
                dat_ray = match_data(dop_range)
                return dat_ray.ravel()

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

    def ToPRD(self):
        """将WSR98D数据转为PRD的数据格式"""
        return PRD(fields=self.fields, scan_type=self.scan_type, time=self.get_scan_time(), \
                          range=self.range, azimuth=self.azimuth, elevation=self.elevation, latitude=self.latitude, \
                          longitude=self.longitude, altitude=self.altitude,
                          sweep_start_ray_index=self.sweep_start_ray_index, \
                          sweep_end_ray_index=self.sweep_end_ray_index, fixed_angle=self.get_fixed_angle(), \
                          bins_per_sweep=self.bins_per_sweep, nyquist_velocity=self.get_NRadar_nyquist_speed(), \
                          frequency=self.frequency, unambiguous_range=self.get_NRadar_unambiguous_range(), \
                          nrays=self.nrays, nsweeps=self.nsweeps, sitename = self.sitename, pyart_radar=self.ToPyartRadar())

    def ToPyartRadar(self):
        """转化为Pyart Radar的对象"""
        dts = self.get_scan_time()
        units = make_time_unit_str(min(dts))
        time = get_metadata('time')
        time['units'] = units
        time['data'] = date2num(dts, units).astype('float32')

        # range
        _range = get_metadata('range')
        # assume that the number of gates and spacing from the first ray is
        # representative of the entire volume
        _range['data'] = self.range
        _range['meters_to_center_of_first_gate'] = self.header['CutConfig']['DopplerResolution'][0]
        _range['meters_between_gates'] = self.header['CutConfig']['DopplerResolution'][0]

        latitude = get_metadata('latitude')
        longitude = get_metadata('longitude')
        altitude = get_metadata('altitude')
        latitude['data'] = np.array([self.latitude], dtype='float64')
        longitude['data'] = np.array([self.longitude], dtype='float64')
        altitude['data'] = np.array([self.altitude], dtype='float64')

        metadata = get_metadata('metadata')
        metadata['original_container'] = 'WSR98D'
        metadata['site_name'] = self.sitename
        metadata['radar_name'] = "WSR98D"

        sweep_start_ray_index = get_metadata('sweep_start_ray_index')
        sweep_end_ray_index = get_metadata('sweep_end_ray_index')
        sweep_start_ray_index['data'] = self.sweep_start_ray_index
        sweep_end_ray_index['data'] = self.sweep_end_ray_index

        sweep_number = get_metadata('sweep_number')
        sweep_number['data'] = np.arange(self.nsweeps, dtype='int32')

        scan_type = self.scan_type

        sweep_mode = get_metadata('sweep_mode')
        if self.scan_type == "ppi":
            sweep_mode['data'] = np.array(self.nsweeps * ['azimuth_surveillance'], dtype='S')
        elif self.scan_type == "rhi":
            sweep_mode['data'] = np.array(self.nsweeps * ['rhi'], dtype='S')
        else:
            sweep_mode['data'] = np.array(self.nsweeps * ['sector'], dtype='S')

        # elevation
        elevation = get_metadata('elevation')
        elevation['data'] = self.elevation

        # azimuth
        azimuth = get_metadata('azimuth')
        azimuth['data'] = self.azimuth

        # fixed_angle
        fixed_angle = get_metadata('fixed_angle')
        fixed_angle['data'] = self.get_fixed_angle()

        # instrument_parameters
        instrument_parameters = self._get_instrument_parameters()

        # fields
        fields = {}
        for field_name_abbr in self.fields.keys():
            field_name = CINRAD_field_mapping[field_name_abbr]
            if field_name is None:
                continue
            field_dic = get_metadata(field_name)
            field_dic['data'] = np.ma.masked_array(self.fields[field_name_abbr],\
                                mask=np.isnan(self.fields[field_name_abbr]), fill_value=get_fillvalue())
            field_dic['_FillValue'] = get_fillvalue()
            fields[field_name] = field_dic
        return Radar(time, _range, fields, metadata, scan_type,
                     latitude, longitude, altitude,
                     sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
                     sweep_end_ray_index,
                     azimuth, elevation,
                     instrument_parameters=instrument_parameters)

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
