# -*- coding: utf-8 -*-
import numpy as np
from .BaseDataProtocol.CCProtocol import dtype_cc
from .util import _prepare_for_read, _unpack_from_buf, make_time_unit_str, get_radar_sitename
import datetime
import pandas as pd
from ..core.NRadar import PRD
from ..configure.pyart_config import get_metadata, get_fillvalue
from ..configure.default_config import CINRAD_field_mapping, _LIGHT_SPEED
from ..core.PyartRadar import Radar
from netCDF4 import date2num

class CCBaseData(object):
    """
        解码CC/CCJ的数据格式
    """

    def __init__(self, filename, station_lon=None, station_lat=None, station_alt=None):
        """
                :param filename:  radar basedata filename
                :param station_lon:  radar station longitude //units: degree east
                :param station_lat:  radar station latitude //units:degree north
                :param station_alt:  radar station altitude //units: meters
        """
        super(CCBaseData, self).__init__()
        self.filename = filename
        self.station_lon = station_lon
        self.station_lat = station_lat
        self.station_alt = station_alt
        self.fid = _prepare_for_read(self.filename)  ##判断是否是压缩文件
        # print(len(self.fid.read()))
        buf_header = self.fid.read(dtype_cc.BaseDataHeaderSize)  ##取出header的buf
        self.header = self._parse_BaseDataHeader(buf_header)
        self._check_cc_basedata()
        self.fid.seek(dtype_cc.BaseDataHeaderSize, 0)  ##移动到径向数据的位置
        self.radial = self._parse_radial()

    def _check_cc_basedata(self):
        """检查雷达数据是否完整"""
        buf_radial_data = self.fid.read()
        assert len(buf_radial_data) == self.nrays * dtype_cc.PerRadialSize, "CC basedata size has problems!"
        return

    def _parse_BaseDataHeader(self, buf_header):
        """
        :param buf_header: 只包含头文件的buf
        :return:
        """
        BaseDataHeader_dict = {}
        ##解码第一部分观测参数
        BaseDataHeader_dict['ObsParam1'], _ = _unpack_from_buf(buf_header, \
                                                               dtype_cc.HeaderSize1_pos,
                                                               dtype_cc.BaseDataHeader['RadarHeader1'])
        ##解码不同仰角的观测参数
        assert BaseDataHeader_dict['ObsParam1']['ucScanMode'] > 100, "only vol support!"
        self.nsweeps = BaseDataHeader_dict['ObsParam1']['ucScanMode'] - 100
        BaseDataHeader_dict['CutConfig'] = np.frombuffer(buf_header, \
                                                         dtype_cc.BaseDataHeader['CutConfigX30'], count=self.nsweeps,
                                                         offset=dtype_cc.CutSize_pos)
        ##解码第二部分观测参数
        BaseDataHeader_dict['ObsParam2'], _ = _unpack_from_buf(buf_header, \
                                                               dtype_cc.HeaderSize2_pos,
                                                               dtype_cc.BaseDataHeader['RadarHeader2'])

        self.nrays = np.sum(BaseDataHeader_dict['CutConfig']['usRecordNumber'])
        self.sweep_end_ray_index_add1 = (np.cumsum(BaseDataHeader_dict['CutConfig']['usRecordNumber'])).astype(int)  ##python格式的结束
        self.sweep_start_ray_index = (self.sweep_end_ray_index_add1 - BaseDataHeader_dict['CutConfig']['usRecordNumber']).astype(int)
        return BaseDataHeader_dict

    def _parse_radial_single(self, buf_radial, radialnumber):
        """解析径向的数据"""
        Radial = {}
        RadialData = np.frombuffer(buf_radial, dtype_cc.RadialData(radialnumber))
        Radial['fields'] = {}
        Radial['fields']['dBZ'] = np.where(RadialData['dBZ'] != -32768, RadialData['dBZ'] / 10.,
                                           np.nan).astype(np.float32)
        Radial['fields']['V'] = np.where(RadialData['V'] != -32768, RadialData['V'] / 10.,
                                         np.nan).astype(np.float32)
        Radial['fields']['W'] = np.where(RadialData['W'] != -32768, RadialData['W'] / 10.,
                                         np.nan).astype(np.float32)
        return Radial

    def _parse_radial(self):
        radial = []
        for isweep in range(self.nsweeps):
            radialnumber = self.header['CutConfig']['usBinNumber'][isweep]
            for _ in range(self.header['CutConfig']['usRecordNumber'][isweep]):
                buf_radial = self.fid.read(dtype_cc.PerRadialSize)
                radial.append(self._parse_radial_single(buf_radial, radialnumber))
        return radial

    def get_nyquist_velocity(self):
        """get nyquist vel per ray
        获取每根径向的不模糊速度
        :return:(nRays)
        """
        nyquist_velocity = np.concatenate([np.array([self.header['CutConfig']['usMaxV'][isweep] / 100.] \
                                                    * self.header['CutConfig']['usRecordNumber'][isweep]) for \
                                           isweep in range(self.nsweeps)])
        return nyquist_velocity.astype(np.float32)

    def get_unambiguous_range(self):
        """
        获取每根径向的不模糊距离
        :return:(nRays)
        """
        return np.concatenate([np.array([self.header['CutConfig']['usMaxL'][isweep] * 10. \
                                         ] * self.header['CutConfig']['usRecordNumber'][isweep]) for \
                               isweep in range(self.nsweeps)])

    def get_scan_time(self):
        """
        获取每根径向的扫描时间
        :return:(nRays)
        """
        params = self.header['ObsParam1']
        start_year = params['ucSYear1'] * 100 + params['ucSYear2']
        end_year = params['ucEYear1'] * 100 + params['ucEYear2']
        start_time = datetime.datetime(year=start_year, month=params['ucSMonth'],
                                       day=params['ucSDay'], hour=params['ucSHour'],
                                       minute=params['ucSMinute'], second=params['ucSSecond'])
        end_time = datetime.datetime(year=end_year, month=params['ucEMonth'],
                                     day=params['ucEDay'], hour=params['ucEHour'],
                                     minute=params['ucEMinute'], second=params['ucESecond'])
        return pd.date_range(start_time, end_time, periods=self.nrays).to_pydatetime()

    def get_sweep_end_ray_index(self):
        """
        获取每个sweep的结束的index，包含在内
        :return:(nsweep)
        """
        return self.sweep_end_ray_index_add1 - 1

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
        return (self.header['CutConfig']['usRecordNumber']).astype(np.int32)

    def get_azimuth(self):
        """
        获取每根径向的方位角
        :return:(nRays)
        """
        return np.concatenate([np.linspace(0, 360, self.header['CutConfig']['usRecordNumber'][isweep]) \
                               for isweep in range(self.nsweeps)], axis=0)

    def get_elevation(self):
        """
        获取每根径向的仰角
        :return: (nRays)
        """
        return np.concatenate([np.array([self.header['CutConfig']['usAngle'][isweep] / 100. \
                                         ] * self.header['CutConfig']['usRecordNumber'][isweep]) for \
                               isweep in range(self.nsweeps)])

    def get_latitude_longitude_altitude_frequency(self):
        """
        获取经纬度高度，雷达频率
        :return:lat, lon, alt, frequency
        """
        lat, lon, alt, frequency =  self.header['ObsParam1']['lLatitudeValue'] / 3600000.,\
                                    self.header['ObsParam1']['lLongitudeValue'] / 3600000., \
                                    self.header['ObsParam1']['lHeight'] / 1000.,\
                                    3 * 10 ** 5 / self.header['ObsParam2']['lWavelength']
        if self.station_lon is not None:
            lon = self.station_lon
        if self.station_lat is not None:
            lat = self.station_lat
        if self.station_alt is not None:
            alt = self.station_alt
        return lat, lon, alt, frequency

    def get_scan_type(self):
        """
        获取扫描的类型
        :return:
        """
        ## only ppi support!
        return "ppi"

    def get_sitename(self):
        return get_radar_sitename(self.filename)

class CC2NRadar(object):
    """到NusitRadar object 的桥梁"""

    def __init__(self, CC):
        self.CC = CC
        self.radial = self.CC.radial
        self.azimuth = self.get_azimuth()
        self.elevation = self.get_elevation()
        self.sweep_start_ray_index = self.get_sweep_start_ray_index()
        self.sweep_end_ray_index = self.get_sweep_end_ray_index()
        self.nrays = self.CC.nrays
        self.nsweeps = self.CC.nsweeps
        self.scan_type = self.CC.get_scan_type()
        self.latitude, self.longitude, self.altitude, self.frequency = \
            self.CC.get_latitude_longitude_altitude_frequency()
        self.bins_per_sweep = self.get_nbins_per_sweep()
        self.max_bins = self.bins_per_sweep.max()
        self.range = self.get_range_per_radial(self.max_bins)
        self.fields = self._get_fields()
        self.sitename = self.CC.get_sitename()

    def get_azimuth(self):
        """
        获取每根径向的方位角
        :return:(nRays)
        """
        return self.CC.get_azimuth()

    def get_elevation(self):
        """
        获取每根径向的仰角
        :return: (nRays)
        """
        return self.CC.get_elevation()

    def get_rays_per_sweep(self):
        """
        获取每个sweep的径向数
        :return:(nsweep)
        """
        return (self.CC.header['CutConfig']['usRecordNumber']).astype(int)

    def get_scan_time(self):
        """
        获取每根径向的扫描时间
        :return:(nRays)
        """
        return self.CC.get_scan_time()

    def get_nyquist_velocity(self):
        """get nyquist vel per ray
        获取每根径向的不模糊速度
        :return:(nRays)
        """
        return self.CC.get_nyquist_velocity()

    def get_unambiguous_range(self):
        """
        获取每根径向的不模糊距离
        :return:(nRays)
        """
        return self.CC.get_unambiguous_range()

    def get_sweep_end_ray_index(self):
        """
        获取每个sweep的结束的index，包含在内
        :return:(nsweep)
        """
        return self.CC.sweep_end_ray_index_add1 - 1

    def get_sweep_start_ray_index(self):
        """
        获取每个sweep的开始的index
        :return:(nsweep)
        """
        return self.CC.sweep_start_ray_index

    def get_nbins_per_sweep(self):
        """
        确定每个sweep V探测的库数
        :return:
        """
        return (self.CC.header['CutConfig']['usBinNumber']).astype(int)

    def get_range_per_radial(self, length):
        """
        确定径向每个库的距离
        :param length:
        :return:
        """
        Resolution = self.CC.header['CutConfig']['usBindWidth'][0] * 2
        return np.linspace(Resolution, Resolution * length, length)

    def _get_fields(self):
        """将所有的field的数据提取出来"""
        fields = {}
        field_keys = self.radial[0]['fields'].keys()
        for ikey in field_keys:
            fields[ikey] = np.array([(iray['fields'][ikey]).ravel() for iray in self.radial])
        return fields

    def get_NRadar_nyquist_speed(self):
        """array shape (nsweeps)"""
        return self.CC.header['CutConfig']['usMaxV'] / 100.

    def get_NRadar_unambiguous_range(self):
        """array shape (nsweeps)"""
        return self.CC.header['CutConfig']['usMaxL'] * 10.

    def get_fixed_angle(self):
        return self.CC.header['CutConfig']['usAngle'] / 100.

    def ToPRD(self):
        """将WSR98D数据转为PRD 的数据格式"""

        return PRD(fields=self.fields, scan_type=self.scan_type, time=self.get_scan_time(), \
                          range=self.range, azimuth=self.azimuth, elevation=self.elevation, latitude=self.latitude, \
                          longitude=self.longitude, altitude=self.altitude,
                          sweep_start_ray_index=self.sweep_start_ray_index, \
                          sweep_end_ray_index=self.sweep_end_ray_index, fixed_angle=self.get_fixed_angle(), \
                          bins_per_sweep=self.bins_per_sweep, nyquist_velocity=self.get_NRadar_nyquist_speed(), \
                          frequency=self.frequency, unambiguous_range=self.get_NRadar_unambiguous_range(), \
                          nrays=self.nrays, nsweeps=self.nsweeps, sitename = self.sitename, pyart_radar=self.ToPyartRadar())

    def ToPyartRadar(self):

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
        _range['meters_to_center_of_first_gate'] = self.CC.header['CutConfig']['usBindWidth'][0] * 2
        _range['meters_between_gates'] = self.CC.header['CutConfig']['usBindWidth'][0] * 2

        latitude = get_metadata('latitude')
        longitude = get_metadata('longitude')
        altitude = get_metadata('altitude')
        latitude['data'] = np.array([self.latitude], dtype='float64')
        longitude['data'] = np.array([self.longitude], dtype='float64')
        altitude['data'] = np.array([self.altitude], dtype='float64')

        metadata = get_metadata('metadata')
        metadata['original_container'] = 'CINRAD/CC'
        metadata['site_name'] = self.sitename
        metadata['radar_name'] = "CINRAD/CC"

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
        pulse_width['data'] = self.CC.header['CutConfig']['usBindWidth'][0] * 2. / _LIGHT_SPEED  # m->sec
        # assume that the parameters in the first ray represent the beam widths,
        # bandwidth and frequency in the entire volume
        wavelength_hz = self.frequency * 10 ** 9
        # radar_beam_width_h
        radar_beam_width_h = get_metadata('radar_beam_width_h')
        radar_beam_width_h['data'] = np.array([0.703125, ], dtype='float32')
        # radar_beam_width_v
        radar_beam_width_v = get_metadata('radar_beam_width_v')
        radar_beam_width_v['data'] = np.array([0.703125, ], dtype='float32')
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
