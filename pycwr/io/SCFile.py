# -*- coding: utf-8 -*-
import numpy as np
from .BaseDataProtocol.SCProtocol import dtype_sc
from .util import _prepare_for_read, _unpack_from_buf, make_time_unit_str, get_radar_sitename
import pandas as pd
import datetime
from ..core.NRadar import PRD
from ..configure.pyart_config import get_metadata, get_fillvalue
from ..configure.default_config import CINRAD_field_mapping, _LIGHT_SPEED
from ..core.PyartRadar import Radar
from netCDF4 import date2num

class SCBaseData(object):
    """
    解码SC/CD 1.0的数据格式
    """
    def __init__(self, filename, station_lon=None, station_lat=None, station_alt=None):
        """
        :param filename:  radar basedata filename
        :param station_lon:  radar station longitude //units: degree east
        :param station_lat:  radar station latitude //units:degree north
        :param station_alt:  radar station altitude //units: meters
        """
        super(SCBaseData, self).__init__()
        self.filename = filename
        self.station_lon = station_lon
        self.station_lat = station_lat
        self.station_alt = station_alt
        self.fid = _prepare_for_read(self.filename) ##判断是否是压缩文件
        buf_header = self.fid.read(dtype_sc.BaseDataHeaderSize) ##取出header的buf
        self.header = self._parse_BaseDataHeader(buf_header)
        self.MaxV = self.header['LayerParam']['MaxV'][0]/100. ##??可能会存在问题，如果不同仰角采用不用的PRF
        self._check_sc_basedata()
        self.fid.seek(dtype_sc.BaseDataHeaderSize, 0) ##移动到径向数据的位置
        self.radial = self._parse_radial()
        self.fid.close()

    def _check_sc_basedata(self):
        """检查雷达数据是否完整"""
        buf_radial_data = self.fid.read()
        assert len(buf_radial_data) == self.nrays * dtype_sc.PerRadialSize, "SC basedata size has problems!"
        return

    def _parse_BaseDataHeader(self, buf_header):
        """
        :param buf_header: 只包含头文件的buf
        :return:
        """
        BaseDataHeader_dict = {}
        ##解码雷达站点信息
        BaseDataHeader_dict['RadarSite'], _ = _unpack_from_buf(buf_header,\
        dtype_sc.RadarSitePos,dtype_sc.BaseDataHeader['RadarSite'])
        ##解码雷达性能参数
        BaseDataHeader_dict['RadarPerformanceParam'], _ = _unpack_from_buf(buf_header,\
        dtype_sc.RadarPerformanceParamPos,dtype_sc.BaseDataHeader['RadarPerformanceParam'])
        ##解码观测参数
        BaseDataHeader_dict['RadarObserationParam_1'], _ = _unpack_from_buf(buf_header, \
        dtype_sc.RadarObserationParamPos_1, dtype_sc.BaseDataHeader['RadarObserationParam_1'])
        ##目前仅支持VOL扫描格式
        assert BaseDataHeader_dict['RadarObserationParam_1']['stype'] > 100, "only vol support!"
        self.nsweeps = BaseDataHeader_dict['RadarObserationParam_1']['stype'] - 100
        ##解码不同仰角的观测参数
        BaseDataHeader_dict['LayerParam'] = np.frombuffer(buf_header, \
        dtype_sc.BaseDataHeader['LayerParamX30'],count=self.nsweeps, offset=dtype_sc.LayerParamPos)
        ####################sc basedata has wrong bandwidth############
        BaseDataHeader_dict["binWidth"] = np.full_like(BaseDataHeader_dict['LayerParam']["binWidth"], 5000, dtype=np.int32)
        ####################sc basedata has wrong bandwidth############
        ##解码其余一些观测参数
        BaseDataHeader_dict['RadarObserationParam_2'], _ = _unpack_from_buf(buf_header,\
            dtype_sc.RadarObserationParamPos_2, dtype_sc.BaseDataHeader['RadarObserationParam_2'])
        self.nrays = np.sum((BaseDataHeader_dict['LayerParam']['recordnumber']).astype(np.int64))
        self.sweep_end_ray_index_add1 = np.cumsum((BaseDataHeader_dict['LayerParam']['recordnumber']).astype(np.int64)) ##python格式的结束
        self.sweep_start_ray_index = self.sweep_end_ray_index_add1 - \
                                     (BaseDataHeader_dict['LayerParam']['recordnumber']).astype(np.int64)
        return BaseDataHeader_dict

    def _parse_radial(self):
        radial = []
        for isweep in range(self.nsweeps):
            MaxV = self.header['LayerParam']['MaxV'][isweep] / 100.
            for _ in range(self.header['LayerParam']['recordnumber'][isweep]):
                buf_radial = self.fid.read(dtype_sc.PerRadialSize)
                radial.append(self._parse_radial_single(buf_radial, 0, MaxV, -1))
        return radial

    def _parse_radial_single(self, buf_radial, start_pos, MaxV, num_bins=-1):
        """
        :param buf_radial:
        :param start_pos: 开始的pos
        :param num_bins: 库数
        :return:
        """
        radial_dict = {}
        radial_dict_tmp, size_tmp = _unpack_from_buf(buf_radial, start_pos, dtype_sc.RadialHeader())
        radial_dict.update(radial_dict_tmp)
        radial_dict['fields'] = {}
        ##basedata save some error data
        #RadialData = np.frombuffer(buf_radial, dtype_sc.RadialData(),\
        #                                 count=num_bins, offset=start_pos+size_tmp)
        RadialData = np.frombuffer(buf_radial, dtype_sc.RadialData(),\
                                         count=500, offset=start_pos+size_tmp)
        radial_dict['fields']['dBZ'] = np.where(RadialData['dBZ'] != 0,\
                            (RadialData['dBZ'].astype(int) - 64)/2., np.nan).astype(np.float32)
        radial_dict['fields']['dBT'] = np.where(RadialData['dBT'] != 0, \
                        (RadialData['dBT'].astype(int) - 64) / 2., np.nan).astype(np.float32)
        radial_dict['fields']['V'] = np.where(RadialData['V'] != 0, \
                       MaxV * (RadialData['V'].astype(int) - 128) / 128., np.nan).astype(np.float32)
        radial_dict['fields']['W'] = np.where(RadialData['W'] != 0, \
                        MaxV * RadialData['W'].astype(int) /256., np.nan).astype(np.float32)
        return radial_dict

    def get_nyquist_velocity(self):
        """get nyquist vel per ray
        获取每根径向的不模糊速度
        :return:(nRays)
        """
        nyquist_velocity = np.concatenate([np.array([self.header['LayerParam']['MaxV'][isweep] / \
                          100.] * self.header['LayerParam']['recordnumber'][isweep]) for \
                            isweep in range(self.nsweeps)])
        return nyquist_velocity.astype(np.float32)
    def get_unambiguous_range(self):
        """
        获取每根径向的不模糊距离
        :return:(nRays)
        """
        return np.concatenate([np.array([self.header['LayerParam']['MaxL'][isweep] *10 \
                            ] * self.header['LayerParam']['recordnumber'][isweep]) for \
                            isweep in range(self.nsweeps)])

    def get_scan_time(self):
        """
        获取每根径向的扫描时间
        :return:(nRays)
        """
        Start_params = self.header['RadarObserationParam_1']
        End_params = self.header['RadarObserationParam_2']
        start_time = datetime.datetime(year=Start_params['syear'], month=Start_params['smonth'],
                                       day=Start_params['sday'], hour=Start_params['shour'],
                                       minute=Start_params['sminute'], second=Start_params['ssecond'])
        end_time = datetime.datetime(year=End_params['Eyear'], month=End_params['Emonth'],
                                       day=End_params['Eday'], hour=End_params['Ehour'],
                                       minute=End_params['Eminute'], second=End_params['Esecond'])
        return pd.date_range(start_time, end_time, periods=self.nrays).to_pydatetime() - datetime.timedelta(hours=8)

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
        return self.sweep_end_ray_index_add1 - self.sweep_start_ray_index

    def get_azimuth(self):
        """
        获取每根径向的方位角
        :return:(nRays)
        """
        return np.concatenate([np.arange(0,360,1.0), ] * self.nsweeps)

    def get_elevation(self):
        """
        获取每根径向的仰角
        :return: (nRays)
        """
        return np.array([(iray['sStrEl'] + iray['sEndEl'])*180./65536 for iray in self.radial])

    def get_latitude_longitude_altitude_frequency(self):
        """
        获取经纬度高度，雷达频率
        :return:lat, lon, alt, frequency
        """
        lat, lon, alt, frequency =  self.header['RadarSite']['latitudevalue']/100., \
                                    self.header['RadarSite']['longitudevalue']/100., \
                                    self.header['RadarSite']['height'] / 1000., 2.765
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
        if self.header['RadarObserationParam_1']['stype'] == 1:
            return "rhi"
        else:
            return "ppi"
    def get_sitename(self):
        return get_radar_sitename(self.filename)

class SC2NRadar(object):
    """到NusitRadar object 的桥梁"""

    def __init__(self, SC):
        self.SC = SC
        self.radial = self.SC.radial
        self.azimuth = self.get_azimuth()
        self.elevation = self.get_elevation()
        self.sweep_start_ray_index = self.get_sweep_start_ray_index()
        self.sweep_end_ray_index = self.get_sweep_end_ray_index()
        self.nrays = self.SC.nrays
        self.nsweeps = self.SC.nsweeps
        self.scan_type = self.SC.get_scan_type()
        self.latitude, self.longitude, self.altitude, self.frequency = \
            self.SC.get_latitude_longitude_altitude_frequency()
        self.bins_per_sweep = self.get_nbins_per_sweep()
        self.max_bins = self.bins_per_sweep.max()
        self.range = self.get_range_per_radial(self.max_bins)
        self.fields = self._get_fields()
        self.sitename = self.SC.get_sitename()

    def get_azimuth(self):
        """
        获取每根径向的方位角
        :return:(nRays)
        """
        return self.SC.get_azimuth()

    def get_elevation(self):
        """
        获取每根径向的仰角
        :return: (nRays)
        """
        return self.SC.get_elevation()

    def get_rays_per_sweep(self):
        """
        获取每个sweep的径向数
        :return:(nsweep)
        """
        return (self.SC.header['LayerParam']['recordnumber']).astype(int)

    def get_scan_time(self):
        """
        获取每根径向的扫描时间
        :return:(nRays)
        """
        return self.SC.get_scan_time()

    def get_nyquist_velocity(self):
        """get nyquist vel per ray
        获取每根径向的不模糊速度
        :return:(nRays)
        """
        return self.SC.get_nyquist_velocity()

    def get_unambiguous_range(self):
        """
        获取每根径向的不模糊距离
        :return:(nRays)
        """
        return self.SC.get_unambiguous_range()

    def get_sweep_end_ray_index(self):
        """
        获取每个sweep的结束的index，包含在内
        :return:(nsweep)
        """
        return self.SC.sweep_end_ray_index_add1 - 1

    def get_sweep_start_ray_index(self):
        """
        获取每个sweep的开始的index
        :return:(nsweep)
        """
        return self.SC.sweep_start_ray_index

    def get_nbins_per_sweep(self):
        """
        确定每个sweep V探测的库数
        :return:
        """
        return np.array([self.radial[idx]['fields']['dBZ'].size for idx in self.sweep_start_ray_index])

    def get_range_per_radial(self, length):
        """
        确定径向每个库的距离
        :param length:
        :return:
        """
        Resolution = self.SC.header['binWidth'][0]/10.
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
        return self.SC.header['LayerParam']['MaxV'] / 100.

    def get_NRadar_unambiguous_range(self):
        """array shape (nsweeps)"""
        return self.SC.header['LayerParam']['MaxL'] * 10.

    def get_fixed_angle(self):
        return self.SC.header['LayerParam']['Swangles'] / 100.

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
        _range['meters_to_center_of_first_gate'] = self.range[0]
        _range['meters_between_gates'] = self.range[0]

        latitude = get_metadata('latitude')
        longitude = get_metadata('longitude')
        altitude = get_metadata('altitude')
        latitude['data'] = np.array([self.latitude], dtype='float64')
        longitude['data'] = np.array([self.longitude], dtype='float64')
        altitude['data'] = np.array([self.altitude], dtype='float64')

        metadata = get_metadata('metadata')
        metadata['original_container'] = 'CINRAD/SC'
        metadata['site_name'] = self.sitename
        metadata['radar_name'] = "CINRAD/SC"

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
        pulse_width['data'] = np.array([self.range[0] / _LIGHT_SPEED,], dtype='float32')  # nanosec->sec
        # assume that the parameters in the first ray represent the beam widths,
        # bandwidth and frequency in the entire volume
        wavelength_hz = self.frequency * 10 ** 9
        # radar_beam_width_h
        radar_beam_width_h = get_metadata('radar_beam_width_h')
        radar_beam_width_h['data'] = np.array([1., ], dtype='float32')
        # radar_beam_width_v
        radar_beam_width_v = get_metadata('radar_beam_width_v')
        radar_beam_width_v['data'] = np.array([1., ], dtype='float32')
        # frequency
        frequency = get_metadata('frequency')
        frequency['data'] = np.array([wavelength_hz, ], dtype='float32')
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
