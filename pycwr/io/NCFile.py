# -*- coding: utf-8 -*-
import xarray as xr
import numpy as np
from datetime import datetime, timedelta
from .util import make_time_unit_str
from ..core.NRadar import PRD
from ..configure.pyart_config import get_metadata, get_fillvalue
from ..configure.default_config import CINRAD_field_mapping, _LIGHT_SPEED
from ..core.PyartRadar import Radar
from netCDF4 import date2num
import pandas as pd

class NC2NRadar(object):
    """到NusitRadar object 的桥梁"""

    def __init__(self, filename, band = "S", ScanDuration=6):

        assert band in ["S", "C", "X"], "band not in S, C, X!"
        self.NC = xr.open_dataset(filename)
        self.StartScanTime = datetime.fromtimestamp(self.NC.RadTime) - timedelta(hours=8)
        self.azimuth = np.expand_dims(self.NC.Dim2.values, 0).repeat(self.NC.Dim1.shape[0], axis=0).ravel()
        self.elevation = np.expand_dims(self.NC.Dim1.values, 1).repeat(self.NC.Dim2.shape[0], axis=1).ravel()
        self.sweep_end_ray_index = np.cumsum(np.array([self.NC.Dim2.shape[0], ] * self.NC.Dim1.shape[0])).astype(int)
        self.sweep_start_ray_index = self.sweep_end_ray_index - self.NC.Dim2.shape[0]
        self.nrays = self.NC.Dim1.shape[0] * self.NC.Dim2.shape[0]
        self.nsweeps = self.NC.Dim1.shape[0]
        self.scan_type = "ppi"
        self.latitude = self.NC.RadLat.values
        self.longitude = self.NC.RadLon.values
        self.altitude = self.NC.RadHgt.values
        if band == "S":
            self.frequency = 3.0
        elif band == "C":
            self.frequency = 5.7
        elif band == "X":
            self.frequency = 9.5
        self.bins_per_sweep = np.array([self.NC.Dim2.shape[0], ] * self.NC.Dim1.shape[0])
        self.max_bins = self.bins_per_sweep.max()
        self.range = self.NC.Dim3.values
        self.fields = self._get_fields()
        self.sitename = "Unknown"
        self.scan_time = (pd.timedelta_range(start="0 min", end="%d min"%ScanDuration, periods=self.nrays) + self.StartScanTime).to_pydatetime()
        # print(self.scan_time)
        self.fixed_angle = self.NC.Dim1.values
        self.nyquist_velocity = self.NC.VMax.values #np.expand_dims(self.NC.VMax.values, 1).repeat(self.NC.Dim2.shape[0], axis=1).ravel()
        self.unambiguous_range = np.full_like(self.nyquist_velocity, self.range.max())

    def _get_fields(self):
        """将所有的field的数据提取出来"""
        tmp_nc = self.NC.drop_vars("tim")
        field_keys = list([ikey for ikey in tmp_nc.keys() if tmp_nc[ikey].ndim == 3])
        vars2keys = {"ref": "dBZ", "vel": "V", "wid": "W", "zdr": "ZDR", "rhv": "CC",
                     "pdp": "PhiDP", "kdp": "KDP"}
        fields = {}
        for ikey in field_keys:
            # print(tmp_nc[ikey])
            fields[vars2keys[ikey]] = tmp_nc[ikey].values.reshape(-1, tmp_nc.Dim3.shape[0])
        return fields
    def get_nyquist_velocity(self):
        return np.expand_dims(self.NC.VMax.values, 1).repeat(self.NC.Dim2.shape[0], axis=1).ravel()

    def ToPRD(self):
        """将WSR98D数据转为PRD 的数据格式"""

        return PRD(fields=self.fields, scan_type=self.scan_type, time=self.scan_time, \
                          range=self.range, azimuth=self.azimuth, elevation=self.elevation, latitude=self.latitude, \
                          longitude=self.longitude, altitude=self.altitude,
                          sweep_start_ray_index=self.sweep_start_ray_index, \
                          sweep_end_ray_index=self.sweep_end_ray_index, fixed_angle=self.fixed_angle, \
                          bins_per_sweep=self.bins_per_sweep, nyquist_velocity=self.nyquist_velocity, \
                          frequency=self.frequency, unambiguous_range=self.unambiguous_range, \
                          nrays=self.nrays, nsweeps=self.nsweeps, sitename = self.sitename, pyart_radar=self.ToPyartRadar())

    def ToPyartRadar(self):

        dts = self.scan_time
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
        _range['meters_between_gates'] = self.range[1] - self.range[0]

        latitude = get_metadata('latitude')
        longitude = get_metadata('longitude')
        altitude = get_metadata('altitude')
        latitude['data'] = np.array([self.latitude], dtype='float64')
        longitude['data'] = np.array([self.longitude], dtype='float64')
        altitude['data'] = np.array([self.altitude], dtype='float64')

        metadata = get_metadata('metadata')
        metadata['original_container'] = 'CINRAD'
        metadata['site_name'] = self.sitename
        metadata['radar_name'] = "CINRAD/Unknown"

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
        fixed_angle['data'] = self.fixed_angle

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
        pulse_width['data'] = (self.range[1] - self.range[0]) * 2. / _LIGHT_SPEED  # m->sec
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

if __name__ == "__main__":
    ds = xr.open_dataset("/Users/zhengyu/Downloads/Z_RADR_I_70706_20220422152835_O_DOR_YLD_CAP_FMT.nc")
    ds1 = xr.open_dataset("/Users/zhengyu/Downloads/Z_RADR_I_FSM00_20210902123230_O_DOR_ETW_CAP_FMT.nc")