# -*- coding: utf-8 -*-
import numpy as np
from .BaseDataProtocol.CCProtocol import dtype_cc
from .util import _prepare_for_read, _read_all, _read_exact, _unpack_from_buf, make_time_unit_str, get_radar_sitename, date2num
import datetime
import pandas as pd
from ..core.NRadar import PRD
from ..configure.pyart_config import get_metadata, get_fillvalue
from ..configure.default_config import CINRAD_field_mapping, _LIGHT_SPEED
from ..core.PyartRadar import Radar

class CCBaseData(object):
    """Decode CC/CCJ base data."""

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
        self.fid = _prepare_for_read(self.filename)
        buf_header = _read_exact(self.fid, dtype_cc.BaseDataHeaderSize, "CC header")
        self.header = self._parse_BaseDataHeader(buf_header)
        self._radial_buf = self._check_cc_basedata()
        self.radial = self._parse_radial()
        self.fid.close()

    def _check_cc_basedata(self):
        """Check that the radial payload length matches the header metadata."""
        buf_radial_data = _read_all(self.fid, "CC radial payload")
        expected = self.nrays * dtype_cc.PerRadialSize
        if len(buf_radial_data) != expected:
            raise ValueError("CC radial payload size does not match the header metadata.")
        return buf_radial_data

    def _parse_BaseDataHeader(self, buf_header):
        """
        :param buf_header: header-only byte buffer
        :return:
        """
        BaseDataHeader_dict = {}
        # Decode the first observation-parameter block.
        BaseDataHeader_dict['ObsParam1'], _ = _unpack_from_buf(buf_header, \
                                                               dtype_cc.HeaderSize1_pos,
                                                               dtype_cc.BaseDataHeader['RadarHeader1'])
        # Decode per-sweep cut configuration.
        if BaseDataHeader_dict['ObsParam1']['ucScanMode'] <= 100:
            raise ValueError("Only volume-scan CC files are supported.")
        self.nsweeps = BaseDataHeader_dict['ObsParam1']['ucScanMode'] - 100
        BaseDataHeader_dict['CutConfig'] = np.frombuffer(buf_header, \
                                                         dtype_cc.BaseDataHeader['CutConfigX30'], count=self.nsweeps,
                                                         offset=dtype_cc.CutSize_pos)
        # Decode the second observation-parameter block.
        BaseDataHeader_dict['ObsParam2'], _ = _unpack_from_buf(buf_header, \
                                                               dtype_cc.HeaderSize2_pos,
                                                               dtype_cc.BaseDataHeader['RadarHeader2'])

        self.nrays = np.sum(BaseDataHeader_dict['CutConfig']['usRecordNumber'])
        self.sweep_end_ray_index_add1 = (np.cumsum(BaseDataHeader_dict['CutConfig']['usRecordNumber'])).astype(int)
        self.sweep_start_ray_index = (self.sweep_end_ray_index_add1 - BaseDataHeader_dict['CutConfig']['usRecordNumber']).astype(int)
        return BaseDataHeader_dict

    def _parse_radial_single(self, buf_radial, radialnumber):
        """Parse one CC radial."""
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
        payload = self._radial_buf
        pos = 0
        for isweep in range(self.nsweeps):
            radialnumber = self.header['CutConfig']['usBinNumber'][isweep]
            for _ in range(self.header['CutConfig']['usRecordNumber'][isweep]):
                buf_radial = payload[pos:pos + dtype_cc.PerRadialSize]
                pos += dtype_cc.PerRadialSize
                radial.append(self._parse_radial_single(buf_radial, radialnumber))
        return radial

    def get_nyquist_velocity(self):
        """Return the per-ray Nyquist velocity."""
        nyquist_velocity = np.concatenate([np.array([self.header['CutConfig']['usMaxV'][isweep] / 100.] \
                                                    * self.header['CutConfig']['usRecordNumber'][isweep]) for \
                                           isweep in range(self.nsweeps)])
        return nyquist_velocity.astype(np.float32)

    def get_unambiguous_range(self):
        """Return the per-ray unambiguous range."""
        return np.concatenate([np.array([self.header['CutConfig']['usMaxL'][isweep] * 10. \
                                         ] * self.header['CutConfig']['usRecordNumber'][isweep]) for \
                               isweep in range(self.nsweeps)])

    def get_scan_time(self):
        """Return the acquisition time for each ray."""
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
        """Return the inclusive end index of each sweep."""
        return self.sweep_end_ray_index_add1 - 1

    def get_sweep_start_ray_index(self):
        """Return the start index of each sweep."""
        return self.sweep_start_ray_index

    def get_rays_per_sweep(self):
        """Return the number of rays in each sweep."""
        return (self.header['CutConfig']['usRecordNumber']).astype(np.int32)

    def get_azimuth(self):
        """Return the azimuth angle for each ray."""
        return np.concatenate([np.linspace(0, 360, self.header['CutConfig']['usRecordNumber'][isweep]) \
                               for isweep in range(self.nsweeps)], axis=0)

    def get_elevation(self):
        """Return the elevation angle for each ray."""
        return np.concatenate([np.array([self.header['CutConfig']['usAngle'][isweep] / 100. \
                                         ] * self.header['CutConfig']['usRecordNumber'][isweep]) for \
                               isweep in range(self.nsweeps)])

    def get_latitude_longitude_altitude_frequency(self):
        """Return latitude, longitude, altitude, and radar frequency."""
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
        """Return the scan type string."""
        # Only PPI is supported in the current CC reader.
        return "ppi"

    def get_sitename(self):
        return get_radar_sitename(self.filename)

class CC2NRadar(object):
    """Bridge from raw CC data to an NRadar object."""

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
        """Return the azimuth angle for each ray."""
        return self.CC.get_azimuth()

    def get_elevation(self):
        """Return the elevation angle for each ray."""
        return self.CC.get_elevation()

    def get_rays_per_sweep(self):
        """Return the number of rays in each sweep."""
        return (self.CC.header['CutConfig']['usRecordNumber']).astype(int)

    def get_scan_time(self):
        """Return the acquisition time for each ray."""
        return self.CC.get_scan_time()

    def get_nyquist_velocity(self):
        """Return the per-ray Nyquist velocity."""
        return self.CC.get_nyquist_velocity()

    def get_unambiguous_range(self):
        """Return the per-ray unambiguous range."""
        return self.CC.get_unambiguous_range()

    def get_sweep_end_ray_index(self):
        """Return the inclusive end index of each sweep."""
        return self.CC.sweep_end_ray_index_add1 - 1

    def get_sweep_start_ray_index(self):
        """Return the start index of each sweep."""
        return self.CC.sweep_start_ray_index

    def get_nbins_per_sweep(self):
        """Return the gate count for each sweep."""
        return (self.CC.header['CutConfig']['usBinNumber']).astype(int)

    def get_range_per_radial(self, length):
        """Return gate-center ranges for a radial of ``length`` bins."""
        Resolution = self.CC.header['CutConfig']['usBindWidth'][0] * 2
        return np.linspace(Resolution, Resolution * length, length)

    def _get_fields(self):
        """Assemble all fields into dense 2-D arrays."""
        fields = {}
        field_keys = self.radial[0]['fields'].keys()
        for ikey in field_keys:
            fields[ikey] = np.array([(iray['fields'][ikey]).ravel() for iray in self.radial])
        return fields

    def get_nradar_nyquist_speed(self):
        """array shape (nsweeps)"""
        return self.CC.header['CutConfig']['usMaxV'] / 100.

    def get_NRadar_nyquist_speed(self):
        """Backward-compatible alias for ``get_nradar_nyquist_speed``."""
        return self.get_nradar_nyquist_speed()

    def get_nradar_unambiguous_range(self):
        """array shape (nsweeps)"""
        return self.CC.header['CutConfig']['usMaxL'] * 10.

    def get_NRadar_unambiguous_range(self):
        """Backward-compatible alias for ``get_nradar_unambiguous_range``."""
        return self.get_nradar_unambiguous_range()

    def get_fixed_angle(self):
        return self.CC.header['CutConfig']['usAngle'] / 100.

    def to_prd(self, effective_earth_radius=None):
        """Build an ``NRadar.PRD`` object from the decoded CC volume."""

        return PRD(fields=self.fields, scan_type=self.scan_type, time=self.get_scan_time(), \
                          range=self.range, azimuth=self.azimuth, elevation=self.elevation, latitude=self.latitude, \
                          longitude=self.longitude, altitude=self.altitude,
                          sweep_start_ray_index=self.sweep_start_ray_index, \
                          sweep_end_ray_index=self.sweep_end_ray_index, fixed_angle=self.get_fixed_angle(), \
                          bins_per_sweep=self.bins_per_sweep, nyquist_velocity=self.get_nradar_nyquist_speed(), \
                          frequency=self.frequency, unambiguous_range=self.get_nradar_unambiguous_range(), \
                          nrays=self.nrays, nsweeps=self.nsweeps, sitename = self.sitename,
                          pyart_radar=None, effective_earth_radius=effective_earth_radius,
                          metadata={"original_container": "CINRAD/CC", "radar_name": "CINRAD/CC"})

    def ToPRD(self, effective_earth_radius=None):
        """Backward-compatible alias for ``to_prd``."""
        return self.to_prd(effective_earth_radius=effective_earth_radius)

    def to_pyart_radar(self, effective_earth_radius=None, **kwargs):
        """Export the decoded CC volume through the PRD Py-ART adapter."""
        return self.to_prd(effective_earth_radius=effective_earth_radius).to_pyart_radar(**kwargs)

    def ToPyartRadar(self, effective_earth_radius=None, **kwargs):
        """Backward-compatible alias for ``to_pyart_radar``."""
        return self.to_pyart_radar(effective_earth_radius=effective_earth_radius, **kwargs)

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
