# -*- coding: utf-8 -*-
import numpy as np
from .BaseDataProtocol.SCProtocol import dtype_sc
from .util import _prepare_for_read, _read_all, _read_exact, _unpack_from_buf, make_time_unit_str, get_radar_sitename, date2num
import pandas as pd
import datetime
from ..core.NRadar import PRD
from ..configure.pyart_config import get_metadata, get_fillvalue
from ..configure.default_config import CINRAD_field_mapping, _LIGHT_SPEED
from ..core.PyartRadar import Radar

class SCBaseData(object):
    """Decode SC/CD 1.0 base data."""
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
        self.fid = _prepare_for_read(self.filename)
        buf_header = _read_exact(self.fid, dtype_sc.BaseDataHeaderSize, "SC header")
        self.header = self._parse_BaseDataHeader(buf_header)
        self.MaxV = self.header['LayerParam']['MaxV'][0] / 100.
        self._radial_buf = self._check_sc_basedata()
        self.radial = self._parse_radial()
        self.fid.close()

    def _check_sc_basedata(self):
        """Check that the radial payload length matches the header metadata."""
        buf_radial_data = _read_all(self.fid, "SC radial payload")
        expected = self.nrays * dtype_sc.PerRadialSize
        if len(buf_radial_data) != expected:
            raise ValueError("SC radial payload size does not match the header metadata.")
        return buf_radial_data

    def _parse_BaseDataHeader(self, buf_header):
        """
        :param buf_header: header-only byte buffer
        :return:
        """
        BaseDataHeader_dict = {}
        # Decode radar site information.
        BaseDataHeader_dict['RadarSite'], _ = _unpack_from_buf(buf_header,\
        dtype_sc.RadarSitePos,dtype_sc.BaseDataHeader['RadarSite'])
        # Decode radar performance parameters.
        BaseDataHeader_dict['RadarPerformanceParam'], _ = _unpack_from_buf(buf_header,\
        dtype_sc.RadarPerformanceParamPos,dtype_sc.BaseDataHeader['RadarPerformanceParam'])
        # Decode the first observation-parameter block.
        BaseDataHeader_dict['RadarObserationParam_1'], _ = _unpack_from_buf(buf_header, \
        dtype_sc.RadarObserationParamPos_1, dtype_sc.BaseDataHeader['RadarObserationParam_1'])
        # Only volume scans are supported here.
        if BaseDataHeader_dict['RadarObserationParam_1']['stype'] <= 100:
            raise ValueError("Only volume-scan SC files are supported.")
        self.nsweeps = BaseDataHeader_dict['RadarObserationParam_1']['stype'] - 100
        # Decode per-sweep layer parameters.
        BaseDataHeader_dict['LayerParam'] = np.frombuffer(buf_header, \
        dtype_sc.BaseDataHeader['LayerParamX30'],count=self.nsweeps, offset=dtype_sc.LayerParamPos)
        # SC base data carries incorrect bin width metadata, keep the historical override.
        BaseDataHeader_dict["binWidth"] = np.full_like(BaseDataHeader_dict['LayerParam']["binWidth"], 5000, dtype=np.int32)
        # Decode the remaining observation parameters.
        BaseDataHeader_dict['RadarObserationParam_2'], _ = _unpack_from_buf(buf_header,\
            dtype_sc.RadarObserationParamPos_2, dtype_sc.BaseDataHeader['RadarObserationParam_2'])
        self.nrays = np.sum((BaseDataHeader_dict['LayerParam']['recordnumber']).astype(np.int64))
        self.sweep_end_ray_index_add1 = np.cumsum((BaseDataHeader_dict['LayerParam']['recordnumber']).astype(np.int64))
        self.sweep_start_ray_index = self.sweep_end_ray_index_add1 - \
                                     (BaseDataHeader_dict['LayerParam']['recordnumber']).astype(np.int64)
        return BaseDataHeader_dict

    def _parse_radial(self):
        radial = []
        payload = self._radial_buf
        pos = 0
        for isweep in range(self.nsweeps):
            MaxV = self.header['LayerParam']['MaxV'][isweep] / 100.
            for _ in range(self.header['LayerParam']['recordnumber'][isweep]):
                buf_radial = payload[pos:pos + dtype_sc.PerRadialSize]
                pos += dtype_sc.PerRadialSize
                radial.append(self._parse_radial_single(buf_radial, 0, MaxV, -1))
        return radial

    def _parse_radial_single(self, buf_radial, start_pos, MaxV, num_bins=-1):
        """
        :param buf_radial:
        :param start_pos: starting byte offset.
        :param num_bins: number of bins.
        :return:
        """
        radial_dict = {}
        radial_dict_tmp, size_tmp = _unpack_from_buf(buf_radial, start_pos, dtype_sc.RadialHeader())
        radial_dict.update(radial_dict_tmp)
        radial_dict['fields'] = {}
        # Historical files may contain invalid trailing radial values.
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
        """Return the per-ray Nyquist velocity."""
        nyquist_velocity = np.concatenate([np.array([self.header['LayerParam']['MaxV'][isweep] / \
                          100.] * self.header['LayerParam']['recordnumber'][isweep]) for \
                            isweep in range(self.nsweeps)])
        return nyquist_velocity.astype(np.float32)
    def get_unambiguous_range(self):
        """Return the per-ray unambiguous range."""
        return np.concatenate([np.array([self.header['LayerParam']['MaxL'][isweep] *10 \
                            ] * self.header['LayerParam']['recordnumber'][isweep]) for \
                            isweep in range(self.nsweeps)])

    def get_scan_time(self):
        """Return the acquisition time for each ray."""
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
        """Return the inclusive end index of each sweep."""
        return self.sweep_end_ray_index_add1 - 1
    def get_sweep_start_ray_index(self):
        """Return the start index of each sweep."""
        return self.sweep_start_ray_index

    def get_rays_per_sweep(self):
        """Return the number of rays in each sweep."""
        return self.sweep_end_ray_index_add1 - self.sweep_start_ray_index

    def get_azimuth(self):
        """Return the azimuth angle for each ray."""
        return np.concatenate([np.arange(0,360,1.0), ] * self.nsweeps)

    def get_elevation(self):
        """Return the elevation angle for each ray."""
        return np.array([(iray['sStrEl'] + iray['sEndEl'])*180./65536 for iray in self.radial])

    def get_latitude_longitude_altitude_frequency(self):
        """Return latitude, longitude, altitude, and radar frequency."""
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
        """Return the scan type string."""
        if self.header['RadarObserationParam_1']['stype'] == 1:
            return "rhi"
        else:
            return "ppi"
    def get_sitename(self):
        return get_radar_sitename(self.filename)

class SC2NRadar(object):
    """Bridge from raw SC data to an NRadar object."""

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
        """Return the azimuth angle for each ray."""
        return self.SC.get_azimuth()

    def get_elevation(self):
        """Return the elevation angle for each ray."""
        return self.SC.get_elevation()

    def get_rays_per_sweep(self):
        """Return the number of rays in each sweep."""
        return (self.SC.header['LayerParam']['recordnumber']).astype(int)

    def get_scan_time(self):
        """Return the acquisition time for each ray."""
        return self.SC.get_scan_time()

    def get_nyquist_velocity(self):
        """Return the per-ray Nyquist velocity."""
        return self.SC.get_nyquist_velocity()

    def get_unambiguous_range(self):
        """Return the per-ray unambiguous range."""
        return self.SC.get_unambiguous_range()

    def get_sweep_end_ray_index(self):
        """Return the inclusive end index of each sweep."""
        return self.SC.sweep_end_ray_index_add1 - 1

    def get_sweep_start_ray_index(self):
        """Return the start index of each sweep."""
        return self.SC.sweep_start_ray_index

    def get_nbins_per_sweep(self):
        """Return the gate count for each sweep."""
        return np.array([self.radial[idx]['fields']['dBZ'].size for idx in self.sweep_start_ray_index])

    def get_range_per_radial(self, length):
        """Return gate-center ranges for a radial of ``length`` bins."""
        Resolution = self.SC.header['binWidth'][0]/10.
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
        return self.SC.header['LayerParam']['MaxV'] / 100.

    def get_NRadar_nyquist_speed(self):
        """Backward-compatible alias for ``get_nradar_nyquist_speed``."""
        return self.get_nradar_nyquist_speed()

    def get_nradar_unambiguous_range(self):
        """array shape (nsweeps)"""
        return self.SC.header['LayerParam']['MaxL'] * 10.

    def get_NRadar_unambiguous_range(self):
        """Backward-compatible alias for ``get_nradar_unambiguous_range``."""
        return self.get_nradar_unambiguous_range()

    def get_fixed_angle(self):
        return self.SC.header['LayerParam']['Swangles'] / 100.

    def to_prd(self, effective_earth_radius=None):
        """Build an ``NRadar.PRD`` object from the decoded SC volume."""

        return PRD(fields=self.fields, scan_type=self.scan_type, time=self.get_scan_time(), \
                          range=self.range, azimuth=self.azimuth, elevation=self.elevation, latitude=self.latitude, \
                          longitude=self.longitude, altitude=self.altitude,
                          sweep_start_ray_index=self.sweep_start_ray_index, \
                          sweep_end_ray_index=self.sweep_end_ray_index, fixed_angle=self.get_fixed_angle(), \
                          bins_per_sweep=self.bins_per_sweep, nyquist_velocity=self.get_nradar_nyquist_speed(), \
                          frequency=self.frequency, unambiguous_range=self.get_nradar_unambiguous_range(), \
                          nrays=self.nrays, nsweeps=self.nsweeps, sitename = self.sitename,
                          pyart_radar=None, effective_earth_radius=effective_earth_radius,
                          metadata={"original_container": "CINRAD/SC", "radar_name": "CINRAD/SC"})

    def ToPRD(self, effective_earth_radius=None):
        """Backward-compatible alias for ``to_prd``."""
        return self.to_prd(effective_earth_radius=effective_earth_radius)

    def to_pyart_radar(self, effective_earth_radius=None, **kwargs):
        """Export the decoded SC volume through the PRD Py-ART adapter."""
        return self.to_prd(effective_earth_radius=effective_earth_radius).to_pyart_radar(**kwargs)

    def ToPyartRadar(self, effective_earth_radius=None, **kwargs):
        """Backward-compatible alias for ``to_pyart_radar``."""
        return self.to_pyart_radar(effective_earth_radius=effective_earth_radius, **kwargs)

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
