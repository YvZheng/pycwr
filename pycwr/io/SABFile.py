# -*- coding: utf-8 -*-
import numpy as np
from .BaseDataProtocol.SABProtocol import dtype_sab
from .util import _prepare_for_read, _read_all, _unpack_from_buf, _validate_offset_length, julian2date, \
    get_radar_info, make_time_unit_str, get_radar_sitename, date2num
from ..core.NRadar import PRD
from ..configure.pyart_config import get_metadata, get_fillvalue
from ..configure.default_config import CINRAD_field_mapping, _LIGHT_SPEED
from ..core.PyartRadar import Radar


class SABBaseData(object):
    """Decode SA/SB/CB/SC2.0 base data into a lightweight radial structure."""

    def __init__(self, filename, station_lon=None, station_lat=None, station_alt=None):
        """
        :param filename:  radar basedata filename
        :param station_lon:  radar station longitude //units: degree east
        :param station_lat:  radar station latitude //units:degree north
        :param station_alt:  radar station altitude //units: meters
        """
        super(SABBaseData, self).__init__()
        self.filename = filename
        self.station_lon = station_lon
        self.station_lat = station_lat
        self.station_alt = station_alt
        self.fid = _prepare_for_read(self.filename)
        self._raw_buf = _read_all(self.fid, "SAB radial payload")
        self.fid.close()
        self.RadialNum, self.nrays = self._determine_radial_record_layout()
        self.radial, self._status, self._azimuth, self._elevation, self._julian_date, self._msends, \
            self._nyquist, self._unambiguous_range = self._parse_radial()
        self._scan_time = None
        self.sweep_start_ray_index = np.where((self._status == 0) | (self._status == 3))[0]
        self.sweep_end_ray_index = np.where((self._status == 2) | (self._status == 4))[0]
        self.nsweeps = len(self.sweep_start_ray_index)
        self._raw_buf = None

    def _determine_radial_record_layout(self):
        """Determine the radial record size and radar family from the raw buffer."""
        if len(self._raw_buf) < 16:
            raise ValueError("SAB file header is truncated.")
        if self._raw_buf[14:16] != b'\x01\x00':
            raise ValueError("File is not a valid SA/SB/CB file.")
        data_len = len(self._raw_buf)
        if (data_len % 2432 != 0) and (data_len % 4132 != 0) and (data_len % 3132 != 0):
            raise ValueError("SAB file size is inconsistent with the supported radial layouts.")
        # Distinguish SA/SB, CB, and SC layouts by fixed radial record size.
        if data_len % 2432 == 0:
            RadialNum = 2432
            self.Type = "SAB"
        elif data_len % 4132 == 0:
            RadialNum = 4132
            self.Type = 'CB'
        else:
            RadialNum = 3132
            self.Type = 'SC'
        return RadialNum, int(data_len / RadialNum)

    def _RadialNum_SAB_CB(self):
        """Backward-compatible alias for ``_determine_radial_record_layout``."""
        return self._determine_radial_record_layout()

    def _parse_radial(self):
        """Parse all radial records from the in-memory SAB buffer."""
        radial = []
        status = []
        azimuth = []
        elevation = []
        julian_date = []
        msends = []
        nyquist = []
        unambiguous_range = []
        buf = self._raw_buf
        for pos in range(0, len(buf), self.RadialNum):
            iray = self._parse_radial_single(buf[pos:pos + self.RadialNum])
            radial.append(iray)
            status.append(iray['RadialStatus'])
            azimuth.append(iray['AZ'] / 8. * 180. / 4096.)
            elevation.append(iray['El'] / 8. * 180. / 4096.)
            julian_date.append(iray['JulianDate'])
            msends.append(iray['mSends'])
            nyquist.append(iray['Nyquist'] / 100.)
            unambiguous_range.append(iray['URange'] / 10.)
        return radial, np.asarray(status), np.asarray(azimuth), np.asarray(elevation), \
            np.asarray(julian_date), np.asarray(msends), np.asarray(nyquist), np.asarray(unambiguous_range)

    def _parse_radial_single(self, radial_buf):
        Radial = {}
        RadialHeader, size_tmp = _unpack_from_buf(radial_buf, 0, dtype_sab.RadialHeader())
        Radial.update(RadialHeader)
        ref_offset, ref_count = _validate_offset_length(
            len(radial_buf),
            int(RadialHeader['PtrOfReflectivity']) + dtype_sab.InfSize,
            int(RadialHeader['GatesNumberOfReflectivity']),
            "SAB reflectivity block",
        )
        dop_offset, dop_count = _validate_offset_length(
            len(radial_buf),
            int(RadialHeader['PtrOfVelocity']) + dtype_sab.InfSize,
            int(RadialHeader['GatesNumberOfDoppler']),
            "SAB velocity block",
        )
        sw_offset, sw_count = _validate_offset_length(
            len(radial_buf),
            int(RadialHeader['PtrOfSpectrumWidth']) + dtype_sab.InfSize,
            int(RadialHeader['GatesNumberOfDoppler']),
            "SAB spectrum-width block",
        )
        dBZ = np.frombuffer(
            radial_buf,
            dtype=np.uint8,
            count=ref_count,
            offset=ref_offset,
        )
        V = np.frombuffer(
            radial_buf,
            dtype=np.uint8,
            count=dop_count,
            offset=dop_offset,
        )
        W = np.frombuffer(
            radial_buf,
            dtype=np.uint8,
            count=sw_count,
            offset=sw_offset,
        )
        Radial['fields'] = {}
        Radial['fields']['dBZ'] = self._decode_field(dBZ, -32.0)
        Radial['fields']['V'] = self._decode_field(V, -63.5)
        Radial['fields']['W'] = self._decode_field(W, -63.5)
        return Radial

    @staticmethod
    def _decode_field(raw, base):
        out = np.full(raw.shape, np.nan, dtype=np.float32)
        valid = raw > 1
        if np.any(valid):
            out[valid] = (raw[valid].astype(np.float32) - 2.0) / 2.0 + base
        return out

    def get_nyquist_velocity(self):
        """Return the per-ray Nyquist velocity."""
        return self._nyquist

    def get_unambiguous_range(self):
        """Return the per-ray unambiguous range in kilometers."""
        return self._unambiguous_range

    def get_scan_time(self):
        """Return the acquisition time for each ray."""
        if self._scan_time is None:
            self._scan_time = np.array(
                [julian2date(julian, ms) for julian, ms in zip(self._julian_date, self._msends)]
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
        lat, lon, alt, frequency = get_radar_info(self.filename)
        if self.station_lon is not None:
            lon = self.station_lon
        if self.station_lat is not None:
            lat = self.station_lat
        if self.station_alt is not None:
            alt = self.station_alt
        return lat, lon, alt, frequency

    def get_scan_type(self):
        """Return the scan type string."""
        return "ppi"

    def get_sitename(self):
        return get_radar_sitename(self.filename)


class SAB2NRadar(object):
    """Bridge from raw SAB data to an NRadar object."""

    def __init__(self, SAB):
        self.SAB = SAB
        self._dbz_index_cache = {}
        self.v_index_alone = self.get_v_idx()
        self.dBZ_index_alone = self.get_dbz_idx()
        self.dBZ_Res = self.SAB.radial[0]["GateSizeOfReflectivity"]
        sab_elevation = self.SAB.get_elevation()
        for index_with_dbz, index_with_v in zip(self.dBZ_index_alone, self.v_index_alone):
            assert abs(sab_elevation[index_with_v] - sab_elevation[index_with_dbz]) < 0.5, "warning! maybe it is a problem."
            self.interp_dBZ(index_with_dbz, index_with_v)
        keep_mask = np.ones(self.SAB.nrays, dtype=bool)
        keep_mask[self.get_remove_radial_indices()] = False
        self.radial = [iray for iray, keep in zip(self.SAB.radial, keep_mask) if keep]
        self._azimuth = self.SAB.get_azimuth()[keep_mask]
        self._elevation = self.SAB.get_elevation()[keep_mask]
        self._julian_date = self.SAB._julian_date[keep_mask]
        self._msends = self.SAB._msends[keep_mask]
        self._nyquist = self.SAB.get_nyquist_velocity()[keep_mask]
        self._unambiguous_range = self.SAB.get_unambiguous_range()[keep_mask]
        self._scan_time = None
        self.nrays = len(self.radial)
        self.nsweeps = self.SAB.nsweeps - self.dBZ_index_alone.size
        status = self.SAB._status[keep_mask]
        self.sweep_start_ray_index = np.where((status == 0) | (status == 3))[0]
        self.sweep_end_ray_index = np.where((status == 2) | (status == 4))[0]
        self.scan_type = self.SAB.get_scan_type()
        self.latitude, self.longitude, self.altitude, self.frequency = \
            self.SAB.get_latitude_longitude_altitude_frequency()
        self.bins_per_sweep = self.get_nbins_per_sweep()
        self.max_bins = self.bins_per_sweep.max()
        self.range = self.get_range_per_radial(self.max_bins)
        self.azimuth = self.get_azimuth()
        self.elevation = self.get_elevation()
        self.extended_fields = self._build_extended_fields()
        self.fields = self._get_fields()
        self.sitename = self.SAB.get_sitename()

    def get_remove_radial_indices(self):
        """Return the radial indices removed after sweep alignment."""
        index_romove = []
        for isweep in self.dBZ_index_alone:
            index_romove.extend(range(self.SAB.sweep_start_ray_index[isweep], \
                                      self.SAB.sweep_end_ray_index[isweep] + 1))
        return index_romove

    def get_reomve_radial_num(self):
        """Backward-compatible alias for ``get_remove_radial_indices``."""
        return self.get_remove_radial_indices()

    def get_v_idx(self):
        """Return sweep indices that contain Doppler moments without reflectivity."""
        flag = np.array([((self.SAB.radial[idx]['fields']["V"].size != 0) and
                          (self.SAB.radial[idx]['fields']["dBZ"].size == 0)) \
                         for idx in self.SAB.sweep_start_ray_index])
        return np.where(flag == 1)[0]

    def get_dbz_idx(self):
        """Return sweep indices that contain reflectivity without Doppler moments."""
        flag = np.array([((self.SAB.radial[idx]['fields']["V"].size == 0) and
                          (self.SAB.radial[idx]['fields']["dBZ"].size != 0)) \
                         for idx in self.SAB.sweep_start_ray_index])
        return np.where(flag == 1)[0]

    def interp_dBZ(self, field_with_dBZ_num, field_without_dBZ_num):
        """
        Copy reflectivity rays onto the matching Doppler-only sweep by azimuth.
        :param field_with_dBZ_num: source sweep index, 0-based.
        :param field_without_dBZ_num: target sweep index, 0-based.
        """
        azimuth = self.SAB.get_azimuth()
        assert (field_with_dBZ_num + 1) == field_without_dBZ_num, "check interp sweep!"
        dbz_az = azimuth[self.SAB.sweep_start_ray_index[field_with_dBZ_num]: \
                         self.SAB.sweep_end_ray_index[field_with_dBZ_num] + 1]
        v_az = azimuth[self.SAB.sweep_start_ray_index[field_without_dBZ_num]: \
                       self.SAB.sweep_end_ray_index[field_without_dBZ_num] + 1]
        dbz_idx = np.argmin(np.abs(dbz_az.reshape(-1, 1) - v_az.reshape(1, -1)), axis=0) + \
                  self.SAB.sweep_start_ray_index[field_with_dBZ_num]
        v_idx = np.arange(self.SAB.sweep_start_ray_index[field_without_dBZ_num], \
                          self.SAB.sweep_end_ray_index[field_without_dBZ_num] + 1)
        for ind_dbz, ind_v in zip(dbz_idx, v_idx):
            self.SAB.radial[ind_v]["fields"]['dBZ'] = self.SAB.radial[ind_dbz]["fields"]['dBZ']

    def get_azimuth(self):
        """Return the azimuth angle for each retained ray."""
        return self._azimuth

    def get_elevation(self):
        """Return the elevation angle for each retained ray."""
        return self._elevation

    def get_rays_per_sweep(self):
        """Return the number of rays in each retained sweep."""
        return self.sweep_end_ray_index - self.sweep_start_ray_index + 1

    def get_scan_time(self):
        """Return the acquisition time for each retained ray."""
        if self._scan_time is None:
            self._scan_time = np.array(
                [julian2date(julian, ms) for julian, ms in zip(self._julian_date, self._msends)]
            )
        return self._scan_time

    def get_nyquist_velocity(self):
        """Return the per-ray Nyquist velocity."""
        return self._nyquist

    def get_unambiguous_range(self):
        """Return the per-ray unambiguous range."""
        return self._unambiguous_range

    def get_sweep_end_ray_index(self):
        """Return the inclusive end index of each retained sweep."""
        return self.sweep_end_ray_index

    def get_sweep_start_ray_index(self):
        """Return the start index of each retained sweep."""
        return self.sweep_start_ray_index

    def get_nbins_per_sweep(self):
        """Return the Doppler gate count for each retained sweep."""
        return np.array([self.radial[idx]['fields']['V'].size for idx in self.sweep_start_ray_index])

    def get_range_per_radial(self, length):
        """Return Doppler gate-center ranges for a radial of ``length`` bins."""
        Resolution = self.radial[0]["GateSizeOfDoppler"]
        return np.linspace(Resolution, Resolution * length, length)

    def get_dbz_range_per_radial(self, length):
        """Return reflectivity gate-center ranges for a radial of ``length`` bins."""
        Resolution = self.dBZ_Res
        start_range = self.radial[0]["GateSizeOfDoppler"]
        return np.linspace(start_range, start_range + Resolution * (length - 1), length)

    def _get_fields(self):
        """Assemble the retained fields into dense 2-D arrays."""
        field_keys = tuple(self.radial[0]['fields'].keys())
        fields = {ikey: np.full((self.nrays, self.max_bins), np.nan, dtype=np.float64) for ikey in field_keys}
        for iray, radial in enumerate(self.radial):
            radial_fields = radial['fields']
            for ikey in field_keys:
                fields[ikey][iray, :] = self._add_or_del_field(radial_fields.get(ikey), ikey)
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
                "azimuth": np.asarray(self._azimuth[start:end + 1], dtype=np.float32),
                "elevation": np.asarray(self._elevation[start:end + 1], dtype=np.float32),
                "aligned_bins": int(aligned_bins),
            }

        return {"dBZ": extended_sweeps} if extended_sweeps else {}

    def _get_dbz_resample_index(self, source_length):
        cached = self._dbz_index_cache.get(source_length)
        if cached is not None:
            return cached

        source = self.get_dbz_range_per_radial(source_length)
        target = self.range
        right = np.searchsorted(source, target, side="left")
        left = np.clip(right - 1, 0, source_length - 1)
        right = np.clip(right, 0, source_length - 1)
        choose_right = np.abs(source[right] - target) < np.abs(target - source[left])
        nearest = np.where(choose_right, right, left)
        valid = (target >= source[0]) & (target <= source[-1])
        self._dbz_index_cache[source_length] = (valid, nearest)
        return valid, nearest

    def _add_or_del_field(self, dat_ray, key):
        """
        Normalize a radial field to the common Doppler-aligned range grid.
        :param dat_ray: radial field data.
        :param key: field name.
        """
        length = self.max_bins
        if key == "dBZ":
            valid, nearest = self._get_dbz_resample_index(dat_ray.size)
            out = np.full((length,), np.nan)
            out[valid] = dat_ray[nearest[valid]]
            return out
        if dat_ray.size >= length:
            return (dat_ray[:length]).ravel()
        else:
            out = np.full((length,), np.nan)
            out[:dat_ray.size] = dat_ray
            return out.ravel()

    def get_nradar_nyquist_speed(self):
        """array shape (nsweeps)"""
        return np.array([self.radial[idx]['Nyquist'] / 100. for idx in self.sweep_start_ray_index])

    def get_NRadar_nyquist_speed(self):
        """Backward-compatible alias for ``get_nradar_nyquist_speed``."""
        return self.get_nradar_nyquist_speed()

    def get_nradar_unambiguous_range(self):
        """array shape (nsweeps)"""
        return np.array([self.radial[idx]['URange'] / 10. for idx in self.sweep_start_ray_index])

    def get_NRadar_unambiguous_range(self):
        """Backward-compatible alias for ``get_nradar_unambiguous_range``."""
        return self.get_nradar_unambiguous_range()

    def get_fixed_angle(self):
        if self.nsweeps == 9:
            fixed_angle = np.array([0.50, 1.45, 2.40, 3.35, 4.30, 6.00, 9.00, 14.6, 19.5])
        elif self.nsweeps == 14:
            fixed_angle = np.array([0.50, 1.45, 2.40, 3.35, 4.30, 5.25, 6.2, 7.5, 8.7, 10, 12, 14, 16.7, 19.5])
        elif self.nsweeps == 6:
            fixed_angle = np.array([0.50, 1.50, 2.50, 2.50, 3.50, 4.50])
        elif self.nsweeps == 4:
            fixed_angle = np.array([0.50, 2.50, 3.50, 4.50])
        else:
            fixed_angle = np.array([self.radial[idx]['El'] / 8. * 180. / 4096. for idx in self.sweep_start_ray_index])
        return fixed_angle

    def to_prd(self, effective_earth_radius=None):
        """Build an ``NRadar.PRD`` object from the decoded SAB volume."""
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
                          metadata={"original_container": "CINRAD/SAB", "radar_name": "CINRAD/SA/SB/CB/SC"})

    def ToPRD(self, effective_earth_radius=None):
        """Backward-compatible alias for ``to_prd``."""
        return self.to_prd(effective_earth_radius=effective_earth_radius)

    def to_pyart_radar(self, effective_earth_radius=None, **kwargs):
        """Export the decoded SAB volume through the PRD Py-ART adapter."""
        return self.to_prd(effective_earth_radius=effective_earth_radius).to_pyart_radar(**kwargs)

    def ToPyartRadar(self, effective_earth_radius=None, **kwargs):
        """Backward-compatible alias for ``to_pyart_radar``."""
        return self.to_pyart_radar(effective_earth_radius=effective_earth_radius, **kwargs)

    def _get_instrument_parameters(self):
        """ Return a dictionary containing instrument parameters. """

        # pulse width
        pulse_width = get_metadata('pulse_width')
        pulse_width['data'] = np.array([self.radial[0]["GateSizeOfDoppler"] / _LIGHT_SPEED,], dtype='float32')  # m->sec

        # assume that the parameters in the first ray represent the beam widths,
        # bandwidth and frequency in the entire volume

        wavelength_hz = self.frequency * 10 ** 9

        # radar_beam_width_h
        radar_beam_width_h = get_metadata('radar_beam_width_h')
        radar_beam_width_h['data'] = np.array([1, ], dtype='float32')

        # radar_beam_width_v
        radar_beam_width_v = get_metadata('radar_beam_width_v')
        radar_beam_width_v['data'] = np.array([1, ], dtype='float32')

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
