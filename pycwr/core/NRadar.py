# -*- coding: utf-8 -*-
"""
The PRD object adapts to Chinese radar volumes where sweep ranges differ
between elevation angles and low-elevation Doppler and reflectivity scans
may be collected separately.
"""
import json
import numpy as np
import xarray as xr
import pyproj
from ..configure.default_config import DEFAULT_METADATA, CINRAD_field_mapping
from ..core.RadarProduct import PRODUCT_REFERENCE_NOTES, derive_et, derive_vil
from ..core.transforms import  antenna_to_cartesian_cwr, cartesian_to_geographic_aeqd,\
    antenna_vectors_to_cartesian_cwr, antenna_vectors_to_cartesian_rhi, cartesian_to_antenna_cwr,\
    antenna_vectors_to_cartesian_vcs, geographic_to_cartesian_aeqd, resolve_effective_earth_radius
from .interop import build_xradar_sweep_datasets, export_pyart_radar, export_xradar_tree
try:
    from .RadarGridC import get_CR_xy, get_CAPPI_xy, get_CAPPI_3d, get_mosaic_CAPPI_3d
except ImportError:
    from .RadarGrid import get_CR_xy, get_CAPPI_xy, get_CAPPI_3d, get_mosaic_CAPPI_3d

class PRD(object):
    """
    Polarimetry Radar Data (PRD)
    A class for storing antenna coordinate radar data.
    Attributes
    ----------
    fields : dict
        Moment fields. with different variables
    scan_type : str
        Type of scan, one of 'ppi', 'rhi', 'sector' or 'other'. If the scan
        volume contains multiple sweep modes this should be 'other'.
    time : datetime object
        Time at the center of each ray.
    range : numpy array //m
        Range to the center of each gate (bin).
    latitude : scalar//units:degree
        Latitude of the instrument.
    longitude: scalar//units:degree
        Longitude of the instrument.
    altitude : scalar//units:m
        Altitude of the instrument, above sea level.
    fixed_angle : (nsweeps) units:degree
        Target angle for thr sweep. Azimuth angle in RHI modes, elevation
        angle in all other modes.
    azimuth : (nrays) units :degree
        Azimuth of antenna, relative to true North. Azimuth angles are
        recommended to be expressed in the range of [0, 360], but other
        representations are not forbidden.
    elevation : (nrays) units :degree
        Elevation of antenna, relative to the horizontal plane. Elevation
        angles are recommended to be expressed in the range of [-180, 180],
        but other representations are not forbidden.
    sweep_start_ray_index : numpy array(nsweeps)
        Index of the first ray in each sweep relative to the start of the
        volume, 0-based.
    sweep_end_ray_index : numpy array(nsweeps)
        Index of the last ray in each sweep relative to the start of the
        volume, 0-based.
    rays_per_sweep : numpy array (nsweeps)
        Number of rays in each sweep. The data key of this attribute is
        create upon first access from the data in the sweep_start_ray_index and
        sweep_end_ray_index attributes. If the sweep locations needs to be
        modified, do this prior to accessing this attribute or use
        :py:func:`init_rays_per_sweep` to reset the attribute.
    bins_per_sweep : numpy array (nsweeps)    !!!##added
        Number of bins in each sweep. The data key of this attribute is
        create upon first access from the data in the sweep_start_ray_index and
        sweep_end_ray_index attributes. If the sweep locations needs to be
        modified, do this prior to accessing this attribute or use
        :py:func:`init_rays_per_sweep` to reset the attribute.
    nyquist_velocity: numpy array (nsweeps) (m/s)
    unambiguous_range:numpy array (nsweeps) (m/s)
    frequency: constant (GHZ)
    extended_fields : dict
        Optional sidecar storage for native-range fields that are longer than
        the aligned sweep range kept in ``fields``.
    nrays : int
        Number of rays in the volume.
    nsweeps : int
        Number of sweep in the volume.

    """

    FIELD_ALIASES = {
        "dBZ": ("dBZ", "Zc"),
        "Zc": ("Zc", "dBZ"),
        "V": ("V", "Vc"),
        "Vc": ("Vc", "V"),
        "W": ("W", "Wc"),
        "Wc": ("Wc", "W"),
        "ZDR": ("ZDR", "ZDRc"),
        "ZDRc": ("ZDRc", "ZDR"),
        "KDP": ("KDP", "KDPc"),
        "KDPc": ("KDPc", "KDP"),
    }

    def __init__(self, fields,  scan_type, time, range, azimuth, elevation,latitude,
                 longitude, altitude, sweep_start_ray_index, sweep_end_ray_index,
                 fixed_angle, bins_per_sweep, nyquist_velocity, frequency, unambiguous_range,
                 nrays, nsweeps, sitename, pyart_radar=None, effective_earth_radius=None,
                 extended_fields=None, metadata=None):
        super(PRD, self).__init__()
        self.effective_earth_radius = resolve_effective_earth_radius(effective_earth_radius)
        self.extended_fields = extended_fields or {}
        scan_type_value = str(np.asarray(scan_type).item()) if np.asarray(scan_type).shape == () else str(scan_type)
        if scan_type_value == "ppi":
            sweep_mode = np.full(nsweeps, "azimuth_surveillance", dtype=object)
        elif scan_type_value == "rhi":
            sweep_mode = np.full(nsweeps, "rhi", dtype=object)
        else:
            sweep_mode = np.full(nsweeps, "sector", dtype=object)
        follow_mode = np.full(nsweeps, "none", dtype=object)
        prt_mode = np.full(nsweeps, "fixed", dtype=object)
        sweep_group_name = np.array(["sweep_%d" % i for i in np.arange(nsweeps, dtype=int)], dtype=object)
        time_values = np.asarray(time)
        if np.issubdtype(time_values.dtype, np.datetime64):
            ray_times_increase = "true" if np.all(np.diff(time_values.astype("datetime64[ns]").astype(np.int64)) >= 0) else "false"
        else:
            ray_times_increase = "true"
        metadata_defaults = {
            "Conventions": "Cf/Radial",
            "version": "Cf/Radial2.1",
            "title": "%s radar volume" % sitename,
            "instrument_name": sitename,
            "institution": "unknown",
            "references": "xradar Cf/Radial2.1 compatible export from pycwr",
            "source": "pycwr.PRD",
            "history": "created by pycwr",
            "comment": "",
            "platform_is_mobile": "false",
            "site_name": sitename,
            "scan_name": scan_type_value,
            "scan_id": 0,
            "ray_times_increase": ray_times_increase,
            "simulated": "false",
            "platform_type": "fixed",
            "instrument_type": "radar",
            "primary_axis": "axis_z",
            "volume_number": 0,
        }
        if metadata is not None:
            metadata_defaults.update(metadata)
        self.metadata = metadata_defaults
        keys = fields.keys()
        self.fields = []
        for idx, (istart, iend) in enumerate(zip(sweep_start_ray_index, sweep_end_ray_index)):
            x, y, z = antenna_vectors_to_cartesian_cwr(range[:bins_per_sweep[idx]], azimuth[istart:iend+1],\
                                                   elevation[istart:iend+1], altitude,
                                                   effective_earth_radius=self.effective_earth_radius)
            lon, lat = cartesian_to_geographic_aeqd(x, y, longitude, latitude)
            isweep_data = xr.Dataset(coords={'azimuth': (['time', ], azimuth[istart:iend+1]),
                                            'elevation': (['time',], elevation[istart:iend+1]),
                                             'x':(['time','range'], x),
                                             'y':(['time', 'range'], y),
                                             'z':(['time', 'range'], z),
                                             'lat':(['time','range'], lat),
                                             'lon':(['time','range'], lon),
                                            'range': range[:bins_per_sweep[idx]], 'time': time[istart:iend+1]})
            isweep_data.azimuth.attrs = DEFAULT_METADATA['azimuth']
            isweep_data.elevation.attrs = DEFAULT_METADATA['elevation']
            isweep_data.range.attrs = DEFAULT_METADATA['range']
            isweep_data.time.attrs = DEFAULT_METADATA['time']
            isweep_data.x.attrs = DEFAULT_METADATA['x']
            isweep_data.y.attrs = DEFAULT_METADATA['y']
            isweep_data.z.attrs = DEFAULT_METADATA['z']
            isweep_data.lon.attrs = DEFAULT_METADATA['lon']
            isweep_data.lat.attrs = DEFAULT_METADATA['lat']
            for ikey in keys:
                isweep_data[ikey] = (['time','range'], fields[ikey][istart:iend+1, :bins_per_sweep[idx]])
                isweep_data[ikey].attrs = DEFAULT_METADATA[CINRAD_field_mapping[ikey]]
            self.fields.append(isweep_data)
        self.scan_info = xr.Dataset(data_vars={"latitude":latitude,"longitude":longitude,
                        "altitude":altitude,"scan_type":scan_type,  "frequency":frequency,
                        "start_time":time[0], "end_time":time[-1],
                        "volume_number": np.int32(self.metadata["volume_number"]),
                        "effective_earth_radius": self.effective_earth_radius,
                         "nyquist_velocity":(['sweep',], nyquist_velocity[:nsweeps]),
                        "unambiguous_range":(['sweep',], unambiguous_range[:nsweeps]),
                        "rays_per_sweep": (['sweep',], (sweep_end_ray_index-sweep_start_ray_index+1)[:nsweeps]),
                        "fixed_angle": (["sweep",], fixed_angle[:nsweeps]),
                        "beam_width":(["sweep",], 360./(sweep_end_ray_index-sweep_start_ray_index+1)[:nsweeps]),
                        "sweep_mode": (["sweep"], sweep_mode),
                        "follow_mode": (["sweep"], follow_mode),
                        "prt_mode": (["sweep"], prt_mode),
                        "sweep_group_name": (["sweep"], sweep_group_name)},
                        coords={"sweep": np.arange(nsweeps, dtype=int)})
        self.scan_info['latitude'].attrs = DEFAULT_METADATA['latitude']
        self.scan_info['longitude'].attrs = DEFAULT_METADATA['longitude']
        self.scan_info['altitude'].attrs = DEFAULT_METADATA['altitude']
        self.scan_info['scan_type'].attrs = DEFAULT_METADATA['scan_type']
        self.scan_info['frequency'].attrs = DEFAULT_METADATA['frequency']
        self.scan_info['effective_earth_radius'].attrs = {
            'units': 'meters',
            'standard_name': 'effective_earth_radius',
            'long_name': 'effective earth radius used by beam geometry',
        }
        self.scan_info['volume_number'].attrs = {'long_name': 'volume number'}
        self.scan_info['nyquist_velocity'].attrs = DEFAULT_METADATA['nyquist_velocity']
        self.scan_info['unambiguous_range'].attrs = DEFAULT_METADATA['unambiguous_range']
        self.scan_info['rays_per_sweep'].attrs = DEFAULT_METADATA['rays_per_sweep']
        self.scan_info['fixed_angle'].attrs = DEFAULT_METADATA['fixed_angle']
        self.scan_info['start_time'].attrs = DEFAULT_METADATA['start_time']
        self.scan_info['end_time'].attrs = DEFAULT_METADATA['end_time']
        self.scan_info['beam_width'].attrs = DEFAULT_METADATA['beam_width']
        self.scan_info['sweep_mode'].attrs = {'long_name': 'sweep mode'}
        self.scan_info['follow_mode'].attrs = {'long_name': 'follow mode'}
        self.scan_info['prt_mode'].attrs = {'long_name': 'PRT mode'}
        self.scan_info['sweep_group_name'].attrs = {'long_name': 'sweep group name'}
        self.nsweeps = nsweeps
        self.nrays = nrays
        self.sitename = sitename
        self._fixed_angle_sort_index = self.scan_info["fixed_angle"].argsort().values
        self._ordered_az_cache = None
        self._vol_cache = {}
        self._pyart_cache = {}
        self._xradar_cache = {}
        self._native_field_cache = {}
        self._sweep_field_cache = {}
        self._sorted_sweep_cache = {}
        self._summary_cache = None
        self._site_projection = None
        self.get_vol_data()
        self.product = xr.Dataset()
        self.PyartRadar = pyart_radar

    @staticmethod
    def _field_prefers_native_range(field_name):
        """Return whether a field should default to the native reflectivity range."""
        return str(field_name) in {"dBZ", "Zc"}

    @classmethod
    def field_alias_candidates(cls, field_name):
        """Return the preferred aliases for one logical radar field name."""
        name = str(field_name)
        return cls.FIELD_ALIASES.get(name, (name,))

    def resolve_field_name(self, field_name, sweep=None, range_mode=None, required=True):
        """Resolve a logical field name to an available field on one sweep or volume."""
        if sweep is None:
            available = set(self.available_fields(range_mode=range_mode))
        else:
            available = set(self.available_fields(sweep=int(sweep), range_mode=range_mode))
        for candidate in self.field_alias_candidates(field_name):
            if candidate in available:
                return candidate
        if required:
            raise KeyError(field_name)
        return None

    def _resolve_field_range_mode(self, field_name, range_mode=None):
        """Resolve the effective range mode for a single field request."""
        if range_mode is None:
            return "native" if self._field_prefers_native_range(field_name) else "aligned"
        self._validate_range_mode(range_mode)
        return range_mode

    def _resolve_export_range_mode(self, field_names=None, range_mode=None):
        """Resolve the effective range mode for mixed-field export requests."""
        if range_mode is not None:
            self._validate_range_mode(range_mode)
            return range_mode
        if field_names:
            normalized = [str(name) for name in field_names]
            if normalized and all(self._field_prefers_native_range(name) for name in normalized):
                return "native"
        return "aligned"

    def to_pyart_radar(
        self,
        range_mode=None,
        field_names=None,
        use_external=None,
        strict=False,
        force_rebuild=False,
    ):
        """Export the PRD volume as a Py-ART Radar object."""
        range_mode = self._resolve_export_range_mode(field_names=field_names, range_mode=range_mode)
        cache_key = (
            range_mode,
            None if field_names is None else tuple(field_names),
            use_external,
            bool(strict),
        )
        if not force_rebuild and cache_key in self._pyart_cache:
            return self._pyart_cache[cache_key]
        radar = export_pyart_radar(
            self,
            existing_radar=self.PyartRadar,
            range_mode=range_mode,
            field_names=field_names,
            use_external=use_external,
            strict=strict,
        )
        self._pyart_cache[cache_key] = radar
        if range_mode == "aligned" and field_names is None and use_external in (None, False):
            self.PyartRadar = radar
        return radar

    def ToPyartRadar(self, **kwargs):
        """Backward-compatible alias for ``to_pyart_radar``."""
        return self.to_pyart_radar(**kwargs)

    def to_wsr98d(self, filename, **kwargs):
        """Write the current PRD volume to WSR-98D base-data format."""
        from ..io import write_wsr98d

        return write_wsr98d(self, filename, **kwargs)

    def to_nexrad_level2_msg31(self, filename, **kwargs):
        """Write the current PRD volume to NEXRAD Level II MSG31 format."""
        from ..io import write_nexrad_level2_msg31

        return write_nexrad_level2_msg31(self, filename, **kwargs)

    def to_nexrad_level2_msg1(self, filename, **kwargs):
        """Write the current PRD volume to NEXRAD Level II MSG1 format."""
        from ..io import write_nexrad_level2_msg1

        return write_nexrad_level2_msg1(self, filename, **kwargs)

    def apply_dualpol_qc(
        self,
        sweeps=None,
        inplace=False,
        band="C",
        use_existing_kdp=True,
        clear_air_mode="label",
        clear_air_max_ref=15.0,
        clear_air_max_rhohv=0.97,
        clear_air_max_phidp_texture=10.0,
        clear_air_max_snr=20.0,
    ):
        """
        Apply dual-polarization quality control to this PRD volume.
        """
        from ..qc.pipeline import apply_dualpol_qc

        return apply_dualpol_qc(
            self,
            sweeps=sweeps,
            inplace=inplace,
            band=band,
            use_existing_kdp=use_existing_kdp,
            clear_air_mode=clear_air_mode,
            clear_air_max_ref=clear_air_max_ref,
            clear_air_max_rhohv=clear_air_max_rhohv,
            clear_air_max_phidp_texture=clear_air_max_phidp_texture,
            clear_air_max_snr=clear_air_max_snr,
        )

    def retrieve_vad(self, sweeps=None, field_name=None, range_mode="aligned", **kwargs):
        """
        Retrieve horizontal wind ring-by-ring from one or more PPI sweeps.
        """
        from ..retrieve.WindField import retrieve_vad

        return retrieve_vad(
            self,
            sweeps=sweeps,
            field_name=field_name,
            range_mode=range_mode,
            **kwargs
        )

    def retrieve_vvp(self, sweep, field_name=None, range_mode="aligned", **kwargs):
        """
        Retrieve a local horizontal wind field from one PPI sweep.
        """
        from ..retrieve.WindField import retrieve_vvp

        return retrieve_vvp(
            self,
            sweep=sweep,
            field_name=field_name,
            range_mode=range_mode,
            **kwargs
        )

    def retrieve_vwp(self, sweeps=None, field_name=None, range_mode="aligned", **kwargs):
        """
        Retrieve a robust vertical wind profile aggregated from VAD retrievals.
        """
        from ..retrieve.WindField import retrieve_vwp

        return retrieve_vwp(
            self,
            sweeps=sweeps,
            field_name=field_name,
            range_mode=range_mode,
            **kwargs
        )

    def classify_hydrometeors(self, sweeps=None, inplace=False, **kwargs):
        """
        Add hydrometeor classification fields to selected sweeps.
        """
        from ..retrieve import apply_hydrometeor_classification

        return apply_hydrometeor_classification(
            self,
            sweeps=sweeps,
            inplace=inplace,
            **kwargs
        )

    def add_hydrometeor_classification(self, sweeps=None, **kwargs):
        """
        In-place convenience wrapper for ``classify_hydrometeors``.
        """
        return self.classify_hydrometeors(sweeps=sweeps, inplace=True, **kwargs)

    def to_xradar_sweeps(self, range_mode=None, field_names=None, force_rebuild=False):
        """Export the PRD volume as sweep-level xarray datasets."""
        range_mode = self._resolve_export_range_mode(field_names=field_names, range_mode=range_mode)
        cache_key = ("datasets", range_mode, None if field_names is None else tuple(field_names))
        if not force_rebuild and cache_key in self._xradar_cache:
            return self._xradar_cache[cache_key]
        datasets = build_xradar_sweep_datasets(self, range_mode=range_mode, field_names=field_names)
        self._xradar_cache[cache_key] = datasets
        return datasets

    def to_xradar(self, range_mode=None, field_names=None, strict=True, force_rebuild=False):
        """Export the PRD volume as a DataTree for xradar-style workflows."""
        range_mode = self._resolve_export_range_mode(field_names=field_names, range_mode=range_mode)
        cache_key = ("tree", range_mode, None if field_names is None else tuple(field_names), bool(strict))
        if not force_rebuild and cache_key in self._xradar_cache:
            return self._xradar_cache[cache_key]
        tree = export_xradar_tree(self, range_mode=range_mode, field_names=field_names, strict=strict)
        self._xradar_cache[cache_key] = tree
        return tree

    def has_extended_field(self, sweep, field_name):
        """Return whether a sweep exposes native-range sidecar data for a field."""
        return int(sweep) in self.extended_fields.get(field_name, {})

    def available_fields(self, sweep=None, range_mode=None):
        """Return the available field names for one sweep or the full volume."""
        range_mode = self._resolve_export_range_mode(range_mode=range_mode)
        if sweep is None:
            names = []
            seen = set()
            for sweep_index in range(int(self.nsweeps)):
                for name in self.available_fields(sweep=sweep_index, range_mode=range_mode):
                    if name not in seen:
                        seen.add(name)
                        names.append(name)
            return names
        sweep = int(sweep)
        names = list(self.fields[sweep].data_vars)
        if range_mode == "native":
            for field_name, sweep_map in self.extended_fields.items():
                if sweep in sweep_map and field_name not in names:
                    names.append(field_name)
        return names

    def sweep_summary(self):
        """Return per-sweep summary records for quick inspection and docs."""
        rows = []
        fixed_angles = np.asarray(self.scan_info["fixed_angle"].values, dtype=np.float64)
        rays_per_sweep = np.asarray(self.scan_info["rays_per_sweep"].values, dtype=np.int32)
        for sweep in range(int(self.nsweeps)):
            aligned_fields = list(self.fields[sweep].data_vars)
            native_fields = [
                field_name
                for field_name, sweep_map in self.extended_fields.items()
                if sweep in sweep_map
            ]
            range_values = np.asarray(self.fields[sweep]["range"].values, dtype=np.float64)
            rows.append(
                {
                    "sweep": sweep,
                    "fixed_angle": float(fixed_angles[sweep]),
                    "rays": int(rays_per_sweep[sweep]),
                    "aligned_fields": aligned_fields,
                    "native_fields": native_fields,
                    "aligned_max_range_m": float(range_values[-1]) if range_values.size else 0.0,
                    "native_max_range_m_by_field": {
                        field_name: float(np.asarray(sweep_map[sweep]["range"], dtype=np.float64)[-1])
                        for field_name, sweep_map in self.extended_fields.items()
                        if sweep in sweep_map and np.asarray(sweep_map[sweep]["range"]).size
                    },
                }
            )
        return rows

    def summary(self):
        """Return a lightweight volume summary for users and tests."""
        if self._summary_cache is None:
            self._summary_cache = {
                "sitename": self.sitename,
                "scan_type": str(np.asarray(self.scan_info["scan_type"].values).item()),
                "nsweeps": int(self.nsweeps),
                "nrays": int(self.nrays),
                "fields": self.available_fields(),
                "sweeps": self.sweep_summary(),
                "site": {
                    "longitude": float(np.asarray(self.scan_info["longitude"].values)),
                    "latitude": float(np.asarray(self.scan_info["latitude"].values)),
                    "altitude": float(np.asarray(self.scan_info["altitude"].values)),
                },
            }
        return self._summary_cache

    def _site_longitude(self):
        return float(np.asarray(self.scan_info["longitude"].values))

    def _site_latitude(self):
        return float(np.asarray(self.scan_info["latitude"].values))

    def _site_altitude(self):
        return float(np.asarray(self.scan_info["altitude"].values))

    def _get_local_projection(self):
        if self._site_projection is None:
            self._site_projection = pyproj.Proj(
                {"proj": "aeqd", "lon_0": self._site_longitude(), "lat_0": self._site_latitude()}
            )
        return self._site_projection

    @staticmethod
    def _validate_range_mode(range_mode):
        if range_mode not in ("aligned", "native"):
            raise ValueError("range_mode must be 'aligned' or 'native'")

    def get_native_sweep_field(self, sweep, field_name):
        """
        Return a sweep field on its native range grid when sidecar data exists.

        The current ``self.fields`` datasets remain aligned to the Doppler grid.
        This helper exposes the longer native reflectivity range without
        changing those existing arrays.
        """
        sweep = int(sweep)
        sweep_dataset = self.fields[sweep]
        if field_name not in sweep_dataset:
            raise KeyError(field_name)

        native_field = self.extended_fields.get(field_name, {}).get(sweep)
        if native_field is None:
            return sweep_dataset[field_name]
        cache_key = (sweep, field_name)
        cached = self._native_field_cache.get(cache_key)
        if cached is not None:
            return cached

        attrs = dict(sweep_dataset[field_name].attrs)
        radar_altitude = self._site_altitude()
        radar_lon = self._site_longitude()
        radar_lat = self._site_latitude()
        x, y, z = antenna_vectors_to_cartesian_cwr(
            native_field["range"],
            native_field["azimuth"],
            native_field["elevation"],
            radar_altitude,
            effective_earth_radius=self.effective_earth_radius,
        )
        lon, lat = cartesian_to_geographic_aeqd(x, y, radar_lon, radar_lat)
        data = xr.DataArray(
            native_field["data"],
            dims=("time", "range"),
            coords={
                "time": native_field["time"],
                "azimuth": ("time", native_field["azimuth"]),
                "elevation": ("time", native_field["elevation"]),
                "x": (("time", "range"), x),
                "y": (("time", "range"), y),
                "z": (("time", "range"), z),
                "lon": (("time", "range"), lon),
                "lat": (("time", "range"), lat),
                "range": native_field["range"],
            },
            name=field_name,
            attrs=attrs,
        )
        if "time" in sweep_dataset[field_name].dims:
            current_time = np.asarray(sweep_dataset["time"].values)
            if np.array_equal(current_time, native_field["time"]):
                self._native_field_cache[cache_key] = data
                return data
            data = data.sel(time=current_time)
            self._native_field_cache[cache_key] = data
            return data
        if "azimuth" in sweep_dataset[field_name].dims:
            current_azimuth = np.asarray(sweep_dataset["azimuth"].values)
            data = data.swap_dims({"time": "azimuth"}).sel(azimuth=current_azimuth)
            self._native_field_cache[cache_key] = data
            return data
        self._native_field_cache[cache_key] = data
        return data

    def get_sweep_field(self, sweep, field_name, range_mode=None, sort_by_azimuth=False):
        """Return a sweep field in aligned or native-range mode."""
        range_mode = self._resolve_field_range_mode(field_name, range_mode=range_mode)
        sweep = int(sweep)
        cache_key = (sweep, field_name, range_mode, bool(sort_by_azimuth))
        cached = self._sweep_field_cache.get(cache_key)
        if cached is not None:
            return cached
        if range_mode == "native":
            field = self.get_native_sweep_field(sweep, field_name)
        else:
            field = self.fields[sweep][field_name]
        if sort_by_azimuth:
            field = self._sort_field_by_azimuth(field)
            self._sweep_field_cache[cache_key] = field
            return field
        self._sweep_field_cache[cache_key] = field
        return field

    @staticmethod
    def _product_name(base_name, range_mode):
        return base_name if range_mode == "aligned" else "%s_native" % base_name

    @staticmethod
    def _clip_range_window(ranges, values, max_range_km):
        """Clip range/value arrays to the requested maximum range."""
        if max_range_km is None:
            return ranges, values
        max_range_m = float(max_range_km) * 1000.0
        range_mask = ranges <= max_range_m
        if not np.any(range_mask):
            return ranges[:1], values[:, :1]
        return ranges[range_mask], values[:, range_mask]

    def _extract_sweep_volume(self, sweep, field_name, range_mode=None, max_range_km=None):
        """Return sweep azimuth/range/value arrays for Cartesian gridding."""
        range_mode = self._resolve_field_range_mode(field_name, range_mode=range_mode)
        sweep = int(sweep)
        sweep_dataset = self.fields[sweep]
        if field_name not in sweep_dataset:
            raise KeyError(field_name)

        if range_mode == "native":
            native_field = self.extended_fields.get(field_name, {}).get(sweep)
            if native_field is not None:
                azimuth = np.asarray(native_field["azimuth"], dtype=np.float64)
                ranges = np.asarray(native_field["range"], dtype=np.float64)
                values = np.asarray(native_field["data"], dtype=np.float64)
                sort_index = np.argsort(azimuth)
                azimuth = azimuth[sort_index]
                values = values[sort_index, :]
                ranges, values = self._clip_range_window(ranges, values, max_range_km)
                return azimuth, ranges, values
            if field_name == "Zc" and "Zc" in sweep_dataset:
                native_ref = self.extended_fields.get("dBZ", {}).get(sweep)
                if native_ref is not None:
                    azimuth = np.asarray(native_ref["azimuth"], dtype=np.float64)
                    ranges = np.asarray(native_ref["range"], dtype=np.float64)
                    values = np.asarray(native_ref["data"], dtype=np.float64)
                    sort_index = np.argsort(azimuth)
                    azimuth = azimuth[sort_index]
                    values = values[sort_index, :]
                    corrected = np.asarray(sweep_dataset["Zc"].values, dtype=np.float64)[sort_index, :]
                    copy_len = min(corrected.shape[1], values.shape[1])
                    values[:, :copy_len] = np.where(
                        np.isfinite(corrected[:, :copy_len]),
                        corrected[:, :copy_len],
                        values[:, :copy_len],
                    )
                    ranges, values = self._clip_range_window(ranges, values, max_range_km)
                    return azimuth, ranges, values

        azimuth = np.asarray(sweep_dataset["azimuth"].values, dtype=np.float64)
        values = np.asarray(sweep_dataset[field_name].values, dtype=np.float64)
        if values.ndim != 2:
            raise ValueError("Sweep field %s must be 2-D for Cartesian gridding." % field_name)
        sort_index = np.argsort(azimuth)
        azimuth = azimuth[sort_index]
        values = values[sort_index, :]
        ranges = np.asarray(sweep_dataset["range"].values, dtype=np.float64)
        ranges, values = self._clip_range_window(ranges, values, max_range_km)
        return azimuth, ranges, values

    def _invalidate_cached_views(self):
        """Clear cached derived views after the underlying sweep datasets change."""
        self._ordered_az_cache = None
        self._vol_cache.clear()
        self._pyart_cache.clear()
        self.PyartRadar = None
        self._xradar_cache.clear()
        self._native_field_cache.clear()
        self._sweep_field_cache.clear()
        self._sorted_sweep_cache.clear()
        self._summary_cache = None

    @staticmethod
    def _sort_field_by_azimuth(field):
        """Return a field ordered by azimuth using a stable NumPy indexer."""
        if "time" in field.dims:
            azimuth = np.asarray(field["azimuth"].values, dtype=np.float64)
            order = np.argsort(azimuth, kind="mergesort")
            if order.size and not np.array_equal(order, np.arange(order.size, dtype=order.dtype)):
                field = field.isel(time=order)
            return field.swap_dims({"time": "azimuth"})
        if "azimuth" in field.dims:
            azimuth = np.asarray(field["azimuth"].values, dtype=np.float64)
            order = np.argsort(azimuth, kind="mergesort")
            if order.size and not np.array_equal(order, np.arange(order.size, dtype=order.dtype)):
                field = field.isel(azimuth=order)
        return field

    def _get_sorted_sweep_dataset(self, sweep):
        """Return a cached azimuth-sorted sweep dataset."""
        sweep = int(sweep)
        cached = self._sorted_sweep_cache.get(sweep)
        if cached is not None:
            return cached
        sorted_dataset = self._sort_field_by_azimuth(self.fields[sweep])
        self._sorted_sweep_cache[sweep] = sorted_dataset
        return sorted_dataset

    def _sort_sweep_by_azimuth(self, sweep_dataset, sweep_index=None):
        """Return a sweep dataset ordered by azimuth without mutating the input."""
        if sweep_index is not None:
            return self._get_sorted_sweep_dataset(sweep_index)
        return self._sort_field_by_azimuth(sweep_dataset)

    def _build_azimuth_sorted_fields(self):
        """Build azimuth-sorted sweep datasets for the full volume."""
        return [self._get_sorted_sweep_dataset(isweep) for isweep in self.scan_info.sweep.values]

    def ordered_az(self, inplace=False):
        """
        Return or apply a view where azimuth is the sweep dimension.
        """
        sorted_fields = self._build_azimuth_sorted_fields()
        if inplace:
            self.fields = sorted_fields
            self._invalidate_cached_views()
            return None
        if self._ordered_az_cache is None:
            self._ordered_az_cache = AzimuthSortedPRD(
                scan_info=self.scan_info,
                fields=sorted_fields,
                effective_earth_radius=self.effective_earth_radius,
            )
        return self._ordered_az_cache

    def _resolve_section_sample_spacing(self, field_name, range_mode, sample_spacing):
        """Resolve section sampling spacing in meters."""
        if sample_spacing is not None:
            spacing = float(sample_spacing)
            if spacing <= 0.0:
                raise ValueError("sample_spacing must be positive")
            return spacing
        reference_field = self.get_sweep_field(int(self.scan_info.sweep.values[0]), field_name, range_mode=range_mode)
        reference_range = np.asarray(reference_field.range.values, dtype=np.float64)
        if reference_range.size < 2:
            raise ValueError("section extraction requires at least two range gates.")
        spacing = float(reference_range[1] - reference_range[0])
        if spacing <= 0.0:
            raise ValueError("range spacing must be positive to extract sections.")
        return spacing

    @staticmethod
    def _build_section_line(start_point, end_point, sample_spacing):
        """Build Cartesian sample points along a requested section line."""
        start_point = np.asarray(start_point, dtype=np.float64)
        end_point = np.asarray(end_point, dtype=np.float64)
        if start_point.shape != (2,) or end_point.shape != (2,):
            raise ValueError("start_point and end_point must be 2-element coordinates.")
        line_length = float(np.hypot(*(end_point - start_point)))
        if line_length == 0.0:
            raise ValueError("start_point and end_point must not be identical.")
        npoints = max(2, int(np.ceil(line_length / float(sample_spacing))) + 1)
        distance = np.linspace(0.0, line_length, npoints, dtype=np.float64)
        fraction = distance / line_length
        x_line = start_point[0] + fraction * (end_point[0] - start_point[0])
        y_line = start_point[1] + fraction * (end_point[1] - start_point[1])
        return x_line, y_line, distance

    @staticmethod
    def _deduplicate_sorted_azimuth(azimuth, values):
        """Drop duplicate azimuth rows after sorting to keep interpolation stable."""
        if azimuth.size == 0:
            return azimuth, values
        unique_azimuth, unique_index = np.unique(azimuth, return_index=True)
        if unique_azimuth.size == azimuth.size:
            return azimuth, values
        return unique_azimuth, values[unique_index, :]

    @staticmethod
    def _pad_section_rows(rows, fill_value=np.nan):
        """Pad a list of 1-D rows to a common distance length."""
        row_arrays = [np.asarray(row, dtype=np.float64).reshape(-1) for row in rows]
        max_len = max(row.size for row in row_arrays)
        padded = np.full((len(row_arrays), max_len), fill_value, dtype=np.float64)
        for idx, row in enumerate(row_arrays):
            padded[idx, : row.size] = row
        return padded

    @classmethod
    def _interpolate_ppi_strip_arrays(cls, azimuth, ranges, values, target_azimuth, target_range):
        """Linearly interpolate a sweep from NumPy azimuth/range/value arrays."""
        azimuth = np.asarray(azimuth, dtype=np.float64)
        ranges = np.asarray(ranges, dtype=np.float64)
        values = np.asarray(values, dtype=np.float64)
        target_azimuth = np.mod(np.asarray(target_azimuth, dtype=np.float64), 360.0)
        target_range = np.asarray(target_range, dtype=np.float64)

        out = np.full(target_range.shape, np.nan, dtype=np.float64)
        if values.ndim != 2 or azimuth.size < 2 or ranges.size < 2:
            return out

        azimuth, values = cls._deduplicate_sorted_azimuth(azimuth, values)
        if azimuth.size < 2:
            return out

        azimuth = np.mod(azimuth, 360.0)
        azimuth_ext = np.concatenate(([azimuth[-1] - 360.0], azimuth, [azimuth[0] + 360.0]))
        value_ext = np.concatenate((values[-1:, :], values, values[:1, :]), axis=0)

        valid = np.isfinite(target_azimuth) & np.isfinite(target_range)
        valid &= target_range >= ranges[0]
        valid &= target_range <= ranges[-1]
        if not np.any(valid):
            return out

        az_hi = np.searchsorted(azimuth_ext, target_azimuth, side="right")
        az_hi = np.clip(az_hi, 1, azimuth_ext.size - 1)
        az_lo = az_hi - 1

        range_hi = np.searchsorted(ranges, target_range, side="right")
        range_hi = np.clip(range_hi, 1, ranges.size - 1)
        range_lo = range_hi - 1

        az_0 = azimuth_ext[az_lo]
        az_1 = azimuth_ext[az_hi]
        r_0 = ranges[range_lo]
        r_1 = ranges[range_hi]
        span_valid = (az_1 > az_0) & (r_1 > r_0)
        valid &= span_valid
        if not np.any(valid):
            return out

        mat_00 = value_ext[az_lo, range_lo]
        mat_01 = value_ext[az_lo, range_hi]
        mat_10 = value_ext[az_hi, range_lo]
        mat_11 = value_ext[az_hi, range_hi]

        with np.errstate(invalid="ignore", divide="ignore"):
            all_four = valid & np.isfinite(mat_00) & np.isfinite(mat_01) & np.isfinite(mat_10) & np.isfinite(mat_11)
            out[all_four] = (
                mat_00[all_four] * (az_1[all_four] - target_azimuth[all_four]) * (r_1[all_four] - target_range[all_four])
                + mat_10[all_four] * (target_azimuth[all_four] - az_0[all_four]) * (r_1[all_four] - target_range[all_four])
                + mat_01[all_four] * (az_1[all_four] - target_azimuth[all_four]) * (target_range[all_four] - r_0[all_four])
                + mat_11[all_four] * (target_azimuth[all_four] - az_0[all_four]) * (target_range[all_four] - r_0[all_four])
            ) / ((az_1[all_four] - az_0[all_four]) * (r_1[all_four] - r_0[all_four]))

            range_low = valid & np.isnan(out) & np.isfinite(mat_00) & np.isfinite(mat_01)
            out[range_low] = (
                mat_00[range_low] * (r_1[range_low] - target_range[range_low])
                + mat_01[range_low] * (target_range[range_low] - r_0[range_low])
            ) / (r_1[range_low] - r_0[range_low])

            range_high = valid & np.isnan(out) & np.isfinite(mat_10) & np.isfinite(mat_11)
            out[range_high] = (
                mat_10[range_high] * (r_1[range_high] - target_range[range_high])
                + mat_11[range_high] * (target_range[range_high] - r_0[range_high])
            ) / (r_1[range_high] - r_0[range_high])

            az_low = valid & np.isnan(out) & np.isfinite(mat_00) & np.isfinite(mat_10)
            out[az_low] = (
                mat_00[az_low] * (az_1[az_low] - target_azimuth[az_low])
                + mat_10[az_low] * (target_azimuth[az_low] - az_0[az_low])
            ) / (az_1[az_low] - az_0[az_low])

            az_high = valid & np.isnan(out) & np.isfinite(mat_01) & np.isfinite(mat_11)
            out[az_high] = (
                mat_01[az_high] * (az_1[az_high] - target_azimuth[az_high])
                + mat_11[az_high] * (target_azimuth[az_high] - az_0[az_high])
            ) / (az_1[az_high] - az_0[az_high])
        return out

    @classmethod
    def _interpolate_ppi_strip(cls, field, target_azimuth, target_range):
        """Linearly interpolate a sweep along arbitrary azimuth/range targets."""
        return cls._interpolate_ppi_strip_arrays(
            np.asarray(field.azimuth.values, dtype=np.float64),
            np.asarray(field.range.values, dtype=np.float64),
            np.asarray(field.values, dtype=np.float64),
            target_azimuth,
            target_range,
        )

    def _build_section_dataset(
        self,
        field_name,
        distance,
        field_values,
        x_coords,
        y_coords,
        z_coords,
        lon_coords,
        lat_coords,
        azimuth,
        ranges,
        elevation,
        source_sweep,
        section_type,
        range_mode,
        interpolation,
        extra_attrs=None,
    ):
        """Create a standardized xarray dataset for section extraction."""
        distance = np.asarray(distance, dtype=np.float64).reshape(-1)
        field_values = np.asarray(field_values, dtype=np.float64)
        x_coords = np.asarray(x_coords, dtype=np.float64)
        y_coords = np.asarray(y_coords, dtype=np.float64)
        z_coords = np.asarray(z_coords, dtype=np.float64)
        lon_coords = np.asarray(lon_coords, dtype=np.float64)
        lat_coords = np.asarray(lat_coords, dtype=np.float64)
        azimuth = np.asarray(azimuth, dtype=np.float64)
        ranges = np.asarray(ranges, dtype=np.float64)
        elevation = np.asarray(elevation, dtype=np.float64)
        source_sweep = np.asarray(source_sweep, dtype=np.int32)
        if field_values.ndim != 2:
            raise ValueError("field_values must be a 2-D section array.")
        if distance.size != field_values.shape[1]:
            raise ValueError("distance coordinate must match the section distance dimension.")
        distance_ground = np.hypot(x_coords, y_coords)

        section = xr.Dataset(
            data_vars={
                field_name: (("sweep", "distance"), field_values),
                "x": (("sweep", "distance"), x_coords),
                "y": (("sweep", "distance"), y_coords),
                "z": (("sweep", "distance"), z_coords),
                "lon": (("sweep", "distance"), lon_coords),
                "lat": (("sweep", "distance"), lat_coords),
                "azimuth": (("sweep", "distance"), azimuth),
                "range": (("sweep", "distance"), ranges),
                "elevation": (("sweep", "distance"), elevation),
                "distance_ground": (("sweep", "distance"), distance_ground),
                "source_sweep": (("sweep",), source_sweep),
            },
            coords={
                "sweep": np.arange(field_values.shape[0], dtype=np.int32),
                "distance": distance,
            },
            attrs={
                "field_name": field_name,
                "section_type": section_type,
                "scan_type": str(np.asarray(self.scan_info.scan_type).item()),
                "range_mode": range_mode,
                "interpolation": interpolation,
            },
        )
        if extra_attrs:
            section.attrs.update(extra_attrs)
        section["distance"].attrs = {
            "units": "meters",
            "standard_name": "distance_along_section",
            "long_name": "distance along the extracted section",
        }
        section["distance_ground"].attrs = {
            "units": "meters",
            "standard_name": "ground_distance_from_radar",
            "long_name": "ground distance used on the vertical-section x axis",
        }
        section["source_sweep"].attrs = {"long_name": "source sweep index in the original PRD volume"}
        for coord_name in ("x", "y", "z", "lon", "lat", "azimuth", "range", "elevation"):
            section[coord_name].attrs = dict(DEFAULT_METADATA[coord_name])
        if self.fields and field_name in self.fields[0]:
            section[field_name].attrs = dict(self.fields[0][field_name].attrs)
        return section

    def extract_section(
        self,
        start,
        end,
        field_name="dBZ",
        point_units="km",
        interpolation="linear",
        range_mode=None,
        sample_spacing=None,
    ):
        """Extract a Cartesian vertical section for PPI scans using linear interpolation."""
        range_mode = self._resolve_field_range_mode(field_name, range_mode=range_mode)
        if interpolation != "linear":
            raise ValueError("Only linear interpolation is currently supported.")
        scan_type = str(np.asarray(self.scan_info.scan_type).item())
        if scan_type != "ppi":
            raise ValueError("extract_section is only supported for ppi scans.")
        if point_units not in ("km", "m"):
            raise ValueError("point_units must be 'km' or 'm'.")
        scale = 1000.0 if point_units == "km" else 1.0
        start_xy = np.asarray(start, dtype=np.float64) * scale
        end_xy = np.asarray(end, dtype=np.float64) * scale
        spacing = self._resolve_section_sample_spacing(field_name, range_mode, sample_spacing)
        x_line, y_line, distance = self._build_section_line(start_xy, end_xy, spacing)
        radar_altitude = self._site_altitude()
        radar_lon = self._site_longitude()
        radar_lat = self._site_latitude()
        lon_line, lat_line = cartesian_to_geographic_aeqd(x_line, y_line, radar_lon, radar_lat)

        field_rows = []
        x_rows = []
        y_rows = []
        z_rows = []
        lon_rows = []
        lat_rows = []
        az_rows = []
        range_rows = []
        elevation_rows = []
        source_sweeps = []
        for sweep in np.asarray(self._fixed_angle_sort_index, dtype=np.int32):
            sweep_elevation = float(self.scan_info.fixed_angle.values[sweep])
            sweep_azimuth, sweep_ranges_native, sweep_values = self._extract_sweep_volume(
                sweep,
                field_name,
                range_mode=range_mode,
            )
            azimuth, ranges, height = cartesian_to_antenna_cwr(
                x_line,
                y_line,
                sweep_elevation,
                radar_altitude,
                effective_earth_radius=self.effective_earth_radius,
            )
            field_rows.append(
                self._interpolate_ppi_strip_arrays(sweep_azimuth, sweep_ranges_native, sweep_values, azimuth, ranges)
            )
            x_rows.append(x_line)
            y_rows.append(y_line)
            z_rows.append(height)
            lon_rows.append(lon_line)
            lat_rows.append(lat_line)
            az_rows.append(azimuth)
            range_rows.append(ranges)
            elevation_rows.append(np.full(distance.shape, sweep_elevation, dtype=np.float64))
            source_sweeps.append(sweep)
        return self._build_section_dataset(
            field_name=field_name,
            distance=distance,
            field_values=np.stack(field_rows, axis=0),
            x_coords=np.stack(x_rows, axis=0),
            y_coords=np.stack(y_rows, axis=0),
            z_coords=np.stack(z_rows, axis=0),
            lon_coords=np.stack(lon_rows, axis=0),
            lat_coords=np.stack(lat_rows, axis=0),
            azimuth=np.stack(az_rows, axis=0),
            ranges=np.stack(range_rows, axis=0),
            elevation=np.stack(elevation_rows, axis=0),
            source_sweep=np.asarray(source_sweeps, dtype=np.int32),
            section_type="section",
            range_mode=range_mode,
            interpolation=interpolation,
            extra_attrs={
                "start_point": tuple(float(v) for v in start_xy),
                "end_point": tuple(float(v) for v in end_xy),
                "point_units": "meters",
            },
        )

    def extract_section_lonlat(
        self,
        start_lonlat,
        end_lonlat,
        field_name="dBZ",
        interpolation="linear",
        range_mode=None,
        sample_spacing=None,
    ):
        """Extract a vertical section from geographic endpoints."""
        radar_lon = self._site_longitude()
        radar_lat = self._site_latitude()
        start_x, start_y = geographic_to_cartesian_aeqd(start_lonlat[0], start_lonlat[1], radar_lon, radar_lat)
        end_x, end_y = geographic_to_cartesian_aeqd(end_lonlat[0], end_lonlat[1], radar_lon, radar_lat)
        section = self.extract_section(
            (float(np.asarray(start_x).reshape(-1)[0]), float(np.asarray(start_y).reshape(-1)[0])),
            (float(np.asarray(end_x).reshape(-1)[0]), float(np.asarray(end_y).reshape(-1)[0])),
            field_name=field_name,
            point_units="m",
            interpolation=interpolation,
            range_mode=range_mode,
            sample_spacing=sample_spacing,
        )
        section.attrs["start_lonlat"] = tuple(float(v) for v in start_lonlat)
        section.attrs["end_lonlat"] = tuple(float(v) for v in end_lonlat)
        return section

    def _extract_ppi_rhi(self, azimuth, field_name="dBZ", range_mode=None):
        """Extract an RHI-style azimuth cut from a PPI volume."""
        if azimuth is None:
            raise ValueError("azimuth is required when extracting an RHI from a ppi volume.")
        range_mode = self._resolve_field_range_mode(field_name, range_mode=range_mode)
        radar_altitude = self._site_altitude()
        radar_lon = self._site_longitude()
        radar_lat = self._site_latitude()
        target_azimuth = float(azimuth)

        field_rows = []
        x_rows = []
        y_rows = []
        z_rows = []
        lon_rows = []
        lat_rows = []
        az_rows = []
        range_rows = []
        elevation_rows = []
        distance_rows = []
        source_sweeps = []
        for sweep in np.asarray(self._fixed_angle_sort_index, dtype=np.int32):
            sweep_azimuth, sweep_ranges, sweep_values = self._extract_sweep_volume(
                sweep,
                field_name,
                range_mode=range_mode,
            )
            sweep_elevation = float(self.scan_info.fixed_angle.values[sweep])
            azimuth_row = np.full(sweep_ranges.shape, target_azimuth, dtype=np.float64)
            elevation_row = np.full(sweep_ranges.shape, sweep_elevation, dtype=np.float64)
            field_rows.append(
                self._interpolate_ppi_strip_arrays(sweep_azimuth, sweep_ranges, sweep_values, azimuth_row, sweep_ranges)
            )
            x_coords, y_coords, z_coords = antenna_to_cartesian_cwr(
                sweep_ranges,
                target_azimuth,
                sweep_elevation,
                radar_altitude,
                effective_earth_radius=self.effective_earth_radius,
            )
            lon_coords, lat_coords = cartesian_to_geographic_aeqd(x_coords, y_coords, radar_lon, radar_lat)
            x_rows.append(x_coords)
            y_rows.append(y_coords)
            z_rows.append(z_coords)
            lon_rows.append(lon_coords)
            lat_rows.append(lat_coords)
            az_rows.append(azimuth_row)
            range_rows.append(sweep_ranges)
            elevation_rows.append(elevation_row)
            distance_rows.append(np.hypot(x_coords, y_coords))
            source_sweeps.append(sweep)
        distance_axis = np.asarray(max(distance_rows, key=lambda row: row.size), dtype=np.float64)
        section = self._build_section_dataset(
            field_name=field_name,
            distance=distance_axis,
            field_values=self._pad_section_rows(field_rows),
            x_coords=self._pad_section_rows(x_rows),
            y_coords=self._pad_section_rows(y_rows),
            z_coords=self._pad_section_rows(z_rows),
            lon_coords=self._pad_section_rows(lon_rows),
            lat_coords=self._pad_section_rows(lat_rows),
            azimuth=self._pad_section_rows(az_rows),
            ranges=self._pad_section_rows(range_rows),
            elevation=self._pad_section_rows(elevation_rows),
            source_sweep=np.asarray(source_sweeps, dtype=np.int32),
            section_type="rhi",
            range_mode=range_mode,
            interpolation="linear",
            extra_attrs={"target_azimuth": target_azimuth},
        )
        section["distance_ground"] = (("sweep", "distance"), self._pad_section_rows(distance_rows))
        section["distance_ground"].attrs = {
            "units": "meters",
            "standard_name": "ground_distance_from_radar",
            "long_name": "ground distance used on the vertical-section x axis",
        }
        return section

    def _extract_native_rhi(self, field_name="dBZ", range_mode=None):
        """Extract a native RHI section without forcing PPI-only logic."""
        range_mode = self._resolve_field_range_mode(field_name, range_mode=range_mode)
        radar_distance = []
        field_rows = []
        x_rows = []
        y_rows = []
        z_rows = []
        lon_rows = []
        lat_rows = []
        az_rows = []
        range_rows = []
        elevation_rows = []
        source_sweeps = []
        nominal_distance = None
        for sweep in np.asarray(self.scan_info.sweep.values, dtype=np.int32):
            sweep_field = self.get_sweep_field(sweep, field_name, range_mode=range_mode)
            field_values = np.asarray(sweep_field.values, dtype=np.float64)
            x_coords = np.asarray(sweep_field.x.values, dtype=np.float64)
            y_coords = np.asarray(sweep_field.y.values, dtype=np.float64)
            z_coords = np.asarray(sweep_field.z.values, dtype=np.float64)
            lon_coords = np.asarray(sweep_field.lon.values, dtype=np.float64)
            lat_coords = np.asarray(sweep_field.lat.values, dtype=np.float64)
            azimuth = np.broadcast_to(np.asarray(sweep_field.azimuth.values, dtype=np.float64)[:, None], field_values.shape)
            elevation = np.broadcast_to(np.asarray(sweep_field.elevation.values, dtype=np.float64)[:, None], field_values.shape)
            ranges = np.broadcast_to(np.asarray(sweep_field.range.values, dtype=np.float64)[None, :], field_values.shape)
            distance_ground = np.hypot(x_coords, y_coords)
            field_rows.append(field_values)
            x_rows.append(x_coords)
            y_rows.append(y_coords)
            z_rows.append(z_coords)
            lon_rows.append(lon_coords)
            lat_rows.append(lat_coords)
            az_rows.append(azimuth)
            range_rows.append(ranges)
            elevation_rows.append(elevation)
            radar_distance.append(distance_ground)
            source_sweeps.append(np.full(field_values.shape[0], sweep, dtype=np.int32))
            if nominal_distance is None:
                nominal_distance = distance_ground[0]
        section = self._build_section_dataset(
            field_name=field_name,
            distance=np.asarray(nominal_distance, dtype=np.float64),
            field_values=np.concatenate(field_rows, axis=0),
            x_coords=np.concatenate(x_rows, axis=0),
            y_coords=np.concatenate(y_rows, axis=0),
            z_coords=np.concatenate(z_rows, axis=0),
            lon_coords=np.concatenate(lon_rows, axis=0),
            lat_coords=np.concatenate(lat_rows, axis=0),
            azimuth=np.concatenate(az_rows, axis=0),
            ranges=np.concatenate(range_rows, axis=0),
            elevation=np.concatenate(elevation_rows, axis=0),
            source_sweep=np.concatenate(source_sweeps, axis=0),
            section_type="rhi",
            range_mode=range_mode,
            interpolation="native",
        )
        section["distance_ground"] = (("sweep", "distance"), np.concatenate(radar_distance, axis=0))
        section["distance_ground"].attrs = {
            "units": "meters",
            "standard_name": "ground_distance_from_radar",
            "long_name": "ground distance used on the vertical-section x axis",
        }
        return section

    def extract_rhi(self, azimuth=None, field_name="dBZ", range_mode=None):
        """Extract an RHI-style section from PPI or native-RHI scans."""
        range_mode = self._resolve_field_range_mode(field_name, range_mode=range_mode)
        scan_type = str(np.asarray(self.scan_info.scan_type).item())
        if scan_type == "ppi":
            return self._extract_ppi_rhi(azimuth, field_name=field_name, range_mode=range_mode)
        if scan_type == "rhi":
            return self._extract_native_rhi(field_name=field_name, range_mode=range_mode)
        raise ValueError("extract_rhi only supports ppi and rhi scans.")

    def _grid_reflectivity_volume_xy(self, XRange, YRange, level_heights, range_mode=None, blind_method="mask"):
        fillvalue = -999.0
        GridX, GridY = np.meshgrid(XRange, YRange, indexing="ij")
        level_heights = np.asarray(level_heights, dtype=np.float64)
        range_mode = self._resolve_field_range_mode("dBZ", range_mode=range_mode)
        self.get_vol_data(field_name="dBZ", fillvalue=fillvalue, range_mode=range_mode)
        vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, _, _ = self.vol
        GridV = get_CAPPI_3d(
            vol_azimuth,
            vol_range,
            fix_elevation,
            vol_value,
            radar_height,
            GridX.astype(np.float64),
            GridY.astype(np.float64),
            level_heights,
            fillvalue,
            self.effective_earth_radius,
            blind_method,
        )
        return range_mode, fillvalue, level_heights, GridV

    def _grid_reflectivity_volume_lonlat(self, XLon, YLat, level_heights, range_mode=None, blind_method="mask"):
        fillvalue = -999.0
        GridLon, GridLat = np.meshgrid(XLon, YLat, indexing="ij")
        proj = self._get_local_projection()
        GridX, GridY = proj(GridLon, GridLat, inverse=False)
        level_heights = np.asarray(level_heights, dtype=np.float64)
        range_mode = self._resolve_field_range_mode("dBZ", range_mode=range_mode)
        self.get_vol_data(field_name="dBZ", fillvalue=fillvalue, range_mode=range_mode)
        vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, _, _ = self.vol
        GridV = get_CAPPI_3d(
            vol_azimuth,
            vol_range,
            fix_elevation,
            vol_value,
            radar_height,
            GridX.astype(np.float64),
            GridY.astype(np.float64),
            level_heights,
            fillvalue,
            self.effective_earth_radius,
            blind_method,
        )
        return range_mode, fillvalue, level_heights, GridV

    def add_product_CR_xy(self, XRange, YRange, range_mode=None):
        """
        Compute composite reflectivity on a Cartesian grid.
        :param XRange: np.ndarray, 1d, units:meters
        :param YRange: np.ndarray, 1d, units:meters
        :return:
        """
        GridX, GridY = np.meshgrid(XRange, YRange, indexing="ij")
        range_mode = self._resolve_field_range_mode("dBZ", range_mode=range_mode)
        self.get_vol_data(field_name="dBZ", fillvalue=-999., range_mode=range_mode)
        vol_azimuth, vol_range, fix_elevation, vol_value, radar_height,\
        radar_lon_0, radar_lat_0 = self.vol
        fillvalue = -999.
        GridV = get_CR_xy(vol_azimuth, vol_range, fix_elevation, vol_value,\
                          radar_height, GridX.astype(np.float64), GridY.astype(np.float64), -999.,
                          self.effective_earth_radius)
        product_name = self._product_name("CR", range_mode)
        self.product.coords["x_cr"] = XRange
        self.product.coords["y_cr"] = YRange
        self.product[product_name] = (('x_cr', 'y_cr'), np.where(GridV==fillvalue, np.nan, GridV))
        self.product.coords["x_cr"].attrs = {'units': 'meters',
                                            'standard_name': 'CR_product_x_axis ',
                                            'long_name': 'east_distance_from_radar',
                                            'axis': 'xy_coordinate',
                                            'comment': 'Distance from radar in east'}
        self.product.coords["y_cr"].attrs = {'units': 'meters',
                                             'standard_name': 'CR_product_y_axis ',
                                             'long_name': 'north_distance_from_radar',
                                             'axis': 'xy_coordinate',
                                             'comment': 'Distance from radar in north'}
        self.product[product_name].attrs = {'units': 'dBZ',
                                    'standard_name': 'Composite_reflectivity_factor',
                                    'long_name': 'Composite_reflectivity_factor',
                                    'axis': 'xy_coordinate',
                                    'comment': 'Maximum reflectance of all level',}

    def add_product_CAPPI_xy(self, XRange, YRange, level_height, range_mode=None):
        """
        Compute a CAPPI product on a Cartesian grid.
        :param XRange: np.ndarray, 1d, units:meters
        :param YRange: np.ndarray, 1d, units:meters
        :param level_height: target height, units: meters
        :return:
        """
        GridX, GridY = np.meshgrid(XRange, YRange, indexing="ij")
        range_mode = self._resolve_field_range_mode("dBZ", range_mode=range_mode)
        self.get_vol_data(field_name="dBZ", fillvalue=-999., range_mode=range_mode)
        vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, \
        radar_lon_0, radar_lat_0 = self.vol
        fillvalue = -999.
        GridV = get_CAPPI_xy(vol_azimuth, vol_range, fix_elevation, vol_value, radar_height,
                             GridX.astype(np.float64), GridY.astype(np.float64), level_height, fillvalue,
                             self.effective_earth_radius)
        product_name = self._product_name("CAPPI_%d" % level_height, range_mode)
        self.product.coords["x_cappi_%d"%level_height] = XRange
        self.product.coords["y_cappi_%d"%level_height] = YRange
        self.product[product_name] = (("x_cappi_%d"%level_height, "y_cappi_%d"%level_height),
                                      np.where(GridV==fillvalue, np.nan, GridV))
        self.product.coords["x_cappi_%d"%level_height].attrs = {'units': 'meters',
                                             'standard_name': 'CAPPI_product_x_axis ',
                                             'long_name': 'east_distance_from_radar',
                                             'axis': 'xy_coordinate',
                                             'comment': 'Distance from radar in east'}
        self.product.coords["y_cappi_%d"%level_height].attrs = {'units': 'meters',
                                             'standard_name': 'CAPPI_product_y_axis ',
                                             'long_name': 'north_distance_from_radar',
                                             'axis': 'xy_coordinate',
                                             'comment': 'Distance from radar in north'}
        self.product[product_name].attrs = {'units': 'dBZ',
                                'standard_name': 'Constant_altitude_plan_position_indicator',
                                'long_name': 'Constant_altitude_plan_position_indicator',
                                'axis': 'xy_coordinate',
                                'comment': 'CAPPI of level %d m.'%level_height, }

    def add_product_CAPPI_3d_xy(self, XRange, YRange, level_heights, range_mode=None):
        """
        Compute a 3D CAPPI volume on a Cartesian grid.
        :param XRange: np.ndarray, 1d, units:meters
        :param YRange: np.ndarray, 1d, units:meters
        :param level_heights: np.ndarray, 1d, units:meters
        :return:
        """
        GridX, GridY = np.meshgrid(XRange, YRange, indexing="ij")
        range_mode = self._resolve_field_range_mode("dBZ", range_mode=range_mode)
        self.get_vol_data(field_name="dBZ", fillvalue=-999., range_mode=range_mode)
        vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, \
        radar_lon_0, radar_lat_0 = self.vol
        fillvalue = -999.
        level_heights = np.asarray(level_heights, dtype=np.float64)
        GridV = get_CAPPI_3d(
            vol_azimuth,
            vol_range,
            fix_elevation,
            vol_value,
            radar_height,
            GridX.astype(np.float64),
            GridY.astype(np.float64),
            level_heights,
            fillvalue,
            self.effective_earth_radius,
        )
        self.product.coords["z_cappi"] = level_heights
        self.product.coords["x_cappi_3d"] = XRange
        self.product.coords["y_cappi_3d"] = YRange
        product_name = self._product_name("CAPPI_3D", range_mode)
        self.product[product_name] = (
            ("z_cappi", "x_cappi_3d", "y_cappi_3d"),
            np.where(GridV == fillvalue, np.nan, GridV),
        )
        self.product.coords["z_cappi"].attrs = {
            'units': 'meters',
            'standard_name': 'CAPPI_product_z_axis ',
            'long_name': 'height_above_mean_sea_level',
            'axis': 'z_coordinate',
        }
        self.product.coords["x_cappi_3d"].attrs = {
            'units': 'meters',
            'standard_name': 'CAPPI_product_x_axis ',
            'long_name': 'east_distance_from_radar',
            'axis': 'xy_coordinate',
            'comment': 'Distance from radar in east',
        }
        self.product.coords["y_cappi_3d"].attrs = {
            'units': 'meters',
            'standard_name': 'CAPPI_product_y_axis ',
            'long_name': 'north_distance_from_radar',
            'axis': 'xy_coordinate',
            'comment': 'Distance from radar in north',
        }
        self.product[product_name].attrs = {
            'units': 'dBZ',
            'standard_name': 'Constant_altitude_plan_position_indicator',
            'long_name': '3D_constant_altitude_plan_position_indicator',
            'axis': 'xyz_coordinate',
            'comment': '3D CAPPI volume for the current radar',
        }

    def add_product_VIL_xy(
        self,
        XRange,
        YRange,
        level_heights,
        range_mode=None,
        blind_method="mask",
        min_dbz=18.0,
        max_dbz_cap=56.0,
    ):
        """Compute single-radar VIL on a Cartesian grid."""
        range_mode, fillvalue, level_heights, GridV = self._grid_reflectivity_volume_xy(
            XRange,
            YRange,
            level_heights,
            range_mode=range_mode,
            blind_method=blind_method,
        )
        product_name = self._product_name("VIL", range_mode)
        self.product.coords["x_vil"] = XRange
        self.product.coords["y_vil"] = YRange
        self.product[product_name] = (
            ("x_vil", "y_vil"),
            derive_vil(GridV, level_heights, fillvalue=fillvalue, min_dbz=min_dbz, max_dbz_cap=max_dbz_cap),
        )
        self.product.coords["x_vil"].attrs = {
            'units': 'meters',
            'standard_name': 'VIL_product_x_axis ',
            'long_name': 'east_distance_from_radar',
            'axis': 'xy_coordinate',
            'comment': 'Distance from radar in east',
        }
        self.product.coords["y_vil"].attrs = {
            'units': 'meters',
            'standard_name': 'VIL_product_y_axis ',
            'long_name': 'north_distance_from_radar',
            'axis': 'xy_coordinate',
            'comment': 'Distance from radar in north',
        }
        self.product[product_name].attrs = {
            'units': 'kg m-2',
            'standard_name': 'Vertically_integrated_liquid_water',
            'long_name': 'Vertically_integrated_liquid_water',
            'axis': 'xy_coordinate',
            'source_levels': json.dumps(level_heights.tolist()),
            'min_dbz': float(min_dbz),
            'max_dbz_cap': float(max_dbz_cap),
            'references': json.dumps(PRODUCT_REFERENCE_NOTES["VIL"], ensure_ascii=False),
        }

    def add_product_ET_xy(
        self,
        XRange,
        YRange,
        level_heights,
        range_mode=None,
        blind_method="mask",
        threshold_dbz=18.0,
    ):
        """Compute single-radar echo-top height on a Cartesian grid."""
        range_mode, fillvalue, level_heights, GridV = self._grid_reflectivity_volume_xy(
            XRange,
            YRange,
            level_heights,
            range_mode=range_mode,
            blind_method=blind_method,
        )
        et_value, et_topped = derive_et(
            GridV,
            level_heights,
            fillvalue=fillvalue,
            threshold_dbz=threshold_dbz,
            return_topped=True,
        )
        product_name = self._product_name("ET", range_mode)
        topped_name = self._product_name("ET_TOPPED", range_mode)
        self.product.coords["x_et"] = XRange
        self.product.coords["y_et"] = YRange
        self.product[product_name] = (("x_et", "y_et"), et_value)
        self.product[topped_name] = (("x_et", "y_et"), np.asarray(et_topped, dtype=np.uint8))
        self.product.coords["x_et"].attrs = {
            'units': 'meters',
            'standard_name': 'ET_product_x_axis ',
            'long_name': 'east_distance_from_radar',
            'axis': 'xy_coordinate',
            'comment': 'Distance from radar in east',
        }
        self.product.coords["y_et"].attrs = {
            'units': 'meters',
            'standard_name': 'ET_product_y_axis ',
            'long_name': 'north_distance_from_radar',
            'axis': 'xy_coordinate',
            'comment': 'Distance from radar in north',
        }
        self.product[product_name].attrs = {
            'units': 'meters',
            'standard_name': 'Echo_top_height',
            'long_name': 'Echo_top_height',
            'axis': 'xy_coordinate',
            'source_levels': json.dumps(level_heights.tolist()),
            'threshold_dbz': float(threshold_dbz),
            'references': json.dumps(PRODUCT_REFERENCE_NOTES["ET"], ensure_ascii=False),
        }
        self.product[topped_name].attrs = {
            'units': '1',
            'standard_name': 'Echo_top_topped_flag',
            'long_name': 'Echo_top_topped_flag',
            'axis': 'xy_coordinate',
            'flag_values': np.array([0, 1], dtype=np.uint8),
            'flag_meanings': 'resolved topped_at_highest_sampled_level',
            'source_levels': json.dumps(level_heights.tolist()),
            'threshold_dbz': float(threshold_dbz),
            'references': json.dumps(PRODUCT_REFERENCE_NOTES["ET"], ensure_ascii=False),
        }

    def add_product_VWP(self, sweeps=None, field_name=None, range_mode="aligned", **kwargs):
        """
        Compute and store a vertical wind profile product derived from VAD.
        """
        profile = self.retrieve_vwp(
            sweeps=sweeps,
            field_name=field_name,
            range_mode=range_mode,
            **kwargs
        )
        self.product.coords["z_vwp"] = np.asarray(profile["height"].values, dtype=np.float64)
        self.product["VWP_u"] = (("z_vwp",), np.asarray(profile["u"].values, dtype=np.float64))
        self.product["VWP_v"] = (("z_vwp",), np.asarray(profile["v"].values, dtype=np.float64))
        self.product["VWP_speed"] = (("z_vwp",), np.asarray(profile["wind_speed"].values, dtype=np.float64))
        self.product["VWP_direction"] = (("z_vwp",), np.asarray(profile["wind_direction"].values, dtype=np.float64))
        self.product["VWP_sample_count"] = (("z_vwp",), np.asarray(profile["sample_count"].values, dtype=np.int32))
        self.product["VWP_fit_rmse"] = (("z_vwp",), np.asarray(profile["fit_rmse"].values, dtype=np.float64))
        self.product.coords["z_vwp"].attrs = dict(profile["height"].attrs)
        self.product["VWP_u"].attrs = dict(profile["u"].attrs)
        self.product["VWP_v"].attrs = dict(profile["v"].attrs)
        self.product["VWP_speed"].attrs = dict(profile["wind_speed"].attrs)
        self.product["VWP_direction"].attrs = dict(profile["wind_direction"].attrs)
        self.product["VWP_sample_count"].attrs = dict(profile["sample_count"].attrs)
        self.product["VWP_fit_rmse"].attrs = dict(profile["fit_rmse"].attrs)
        for variable in ("VWP_u", "VWP_v", "VWP_speed", "VWP_direction", "VWP_sample_count", "VWP_fit_rmse"):
            self.product[variable].attrs["source_profile"] = "VWP"
            self.product[variable].attrs["references"] = profile.attrs.get("references", "[]")
        self.product.attrs["VWP_config"] = json.dumps(
            {
                "method": "VWP",
                "range_mode": range_mode,
                "velocity_field_used_by_sweep": profile.attrs.get("velocity_field_used_by_sweep", "{}"),
            },
            ensure_ascii=False,
            sort_keys=True,
        )
        return profile

    def add_product_CR_lonlat(self, XLon, YLat, range_mode=None):
        """
        Compute composite reflectivity on a longitude/latitude grid.
        :param XLon:np.ndarray, 1d, units:degree
        :param YLat:np.ndarray, 1d, units:degree
        :return:
        """
        fillvalue = -999.
        GridLon, GridLat = np.meshgrid(XLon, YLat, indexing="ij")
        proj = self._get_local_projection()
        GridX, GridY = proj(GridLon, GridLat, inverse=False)
        range_mode = self._resolve_field_range_mode("dBZ", range_mode=range_mode)
        self.get_vol_data(field_name="dBZ", fillvalue=fillvalue, range_mode=range_mode)
        vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, \
        radar_lon_0, radar_lat_0 = self.vol
        GridV = get_CR_xy(vol_azimuth, vol_range, fix_elevation, vol_value, \
                          radar_height, GridX.astype(np.float64), GridY.astype(np.float64), -999.,
                          self.effective_earth_radius)
        product_name = self._product_name("CR_geo", range_mode)
        self.product.coords["lon_cr"] = XLon
        self.product.coords["lat_cr"] = YLat
        self.product[product_name] = (('lon_cr', 'lat_cr'), np.where(GridV == fillvalue, np.nan, GridV))
        self.product.coords["lon_cr"].attrs = {'units': 'degrees',
                                             'standard_name': 'CR_product_lon_axis ',
                                             'long_name': 'longitude_cr',
                                             'axis': 'lonlat_coordinate'}
        self.product.coords["lat_cr"].attrs = {'units': 'degrees',
                                             'standard_name': 'CR_product_lat_axis ',
                                             'long_name': 'latitude_cr',
                                             'axis': 'lonlat_coordinate'}
        self.product[product_name].attrs = {'units': 'dBZ',
                                    'standard_name': 'Composite_reflectivity_factor_lonlat_grid',
                                    'long_name': 'Composite_reflectivity_factor',
                                    'axis': 'lonlat_coordinate',
                                    'comment': 'Maximum reflectance of all level', }

    def add_product_CAPPI_lonlat(self, XLon, YLat, level_height, range_mode=None):
        """
        Compute CAPPI on a longitude/latitude grid.
        :param XLon:np.ndarray, 1d, units:degrees
        :param YLat:np.ndarray, 1d, units:degrees
        :param level_height: target height, units: meters
        :return:
        """
        fillvalue = -999.
        GridLon, GridLat = np.meshgrid(XLon, YLat, indexing="ij")
        proj = self._get_local_projection()
        GridX, GridY = proj(GridLon, GridLat, inverse=False)
        range_mode = self._resolve_field_range_mode("dBZ", range_mode=range_mode)
        self.get_vol_data(field_name="dBZ", fillvalue=fillvalue, range_mode=range_mode)
        vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, \
        radar_lon_0, radar_lat_0 = self.vol
        GridV = get_CAPPI_xy(vol_azimuth, vol_range, fix_elevation, vol_value, radar_height,
                             GridX.astype(np.float64), GridY.astype(np.float64), level_height, fillvalue,
                             self.effective_earth_radius)
        product_name = self._product_name("CAPPI_geo_%d" % level_height, range_mode)
        self.product.coords["lon_cappi_%d" % level_height] = XLon
        self.product.coords["lat_cappi_%d" % level_height] = YLat
        self.product[product_name] = (("lon_cappi_%d" % level_height, "lat_cappi_%d" % level_height),
                                      np.where(GridV == fillvalue, np.nan, GridV))
        self.product.coords["lon_cappi_%d" % level_height].attrs = {'units': 'degrees',
                                                                  'standard_name': 'CAPPI_product_lon_axis ',
                                                                  'long_name': 'longitude_CAPPI',
                                                                  'axis': 'lonlat_coordinate',}
        self.product.coords["lat_cappi_%d" % level_height].attrs = {'units': 'degrees',
                                                                  'standard_name': 'CAPPI_product_lat_axis ',
                                                                  'long_name': 'latitude_CAPPI',
                                                                  'axis': 'lonlat_coordinate',}
        self.product[product_name].attrs = {'units': 'dBZ',
                                            'standard_name': 'Constant_altitude_plan_position_indicator',
                                            'long_name': 'Constant_altitude_plan_position_indicator',
                                            'axis': 'lonlat_coordinate',
                                            'comment': 'CAPPI of level %d m' % level_height, }

    def add_product_VIL_lonlat(
        self,
        XLon,
        YLat,
        level_heights,
        range_mode=None,
        blind_method="mask",
        min_dbz=18.0,
        max_dbz_cap=56.0,
    ):
        """Compute single-radar VIL on a longitude/latitude grid."""
        range_mode, fillvalue, level_heights, GridV = self._grid_reflectivity_volume_lonlat(
            XLon,
            YLat,
            level_heights,
            range_mode=range_mode,
            blind_method=blind_method,
        )
        product_name = self._product_name("VIL_geo", range_mode)
        self.product.coords["lon_vil"] = XLon
        self.product.coords["lat_vil"] = YLat
        self.product[product_name] = (
            ("lon_vil", "lat_vil"),
            derive_vil(GridV, level_heights, fillvalue=fillvalue, min_dbz=min_dbz, max_dbz_cap=max_dbz_cap),
        )
        self.product.coords["lon_vil"].attrs = {
            'units': 'degrees',
            'standard_name': 'VIL_product_lon_axis ',
            'long_name': 'longitude_vil',
            'axis': 'lonlat_coordinate',
        }
        self.product.coords["lat_vil"].attrs = {
            'units': 'degrees',
            'standard_name': 'VIL_product_lat_axis ',
            'long_name': 'latitude_vil',
            'axis': 'lonlat_coordinate',
        }
        self.product[product_name].attrs = {
            'units': 'kg m-2',
            'standard_name': 'Vertically_integrated_liquid_water_lonlat_grid',
            'long_name': 'Vertically_integrated_liquid_water',
            'axis': 'lonlat_coordinate',
            'source_levels': json.dumps(level_heights.tolist()),
            'min_dbz': float(min_dbz),
            'max_dbz_cap': float(max_dbz_cap),
            'references': json.dumps(PRODUCT_REFERENCE_NOTES["VIL"], ensure_ascii=False),
        }

    def add_product_ET_lonlat(
        self,
        XLon,
        YLat,
        level_heights,
        range_mode=None,
        blind_method="mask",
        threshold_dbz=18.0,
    ):
        """Compute single-radar echo-top height on a longitude/latitude grid."""
        range_mode, fillvalue, level_heights, GridV = self._grid_reflectivity_volume_lonlat(
            XLon,
            YLat,
            level_heights,
            range_mode=range_mode,
            blind_method=blind_method,
        )
        et_value, et_topped = derive_et(
            GridV,
            level_heights,
            fillvalue=fillvalue,
            threshold_dbz=threshold_dbz,
            return_topped=True,
        )
        product_name = self._product_name("ET_geo", range_mode)
        topped_name = self._product_name("ET_TOPPED_geo", range_mode)
        self.product.coords["lon_et"] = XLon
        self.product.coords["lat_et"] = YLat
        self.product[product_name] = (("lon_et", "lat_et"), et_value)
        self.product[topped_name] = (("lon_et", "lat_et"), np.asarray(et_topped, dtype=np.uint8))
        self.product.coords["lon_et"].attrs = {
            'units': 'degrees',
            'standard_name': 'ET_product_lon_axis ',
            'long_name': 'longitude_et',
            'axis': 'lonlat_coordinate',
        }
        self.product.coords["lat_et"].attrs = {
            'units': 'degrees',
            'standard_name': 'ET_product_lat_axis ',
            'long_name': 'latitude_et',
            'axis': 'lonlat_coordinate',
        }
        self.product[product_name].attrs = {
            'units': 'meters',
            'standard_name': 'Echo_top_height_lonlat_grid',
            'long_name': 'Echo_top_height',
            'axis': 'lonlat_coordinate',
            'source_levels': json.dumps(level_heights.tolist()),
            'threshold_dbz': float(threshold_dbz),
            'references': json.dumps(PRODUCT_REFERENCE_NOTES["ET"], ensure_ascii=False),
        }
        self.product[topped_name].attrs = {
            'units': '1',
            'standard_name': 'Echo_top_topped_flag_lonlat_grid',
            'long_name': 'Echo_top_topped_flag',
            'axis': 'lonlat_coordinate',
            'flag_values': np.array([0, 1], dtype=np.uint8),
            'flag_meanings': 'resolved topped_at_highest_sampled_level',
            'source_levels': json.dumps(level_heights.tolist()),
            'threshold_dbz': float(threshold_dbz),
            'references': json.dumps(PRODUCT_REFERENCE_NOTES["ET"], ensure_ascii=False),
        }

    def get_vol_data(self, field_name="dBZ", fillvalue=-999., range_mode=None, max_range_km=None):
        """
        Prepare azimuth-sorted volume data used by Cartesian products.
        :return:
        """
        range_mode = self._resolve_field_range_mode(field_name, range_mode=range_mode)
        max_range_km = None if max_range_km is None else float(max_range_km)
        cache_key = (field_name, float(fillvalue), range_mode, max_range_km)
        if cache_key in self._vol_cache:
            self.vol = self._vol_cache[cache_key]
            return
        sweep_order = self._fixed_angle_sort_index
        fixed_elevation = self.scan_info["fixed_angle"].values[sweep_order]
        vol_azimuth = []
        vol_range = []
        vol_value = []
        for sweep in sweep_order:
            azimuth, ranges, values = self._extract_sweep_volume(
                sweep,
                field_name,
                range_mode=range_mode,
                max_range_km=max_range_km,
            )
            vol_azimuth.append(azimuth)
            vol_range.append(ranges)
            vol_value.append(np.where(np.isnan(values), fillvalue, values).astype(np.float64))
        radar_height = float(self.scan_info["altitude"].values)
        radar_lon_0 = float(self.scan_info["longitude"].values)
        radar_lat_0 = float(self.scan_info["latitude"].values)
        self.vol = vol_azimuth, vol_range, fixed_elevation.astype(np.float64), vol_value, radar_height, radar_lon_0, radar_lat_0
        self._vol_cache[cache_key] = self.vol

    def get_RHI_data(self, az, field_name="dBZ", range_mode=None):
        """
        Extract an RHI-style vertical cross-section at the nearest azimuth.
        :param az:
        :param field_name:
        :return:
        """
        section = self.extract_rhi(azimuth=az, field_name=field_name, range_mode=range_mode)
        mesh_rhi = []
        mesh_range = []
        mesh_z = []
        if section.attrs.get("scan_type") == "rhi":
            for idx in range(section.sizes["sweep"]):
                mesh_rhi.append(np.asarray(section[field_name].isel(sweep=idx).values, dtype=np.float64).reshape(1, -1))
                mesh_range.append(np.asarray(section["distance_ground"].isel(sweep=idx).values, dtype=np.float64).reshape(1, -1))
                mesh_z.append(np.asarray(section["z"].isel(sweep=idx).values, dtype=np.float64).reshape(1, -1))
            return mesh_range, mesh_z, mesh_rhi

        radar_altitude = float(np.asarray(self.scan_info.altitude.values))
        for idx in range(section.sizes["sweep"]):
            source_sweep = int(section["source_sweep"].values[idx])
            ranges = np.asarray(section["range"].isel(sweep=idx).values, dtype=np.float64)
            azimuth = np.asarray(section["azimuth"].isel(sweep=idx).values[0], dtype=np.float64)
            elevation = np.asarray(section["elevation"].isel(sweep=idx).values[0], dtype=np.float64)
            beam_width = float(self.scan_info.beam_width.values[source_sweep])
            x_coords, y_coords, z_coords = antenna_vectors_to_cartesian_rhi(
                ranges,
                np.asarray([azimuth], dtype=np.float64),
                np.asarray([elevation], dtype=np.float64),
                radar_altitude,
                BeamWidth=beam_width,
                effective_earth_radius=self.effective_earth_radius,
            )
            mesh_range.append(np.sqrt(x_coords ** 2 + y_coords ** 2))
            mesh_z.append(z_coords)
            mesh_rhi.append(np.asarray(section[field_name].isel(sweep=idx).values, dtype=np.float64).reshape(1, -1))
        return mesh_range, mesh_z, mesh_rhi

    def get_vcs_data(self, start_point, end_point, field_name, range_mode=None):
        """
        Extract a vertical cross-section along a Cartesian line segment.
        :param start_point:
        :param end_point:
        :param field_name:
        :return:
        """
        section = self.extract_section(
            start_point,
            end_point,
            field_name=field_name,
            point_units="m",
            interpolation="linear",
            range_mode=range_mode,
        )
        mesh_xy = []
        mesh_z = []
        mesh_vcs = []
        radar_altitude = float(np.asarray(self.scan_info.altitude.values))
        for idx in range(section.sizes["sweep"]):
            source_sweep = int(section["source_sweep"].values[idx])
            beam_width = float(self.scan_info.beam_width.values[source_sweep])
            distance = np.asarray(section["distance"].values, dtype=np.float64)
            azimuth = np.asarray(section["azimuth"].isel(sweep=idx).values, dtype=np.float64)
            ranges = np.asarray(section["range"].isel(sweep=idx).values, dtype=np.float64)
            elevation = np.asarray(section["elevation"].isel(sweep=idx).values, dtype=np.float64)
            _, _, z_coords = antenna_vectors_to_cartesian_vcs(
                ranges,
                azimuth,
                elevation,
                radar_altitude,
                BeamWidth=beam_width,
                effective_earth_radius=self.effective_earth_radius,
            )
            mesh_xy.append(np.stack([distance, distance], axis=0))
            mesh_z.append(z_coords)
            mesh_vcs.append(np.asarray(section[field_name].isel(sweep=idx).values, dtype=np.float64).reshape(1, -1))
        return mesh_xy, mesh_z, mesh_vcs


def grid_3d_network_xy(
    radars,
    XRange,
    YRange,
    level_heights,
    field_name="dBZ",
    fillvalue=-999.,
    lon_0=None,
    lat_0=None,
    range_mode=None,
):
    """
    Compute a 3D multi-radar CAPPI network on a shared Cartesian grid.
    Radar values are combined using maximum-value compositing.
    """
    if not radars:
        raise ValueError("radars must contain at least one PRD instance.")

    XRange = np.asarray(XRange, dtype=np.float64)
    YRange = np.asarray(YRange, dtype=np.float64)
    level_heights = np.asarray(level_heights, dtype=np.float64)

    radar_lons = np.array([float(radar.scan_info["longitude"].values) for radar in radars], dtype=np.float64)
    radar_lats = np.array([float(radar.scan_info["latitude"].values) for radar in radars], dtype=np.float64)
    if lon_0 is None:
        lon_0 = float(np.mean(radar_lons))
    if lat_0 is None:
        lat_0 = float(np.mean(radar_lats))

    proj = pyproj.Proj({"proj": "aeqd", "lon_0": lon_0, "lat_0": lat_0})
    GridX, GridY = np.meshgrid(XRange, YRange, indexing="ij")

    vol_azimuth = []
    vol_range = []
    fix_elevation = []
    vol_value = []
    radar_x = np.empty(len(radars), dtype=np.float64)
    radar_y = np.empty(len(radars), dtype=np.float64)
    radar_height = np.empty(len(radars), dtype=np.float64)
    radar_effective_earth_radius = []

    for iradar, radar in enumerate(radars):
        effective_range_mode = radar._resolve_field_range_mode(field_name, range_mode=range_mode)
        radar.get_vol_data(field_name=field_name, fillvalue=fillvalue, range_mode=effective_range_mode)
        radar_vol_azimuth, radar_vol_range, radar_fix_elevation, radar_vol_value, \
        radar_height_value, radar_lon, radar_lat = radar.vol
        radar_x[iradar], radar_y[iradar] = proj(radar_lon, radar_lat, inverse=False)
        radar_height[iradar] = radar_height_value
        radar_effective_earth_radius.append(getattr(radar, "effective_earth_radius", None))
        vol_azimuth.append(radar_vol_azimuth)
        vol_range.append(radar_vol_range)
        fix_elevation.append(radar_fix_elevation)
        vol_value.append(radar_vol_value)

    GridV = get_mosaic_CAPPI_3d(
        vol_azimuth,
        vol_range,
        fix_elevation,
        vol_value,
        radar_x,
        radar_y,
        radar_height,
        GridX.astype(np.float64),
        GridY.astype(np.float64),
        level_heights,
        fillvalue,
        radar_effective_earth_radius,
    )

    product = xr.Dataset(coords={"z": level_heights, "x": XRange, "y": YRange})
    product["network_3d"] = (("z", "x", "y"), np.where(GridV == fillvalue, np.nan, GridV))
    product["network_3d"].attrs = {
        "units": "dBZ",
        "long_name": "3D_multi_radar_CAPPI_network",
        "comment": "Maximum-value multi-radar CAPPI mosaic on a shared Cartesian grid",
        "grid_origin_lon": lon_0,
        "grid_origin_lat": lat_0,
        "field_name": field_name,
        "range_mode": radars[0]._resolve_field_range_mode(field_name, range_mode=range_mode),
    }
    return product

class AzimuthSortedPRD:
    """Lightweight azimuth-sorted view of a PRD volume."""

    __slots__ = ("scan_info", "fields", "product", "effective_earth_radius")

    def __init__(self, scan_info=None, fields=None, product=None, effective_earth_radius=None):
        self.scan_info = scan_info
        self.fields = [] if fields is None else list(fields)
        self.product = xr.Dataset() if product is None else product
        self.effective_earth_radius = effective_earth_radius


PRD_AZ = AzimuthSortedPRD
