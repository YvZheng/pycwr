# -*- coding: utf-8 -*-
"""
为了适应中国雷达在不同仰角的探测距离不同以及前几层仰角 dop和ref分开扫描的问题
提出NuistRadar Object，以方便后续的算法及绘图
"""
import numpy as np
import xarray as xr
from ..configure.default_config import DEFAULT_METADATA, FILL_VALUE, CINRAD_field_mapping
from ..core.transforms import antenna_vectors_to_cartesian
from scipy import spatial
from ..interp.RadarInterp import radar_interp2d_var, radar_interp2d

class NuistRadar(object):
    """
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
    nrays : int
        Number of rays in the volume.
    nsweeps : int
        Number of sweep in the volume.

    """

    def __init__(self, fields,  scan_type, time, range, azimuth, elevation,latitude,
                 longitude, altitude, sweep_start_ray_index, sweep_end_ray_index,
                 fixed_angle, bins_per_sweep, nyquist_velocity, frequency, unambiguous_range,
                 nrays, nsweeps, sitename):

        super(NuistRadar, self).__init__()
        keys = fields.keys()
        self.fields = []
        for idx, (istart, iend) in enumerate(zip(sweep_start_ray_index, sweep_end_ray_index)):
            isweep_data = xr.Dataset(coords={'azimuth': (['time', ], azimuth[istart:iend+1]),
                                            'elevation': (['time',], elevation[istart:iend+1]),
                                            'range': range[:bins_per_sweep[idx]], 'time': time[istart:iend+1]})
            isweep_data.azimuth.attrs = DEFAULT_METADATA['azimuth']
            isweep_data.elevation.attrs = DEFAULT_METADATA['elevation']
            isweep_data.range.attrs = DEFAULT_METADATA['range']
            isweep_data.time.attrs = DEFAULT_METADATA['time']
            for ikey in keys:
                isweep_data[ikey] = (['time','range'], fields[ikey][istart:iend+1, :bins_per_sweep[idx]])
                isweep_data[ikey].attrs = DEFAULT_METADATA[CINRAD_field_mapping[ikey]]
            self.fields.append(isweep_data)
        self.scan_info = xr.Dataset(data_vars={"latitude":latitude,"longitude":longitude,
                        "altitude":altitude,"scan_type":scan_type,  "frequency":frequency,
                         "nyquist_velocity":(['sweep',], nyquist_velocity),
                        "unambiguous_range":(['sweep',], unambiguous_range),
                        "rays_per_sweep": (['sweep',], sweep_end_ray_index-sweep_start_ray_index+1),
                        "fixed_angle": (["sweep",], fixed_angle)},
                        coords={"sweep": np.arange(nsweeps, dtype=int)})

        self.scan_info['latitude'].attrs = DEFAULT_METADATA['latitude']
        self.scan_info['longitude'].attrs = DEFAULT_METADATA['longitude']
        self.scan_info['altitude'].attrs = DEFAULT_METADATA['altitude']
        self.scan_info['scan_type'].attrs = DEFAULT_METADATA['scan_type']
        self.scan_info['frequency'].attrs = DEFAULT_METADATA['frequency']
        self.scan_info['nyquist_velocity'].attrs = DEFAULT_METADATA['nyquist_velocity']
        self.scan_info['unambiguous_range'].attrs = DEFAULT_METADATA['unambiguous_range']
        self.scan_info['rays_per_sweep'].attrs = DEFAULT_METADATA['rays_per_sweep']
        self.scan_info['fixed_angle'].attrs = DEFAULT_METADATA['fixed_angle']
        self.nsweeps = nsweeps
        self.nrays = nrays
        self.sitename = sitename

    def get_vertical_section(self, start_point, end_point, field_name):
        """
        :param start_point: units:m
        :param end_point:  units:m
        :return:
        """
        start_x, start_y = start_point
        end_x, end_y = end_point
        bins_res = (self.fields[0].range[1] - self.fields[0].range[0]).values
        start_end_dis = np.sqrt((start_x-end_x)**2 + (start_y-end_y)**2)
        npoints = int(start_end_dis/(bins_res/2.) + 1)
        x_line = np.linspace(start_x, end_x, npoints)
        y_line = np.linspace(start_y, end_y, npoints)
        target = np.c_[x_line, y_line]
        xy_line_1d = np.linspace(0, start_end_dis, npoints)
        z_values = []
        field_values = []
        #先对剖线取最邻近点
        for ifield in self.fields:
            _x, _y, _z = antenna_vectors_to_cartesian(ifield.range.values, \
                                ifield.azimuth.values, ifield.elevation.values)
            kdtree = spatial.cKDTree(np.c_[_x.ravel(), _y.ravel()])
            _distance, _idx = kdtree.query(target, k=1, n_jobs=-1)
            _z_value = np.where(_distance > bins_res*2, np.nan, _z.ravel()[_idx])
            _field_value = np.where(_distance > bins_res*2, np.nan, ifield[field_name].values.ravel()[_idx])
            z_values.append(_z_value)
            field_values.append(_field_value)
        z_values = np.asarray(z_values)
        field_values = np.asarray(field_values)
        xy_line_values = np.stack([xy_line_1d,]*len(self.fields), axis=0)
        mask_flag = (np.isnan(z_values.ravel()) | np.isnan(field_values.ravel()))
        z_values_ravel = z_values.ravel()[~mask_flag]
        xy_line_values_ravel = xy_line_values.ravel()[~mask_flag]
        field_values_ravel = field_values.ravel()[~mask_flag]
        mesh_xy, mesh_z = np.mgrid[0:start_end_dis:npoints*1j, 0:np.nanmax(z_values):bins_res/2.] #生成剖面的网格
        grid_field = radar_interp2d_var(np.c_[xy_line_values_ravel, z_values_ravel], field_values_ravel,\
                              (mesh_xy, mesh_z), bandwidth=360/self.fields[0].azimuth.size)
        return mesh_xy, mesh_z, grid_field
