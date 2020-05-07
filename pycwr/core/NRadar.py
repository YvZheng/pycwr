# -*- coding: utf-8 -*-
"""
为了适应中国雷达在不同仰角的探测距离不同以及前几层仰角 dop和ref分开扫描的问题
提出PRD Object，以方便后续的算法及绘图
"""
import numpy as np
import xarray as xr
import pyproj
from ..configure.default_config import DEFAULT_METADATA, CINRAD_field_mapping
from ..core.transforms import  cartesian_to_geographic_aeqd,\
    antenna_vectors_to_cartesian_cwr, antenna_vectors_to_cartesian_rhi, cartesian_to_antenna_cwr,\
    antenna_vectors_to_cartesian_vcs
from .RadarGridC import get_CR_xy, get_CAPPI_xy

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
    nrays : int
        Number of rays in the volume.
    nsweeps : int
        Number of sweep in the volume.

    """

    def __init__(self, fields,  scan_type, time, range, azimuth, elevation,latitude,
                 longitude, altitude, sweep_start_ray_index, sweep_end_ray_index,
                 fixed_angle, bins_per_sweep, nyquist_velocity, frequency, unambiguous_range,
                 nrays, nsweeps, sitename, pyart_radar=None):
        super(PRD, self).__init__()
        keys = fields.keys()
        self.fields = []
        for idx, (istart, iend) in enumerate(zip(sweep_start_ray_index, sweep_end_ray_index)):
            x, y, z = antenna_vectors_to_cartesian_cwr(range[:bins_per_sweep[idx]], azimuth[istart:iend+1],\
                                                   elevation[istart:iend+1], altitude)
            lon, lat = cartesian_to_geographic_aeqd(x, y, longitude, latitude)
            isweep_data = xr.Dataset(coords={'azimuth': (['time', ], azimuth[istart:iend+1]),
                                            'elevation': (['time',], elevation[istart:iend+1]),
                                             'x':(['time','range'], x),
                                             'y':(['time', 'range'], y),
                                             'z':(['time', 'range'], z+altitude),
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
                         "nyquist_velocity":(['sweep',], nyquist_velocity),
                        "unambiguous_range":(['sweep',], unambiguous_range),
                        "rays_per_sweep": (['sweep',], sweep_end_ray_index-sweep_start_ray_index+1),
                        "fixed_angle": (["sweep",], fixed_angle),
                        "beam_width":(["sweep",], 360./(sweep_end_ray_index-sweep_start_ray_index+1))},
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
        self.scan_info['start_time'].attrs = DEFAULT_METADATA['start_time']
        self.scan_info['end_time'].attrs = DEFAULT_METADATA['end_time']
        self.scan_info['beam_width'].attrs = DEFAULT_METADATA['beam_width']
        self.nsweeps = nsweeps
        self.nrays = nrays
        self.sitename = sitename
        self.get_vol_data()
        self.product = xr.Dataset()
        self.PyartRadar = pyart_radar

    def ToPyartRadar(self):
        return self.PyartRadar

    def ordered_az(self, inplace=False):
        """
        regrid radar object by azimuth
        :return:
        """
        if inplace:
            for isweep in self.scan_info.sweep.values:
                self.fields[isweep] = self.fields[isweep].swap_dims({"time": "azimuth"}).sortby("azimuth") ##对数据重新排序
            return None
        else:
            prd_dat = PRD_AZ()
            prd_dat.scan_info = self.scan_info
            for isweep in self.scan_info.sweep.values:
                prd_dat.fields.append(self.fields[isweep].swap_dims({"time": "azimuth"}).sortby("azimuth"))
            return prd_dat

    def add_product_CR_xy(self, XRange, YRange):
        """
        计算给定范围的组合反射率
        :param XRange: np.ndarray, 1d, units:meters
        :param YRange: np.ndarray, 1d, units:meters
        :return:
        """
        GridX, GridY = np.meshgrid(XRange, YRange, indexing="ij")
        vol_azimuth, vol_range, fix_elevation, vol_value, radar_height,\
        radar_lon_0, radar_lat_0 = self.vol
        fillvalue = -999.
        GridV = get_CR_xy(vol_azimuth, vol_range, fix_elevation, vol_value,\
                          radar_height, GridX.astype(np.float64), GridY.astype(np.float64), -999.)
        self.product.coords["x_cr"] = XRange
        self.product.coords["y_cr"] = YRange
        self.product["CR"] = (('x_cr', 'y_cr'), np.where(GridV==fillvalue, np.nan, GridV))
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
        self.product["CR"].attrs = {'units': 'dBZ',
                                    'standard_name': 'Composite_reflectivity_factor',
                                    'long_name': 'Composite_reflectivity_factor',
                                    'axis': 'xy_coordinate',
                                    'comment': 'Maximum reflectance of all level',}

    def add_product_CAPPI_xy(self, XRange, YRange, level_height):
        """
        计算给定范围的CAPPI的图像
        :param XRange: np.ndarray, 1d, units:meters
        :param YRange: np.ndarray, 1d, units:meters
        :param level_height: 要插值的高度，常量, units:meters
        :return:
        """
        GridX, GridY = np.meshgrid(XRange, YRange, indexing="ij")
        vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, \
        radar_lon_0, radar_lat_0 = self.vol
        fillvalue = -999.
        GridV = get_CAPPI_xy(vol_azimuth, vol_range, fix_elevation, vol_value, radar_height,
                             GridX.astype(np.float64), GridY.astype(np.float64), level_height, fillvalue)
        self.product.coords["x_cappi_%d"%level_height] = XRange
        self.product.coords["y_cappi_%d"%level_height] = YRange
        self.product["CAPPI_%d"%level_height] = (("x_cappi_%d"%level_height, "y_cappi_%d"%level_height),
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
        self.product["CAPPI_%d"%level_height].attrs = {'units': 'dBZ',
                                    'standard_name': 'Constant_altitude_plan_position_indicator',
                                    'long_name': 'Constant_altitude_plan_position_indicator',
                                    'axis': 'xy_coordinate',
                                    'comment': 'CAPPI of level %d m.'%level_height, }

    def add_product_CR_lonlat(self, XLon, YLat):
        """
        计算给定经纬度范围的组合反射率
        :param XLon:np.ndarray, 1d, units:degree
        :param YLat:np.ndarray, 1d, units:degree
        :return:
        """
        fillvalue = -999.
        GridLon, GridLat = np.meshgrid(XLon, YLat, indexing="ij")
        projparams = {"proj": "aeqd", "lon_0": self.scan_info["longitude"].values,
                      "lat_0": self.scan_info["latitude"].values}
        proj = pyproj.Proj(projparams)
        GridX, GridY = proj(GridLon, GridLat, inverse=False)
        vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, \
        radar_lon_0, radar_lat_0 = self.vol
        GridV = get_CR_xy(vol_azimuth, vol_range, fix_elevation, vol_value, \
                          radar_height, GridX.astype(np.float64), GridY.astype(np.float64), -999.)
        self.product.coords["lon_cr"] = XLon
        self.product.coords["lat_cr"] = YLat
        self.product["CR_geo"] = (('lon_cr', 'lat_cr'), np.where(GridV == fillvalue, np.nan, GridV))
        self.product.coords["lon_cr"].attrs = {'units': 'degrees',
                                             'standard_name': 'CR_product_lon_axis ',
                                             'long_name': 'longitude_cr',
                                             'axis': 'lonlat_coordinate'}
        self.product.coords["lat_cr"].attrs = {'units': 'degrees',
                                             'standard_name': 'CR_product_lat_axis ',
                                             'long_name': 'latitude_cr',
                                             'axis': 'lonlat_coordinate'}
        self.product["CR_geo"].attrs = {'units': 'dBZ',
                                    'standard_name': 'Composite_reflectivity_factor_lonlat_grid',
                                    'long_name': 'Composite_reflectivity_factor',
                                    'axis': 'lonlat_coordinate',
                                    'comment': 'Maximum reflectance of all level', }

    def add_product_CAPPI_lonlat(self, XLon, YLat, level_height):
        """
        计算给定经纬度范围的CAPPI
        :param XLon:np.ndarray, 1d, units:degrees
        :param YLat:np.ndarray, 1d, units:degrees
        :param level_height:常量，要计算的高度
        :return:
        """
        fillvalue = -999.
        GridLon, GridLat = np.meshgrid(XLon, YLat, indexing="ij")
        projparams = {"proj": "aeqd", "lon_0": self.scan_info["longitude"].values,
                      "lat_0": self.scan_info["latitude"].values}
        proj = pyproj.Proj(projparams)
        GridX, GridY = proj(GridLon, GridLat, inverse=False)
        vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, \
        radar_lon_0, radar_lat_0 = self.vol
        GridV = get_CAPPI_xy(vol_azimuth, vol_range, fix_elevation, vol_value, radar_height,
                             GridX.astype(np.float64), GridY.astype(np.float64), level_height, fillvalue)
        self.product.coords["lon_cappi_%d" % level_height] = XLon
        self.product.coords["lat_cappi_%d" % level_height] = YLat
        self.product["CAPPI_geo_%d" % level_height] = (("lon_cappi_%d" % level_height, "lat_cappi_%d" % level_height),
                                                   np.where(GridV == fillvalue, np.nan, GridV))
        self.product.coords["lon_cappi_%d" % level_height].attrs = {'units': 'degrees',
                                                                  'standard_name': 'CAPPI_product_lon_axis ',
                                                                  'long_name': 'longitude_CAPPI',
                                                                  'axis': 'lonlat_coordinate',}
        self.product.coords["lat_cappi_%d" % level_height].attrs = {'units': 'degrees',
                                                                  'standard_name': 'CAPPI_product_lat_axis ',
                                                                  'long_name': 'latitude_CAPPI',
                                                                  'axis': 'lonlat_coordinate',}
        self.product["CAPPI_geo_%d" % level_height].attrs = {'units': 'dBZ',
                                                         'standard_name': 'Constant_altitude_plan_position_indicator',
                                                         'long_name': 'Constant_altitude_plan_position_indicator',
                                                         'axis': 'lonlat_coordinate',
                                                         'comment': 'CAPPI of level %d m' % level_height, }

    def get_vol_data(self, field_name="dBZ", fillvalue=-999.):
        """
        获取用于插值的雷达体扫数据
        :return:
        """
        order_vol = self.ordered_az()
        order_vol.fields = [order_vol.fields[i] for i in order_vol.scan_info["fixed_angle"].argsort().values]
        order_vol.scan_info["fixed_angle"] = order_vol.scan_info["fixed_angle"].sortby("sweep")
        vol_azimuth = [ppi.azimuth.values for ppi in order_vol.fields]
        vol_range = [ppi.range.values for ppi in order_vol.fields]
        fix_elevation = order_vol.scan_info["fixed_angle"].values
        vol_value = [np.where(np.isnan(ppi[field_name].values), fillvalue, ppi[field_name].values).astype(np.float64) for ppi in order_vol.fields]
        radar_height = float(order_vol.scan_info["altitude"].values)
        radar_lon_0 = float(order_vol.scan_info["longitude"].values)
        radar_lat_0 = float(order_vol.scan_info["latitude"].values)
        self.vol = vol_azimuth, vol_range, fix_elevation.astype(np.float64), vol_value, radar_height, radar_lon_0, radar_lat_0

    def get_RHI_data(self, az, field_name="dBZ"):
        """
        获取RHI剖面数据
        :param az:
        :param field_name:
        :return:
        """
        mesh_RHI = []
        mesh_RANGE = []
        mesh_Z = []
        order_dat = self.ordered_az()
        for isweep in self.scan_info.sweep.values:
            if (isweep>0) and (self.scan_info.fixed_angle.values[isweep] < self.scan_info.fixed_angle.values[isweep-1]):
                continue  ##remove VCP26 Type data
            isweep_data = order_dat.fields[isweep].sel(azimuth=az, method='nearest')[field_name]
            x, y, z = antenna_vectors_to_cartesian_rhi(isweep_data.range, isweep_data.azimuth,\
                                                       isweep_data.elevation, self.scan_info.altitude.values)
            mesh_xy = np.sqrt(x**2 + y**2)
            mesh_RHI.append(isweep_data.values.reshape(1,-1))
            mesh_RANGE.append(mesh_xy)
            mesh_Z.append(z)
        return mesh_RANGE, mesh_Z, mesh_RHI

    def get_vcs_data(self, start_point, end_point, field_name):
        """
        :param start_point:
        :param end_point:
        :param field_name:
        :return:
        """
        order_dat = self.ordered_az() ##求排序后的数据
        start_x, start_y = start_point
        end_x, end_y = end_point
        bins_res = (self.fields[0].range[1] - self.fields[0].range[0]).values
        start_end_dis = np.sqrt((start_x - end_x) ** 2 + (start_y - end_y) ** 2)
        npoints = int(start_end_dis/bins_res + 1)
        x_line = np.linspace(start_x, end_x, npoints)
        y_line = np.linspace(start_y, end_y, npoints)
        xy_line_1d = np.linspace(0, start_end_dis, npoints)
        xy = np.stack([xy_line_1d, xy_line_1d], axis=0)
        mesh_Z = []
        mesh_vcs = []
        mesh_xy = []
        # 先对剖线取最邻近点
        for isweep, ifield in enumerate(order_dat.fields):
            if (isweep>0) and (self.scan_info.fixed_angle.values[isweep] < self.scan_info.fixed_angle.values[isweep-1]):
                continue  ##remove VCP26 Type data
            az, ranges, _ = cartesian_to_antenna_cwr(x_line, y_line, self.scan_info.fixed_angle.values[isweep],\
                                                     self.scan_info.altitude.values)
            vcs_data = ifield[field_name].sel(azimuth=xr.DataArray(az, dims="vcs_r"),\
                                              range=xr.DataArray(ranges, dims="vcs_r"), method="nearest") ##选取插值的点
            _, _, z = antenna_vectors_to_cartesian_vcs(vcs_data.range, vcs_data.azimuth, vcs_data.elevation,\
                                             order_dat.scan_info.altitude.values, \
                                             self.scan_info.beam_width.values[isweep])
            mesh_xy.append(xy)
            mesh_Z.append(z)
            mesh_vcs.append(vcs_data.values.reshape(1,-1))
        return mesh_xy, mesh_Z, mesh_vcs

class PRD_AZ:
    """
    data obj for radar data, AZ as dims!
    """
    def __init__(self):
        self.scan_info = None
        self.fields = []
        self.product = xr.Dataset()





