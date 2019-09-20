# -*- coding: utf-8 -*-
"""
为了适应中国雷达在不同仰角的探测距离不同以及前几层仰角 dop和ref分开扫描的问题
提出NuistRadar Object，以方便后续的算法及绘图
"""
import sys
import numpy as np
import xarray as xr
from default_config import DEFAULT_METADATA, FILL_VALUE

class NuistRadar(object):
    """
    A class for storing antenna coordinate radar data.
    Attributes
    ----------
    fields : List, nsweep,
        (time/elevation/azimuth, bins)
        Moment fields.
    metadata : dict ##描述站点和仪器的一些信息
        Metadata describing the instrument and data.
    scan_type : str
        Type of scan, one of 'ppi', 'rhi', 'sector' or 'other'. If the scan
        volume contains multiple sweep modes this should be 'other'.
    latitude : DataArray
        Latitude of the instrument.
    longitude: DataArray
        Longitude of the instrument.
    altitude : DataArray
        Altitude of the instrument, above sea level.
    fixed_angle : DataArray
        Target angle for thr sweep. Azimuth angle in RHI modes, elevation
        angle in all other modes.
    rays_per_sweep : DataArray
        Number of rays in each sweep. The data key of this attribute is
        create upon first access from the data in the sweep_start_ray_index and
        sweep_end_ray_index attributes. If the sweep locations needs to be
        modified, do this prior to accessing this attribute or use
        :py:func:`init_rays_per_sweep` to reset the attribute.
    bins_per_sweep : DataArray     !!!##added
        Number of bins in each sweep. The data key of this attribute is
        create upon first access from the data in the sweep_start_ray_index and
        sweep_end_ray_index attributes. If the sweep locations needs to be
        modified, do this prior to accessing this attribute or use
        :py:func:`init_rays_per_sweep` to reset the attribute.
    instrument_parameters : dict{'nyquist_velocity', 'unambiguous_range', 'frequency'}
        Instrument parameters, if not provided this attribute is set to None,
        indicating these parameters are not avaiable. This dictionary also
        includes variables in the radar_parameters CF/Radial subconvention.
    nrays : int
        Number of rays in the volume.
    nsweeps : int
        Number of sweep in the volume.
    """

    def __init__(self):
        super(NuistRadar, self).__init__()


