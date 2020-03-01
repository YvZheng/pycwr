import numpy as np
from .transforms import cartesian_to_antenna_cwr, cartesian_xyz_to_antenna

def interp_ppi(az, r, az_0, az_1, r_0, r_1, mat_00, mat_01, mat_10, mat_11, fillvalue=-999.):
    """
    利用雷达扫描的周围四个点插值中间的点(az, r)
    interp radar ppi scan data
    az : target azimuth, units:degree
    r : target range, units:meters
    az_0 : grid start azimuth, units:degree
    az_1 : grid end azimuth, units:degree
    r_0 : grid start range , units : meters
    r_1 : grid end range, units: meters
    mat_00: data for [az_0, r_0]
    mat_01: data for [az_0, r_1]
    mat_10: data for [az_1, r_0]
    mat_11: data for [az_1, r_1]
    fillvalue: fillvalue for mat
    return target value interped, units: like mat
    """
    if ((mat_00 != fillvalue) and (mat_01 != fillvalue)) and ((mat_10 != fillvalue) and (mat_11 != fillvalue)):
        interped = (mat_00 * (az_1 - az) * (r_1 - r) + mat_10 * (az - az_0) * (r_1 - r) + mat_01 * (az_1 - az) * (r - r_0) + mat_11 * (az - az_0) * (r - r_0))/(r_1 - r_0)/(az_1 - az_0)
    elif (mat_00 != fillvalue) and (mat_01 != fillvalue):
        interped = (mat_00 * (r_1 - r) + mat_01 * (r - r_0))/(r_1 - r_0)
    elif ((mat_10 != fillvalue) and (mat_11 != fillvalue)):
        interped = (mat_10 * (r_1 - r) + mat_11 * (r - r_0))/(r_1 - r_0)
    elif ((mat_00 != fillvalue) and (mat_10 != fillvalue)):
        interped = (mat_00 * (az_1 - az) + mat_10 * (az - az_0))/(az_1 - az_0)
    elif ((mat_01 != fillvalue) and (mat_11 != fillvalue)):
        interped = (mat_01 * (az_1 - az) + mat_11 * (az - az_0))/(az_1 - az_0)
    else:
        interped = fillvalue
    return interped

def interp_azimuth(az, az_0, az_1, dat_0, dat_1, fillvalue=-999.):
    """
    在两个方位角或者距离之间进行插值
    """
    if (dat_0 != fillvalue) and (dat_1 != fillvalue):
        return ((az_1 - az)*dat_0 + (az - az_0) * dat_1)/(az_1 - az_0)
    elif dat_0 == fillvalue:
        return dat_1
    else:
        return dat_0

def ppi_to_grid(azimuth, ranges, elevation, mat_ppi, radar_height, GridX, GridY, fillvalue=-999.):
    """
    将ppi转化为grid
    :param azimuth: mat_ppi第一个维度对应的方位角, np.ndarray (1d), units:degree
    :param ranges: mat_ppi第二个维度对应的斜距, np.ndarray (1d), units:meters
    :param elevation: 该层ppi扫描的仰角, const, 常量, units:degree
    :param mat_ppi: 待插值的格点数据， np.ndarray(2d), units: //
    :param radar_height: 雷达距离海平面的高度，const，常量, units:meters
    :param GridX: 待插值的二维格点，X坐标, np.ndarray(2d), units:meters
    :param GridY: 待插值的二维格点，Y坐标, np.ndarray(2d), units:meters
    :param fillvalue: 缺测值对应于mat_ppi
    :return:  GridValue(插值后产品的结果)
    """
    assert np.dtype(np.float64).type == GridX.dtype, "GridX 的类型须是np.float64， 以免超出最大范围"
    Grid_az, Grid_range, Grid_Z = cartesian_to_antenna_cwr(GridX, GridY, elevation, radar_height)
    GridValue = np.full_like(GridX, fillvalue, dtype=np.float64)
    Nx, Ny = GridX.shape
    Naz, Nr = mat_ppi.shape
    for ix in range(Nx):
        for iy in range(Ny):
            if Grid_range[ix, iy] <= ranges[-1] and Grid_range[ix, iy] >= ranges[0]:
                for iaz in range(Naz):
                    if Grid_az[ix, iy] < azimuth[iaz]:
                        break
                else:
                    iaz = 0
                    Grid_az[ix, iy] = Grid_az[ix, iy] - 360.
                if iaz == 0:
                    az_last = azimuth[iaz - 1] - 360.
                else:
                    az_last = azimuth[iaz - 1]

                for ir in range(1, Nr):
                    if Grid_range[ix, iy] < ranges[ir]:
                        break
                GridValue[ix, iy] = interp_ppi(Grid_az[ix, iy], Grid_range[ix,iy], az_last, azimuth[iaz],
                                               ranges[ir - 1], ranges[ir], mat_ppi[iaz - 1, ir - 1],
                                               mat_ppi[iaz - 1, ir], mat_ppi[iaz, ir - 1], mat_ppi[iaz, ir],
                                               fillvalue)
    return GridX, GridY, GridValue

def get_CR_xy(vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, GridX, GridY, fillvalue=-999.):
    """
    计算组合反射率，利用雷达体扫的数据
    :param vol_azimuth:存放多个仰角体扫方位角的列表, list, units:degree
    :param vol_range:存放多个仰角体扫距离的列表, list, units:meters
    :param fix_elevation:每个仰角体扫对应的仰角， np.ndarray， 1d
    :param vol_value: 存放多个仰角体扫数据的列表, list
    :param radar_height: 常量, 雷达距离海平面的高度， units:meters
    :param GridX: 组合反射率的二维格点的X的值, units:meters
    :param GridY: 组合反射率的二维格点的Y的值, units:meters
    :param fillvalue: 缺测值
    :return:
    """
    Ne = fix_elevation.shape[0]
    Nx, Ny = GridX.shape
    GridValue = np.zeros([Ne, Nx, Ny], dtype=np.float64)
    for ie in range(Ne):
        GridX, GridY, GridValue[ie] = ppi_to_grid(vol_azimuth[ie], vol_range[ie], fix_elevation[ie], vol_value[ie],
                                                  radar_height, GridX, GridY, fillvalue)
    GridValue = np.nanmax(np.where(GridValue == fillvalue, np.nan, GridValue), axis=0)
    return np.where(np.isnan(GridValue), fillvalue, GridValue)

def get_CAPPI_xy(vol_azimuth, vol_range, fix_elevation, vol_value, radar_height,
                 GridX, GridY, level_height, fillvalue=-999.):
    """
    由雷达体扫数据，插值CAPPI
    :param vol_azimuth:存放多个仰角体扫方位角的列表, list, units:degree
    :param vol_range:存放多个仰角体扫距离的列表, list, units:meters
    :param fix_elevation:每个仰角体扫对应的仰角， np.ndarray， 1d
    :param vol_value: 存放多个仰角体扫数据的列表, list
    :param radar_height:常量, 雷达距离海平面的高度， units:meters
    :param GridX:要插值的二维格点X, np.ndarray, 2d, units:meters
    :param GridY:要插值的二维格点Y, np.ndarray, 2d, units:meters
    :param level_height:常量，待插值的高度， units:meters
    :param fillvalue:常量，缺测值
    :return:
    """
    Ne = fix_elevation.shape[0]
    Nx, Ny = GridX.shape
    GridValue = np.full_like(GridX, fillvalue,dtype=np.float64)
    Grid_az, Grid_range, Grid_el = cartesian_xyz_to_antenna(GridX, GridY, level_height, radar_height)
    for ix in range(Nx):
        for iy in range(Ny):
            if (Grid_el[ix, iy] >= fix_elevation[0]) and (Grid_el[ix, iy] <= fix_elevation[-1]):
                for ie in range(Ne):
                    if Grid_el[ix, iy] < fix_elevation[ie]:
                        break
                for iaz_0, temp_az in enumerate(vol_azimuth[ie - 1]):
                    if Grid_az[ix, iy] < temp_az:
                        break
                else:
                    iaz_0 = 0

                for range_0, temp_range in enumerate(vol_range[ie - 1]):
                    if Grid_range[ix, iy] < temp_range:
                        break

                if iaz_0 == 0:
                    az_last_0 = vol_azimuth[ie - 1][iaz_0 - 1] - 360.
                else:
                    az_last_0 = vol_azimuth[ie - 1][iaz_0 - 1]

                for iaz_1, temp_az in enumerate(vol_azimuth[ie]):
                    if Grid_az[ix, iy] < temp_az:
                        break
                else:
                    iaz_1 = 0

                for range_1, temp_range in enumerate(vol_range[ie]):
                    if Grid_range[ix, iy] < temp_range:
                        break

                if iaz_1 == 0:
                    az_last_1 = vol_azimuth[ie][iaz_1 - 1] - 360.
                else:
                    az_last_1 = vol_azimuth[ie][iaz_1 - 1]
                if (Grid_range[ix, iy] <= vol_range[ie][-1]) and (Grid_range[ix, iy] <= vol_range[ie - 1][-1]):
                    ER00 = interp_azimuth(Grid_az[ix, iy], az_last_0, vol_azimuth[ie - 1][iaz_0],
                                          vol_value[ie - 1][iaz_0 - 1, range_0 - 1],
                                          vol_value[ie - 1][iaz_0, range_0 - 1], fillvalue)
                    ER01 = interp_azimuth(Grid_az[ix, iy], az_last_0, vol_azimuth[ie - 1][iaz_0],
                                          vol_value[ie - 1][iaz_0 - 1, range_0], vol_value[ie - 1][iaz_0, range_0], fillvalue)
                    ER10 = interp_azimuth(Grid_az[ix, iy], az_last_1, vol_azimuth[ie][iaz_1], vol_value[ie][iaz_1 - 1, range_1 - 1],
                                          vol_value[ie][iaz_1, range_1 - 1], fillvalue)
                    ER11 = interp_azimuth(Grid_az[ix, iy], az_last_1, vol_azimuth[ie][iaz_1], vol_value[ie][iaz_1 - 1, range_1],
                                          vol_value[ie][iaz_1, range_1], fillvalue)
                    IER0 = interp_azimuth(Grid_range[ix, iy], vol_range[ie - 1][range_0 - 1], vol_range[ie - 1][range_0], ER00, ER01,
                                          fillvalue)
                    IER1 = interp_azimuth(Grid_range[ix, iy], vol_range[ie][range_1 - 1], vol_range[ie][range_1], ER10, ER11, fillvalue)
                    GridValue[ix, iy] = interp_azimuth(Grid_el[ix, iy], fix_elevation[ie - 1], fix_elevation[ie], IER0, IER1, fillvalue)
    return GridValue
