from libc.math cimport sin, cos, asin, acos, tan, atan
import numpy as np
cimport numpy as cnp
import cython

def antenna_to_cartesian(double ranges, double azimuth, double elevation, double h):
    """
    将天线坐标系转换为笛卡尔坐标系
    """
    cdef double PI = 3.141592653589793
    cdef double R = 8494666.6666666661
    cdef double theta_a = azimuth/180.*PI
    cdef double theta_e = elevation/180.*PI
    cdef double x, y, z, s
    z = ((ranges * cos(theta_e))**2 + (R + h + ranges*sin(theta_e))**2)**0.5 - R
    s = R * asin(ranges*cos(theta_e)/(R+z))
    x = s * sin(theta_a)
    y = s * cos(theta_a)
    return x, y, z

def xye_to_antenna(double x, double y, double elevation, double h):
    """
    将直角坐标系和仰角转换为天线坐标系
    """
    cdef double PI = 3.141592653589793
    cdef double R = 8494666.6666666661
    cdef double s, theta_e, azimuth, ranges, z
    s = (x**2 + y**2)**0.5
    theta_e = elevation/180.*PI
    ranges = tan(s/R) * (R+h)/cos(theta_e)
    z = (R+h)/cos(theta_e + s/R) * cos(theta_e) - R
    azimuth = xy_to_azimuth(x, y)
    return azimuth, ranges, z

def cartesian_to_antenna(double x, double y, double z, double h):
    """
    直角坐标系转换为天线坐标系
    ranges: 天线坐标系距离雷达天线的距离
    """
    cdef double PI = 3.141592653589793
    cdef double R = 8494666.6666666661
    cdef double elevation, ranges, azimuth
    ranges = ((R+h)**2 + (R+z)**2 - 2*(R+h)*(R+z)*cos((x**2 + y**2)**0.5/R))**0.5
    elevation = acos((R+z)*sin((x**2 + y**2)**0.5/R)/ranges) * 180./PI
    azimuth = xy_to_azimuth(x, y)
    return azimuth, ranges, elevation

@cython.cdivision(True)
def xy_to_azimuth(double x, double y):
    """
    using x and y to cal azimuth
    input
    x : units : meters
    y : units : meters
    return
    azimuth: units:degree
    """
    cdef double PI = 3.141592653589793
    cdef double az = PI/2.0 - atan(y/x)
    if az >= 0:
        az = az * 180. / PI
    else:
        az = (2*PI + az) * 180./PI
    if x<0:
        az = 180 + az
    return az

def interp_ppi(double az, double r, double az_0, double az_1, double r_0, double r_1, double mat_00, double mat_01, double mat_10, double mat_11, double fillvalue):
    """
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
    cdef double interped
    if ((mat_00 != fillvalue) and (mat_01 != fillvalue)) and ((mat_10 != fillvalue) and(mat_11 != fillvalue)):
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

@cython.boundscheck(False)  # Deactivate bounds checking
def ppi_to_grid(cnp.ndarray[cnp.float64_t, ndim=1] azimuth, cnp.ndarray[cnp.float64_t, ndim=1] ranges,
                double elevation, cnp.ndarray[cnp.float64_t, ndim=2] mat_ppi, double radar_height,
                cnp.ndarray[cnp.float64_t, ndim=2] GridX, cnp.ndarray[cnp.float64_t, ndim=2] GridY,
                double fillvalue):
    """
    将PPI扫描格点化
    :param azimuth:mat_ppi第一个维度对应的方位角, np.ndarray (1d), units:degree
    :param ranges:mat_ppi第二个维度对应的斜距, np.ndarray (1d), units:meters
    :param elevation:该层ppi扫描的仰角, const, 常量, units:degree
    :param mat_ppi:待插值的格点数据， np.ndarray(2d), units: //
    :param radar_height:雷达距离海平面的高度，const，常量, units:meters
    :param GridX:待插值的二维格点，X坐标, np.ndarray(2d), units:meters
    :param GridY:待插值的二维格点，Y坐标, np.ndarray(2d), units:meters
    :param fillvalue:缺测值对应于mat_ppi
    :return:
    """
    cdef int naz = azimuth.shape[0]
    cdef int nrange = ranges.shape[0]
    cdef int Nx = GridX.shape[0]
    cdef int Ny = GridY.shape[1]
    cdef cnp.ndarray GridValue = np.zeros([Nx, Ny], dtype=np.float64)
    cdef double az, r, z, az_last
    cdef int ix, iy, iaz, ir
    for ix in range(Nx):
        for iy in range(Ny):
            az, r, z  = xye_to_antenna(GridX[ix, iy], GridY[ix, iy], elevation, radar_height)
            if r > ranges[-1]:
                GridValue[ix, iy] = fillvalue
            elif r < ranges[0]:
                GridValue[ix, iy] = fillvalue
            else:
                for iaz in range(naz):
                    if az < azimuth[iaz]:
                        break
                else:
                    iaz = 0
                    az = az - 360.

                if iaz == 0:
                    az_last = azimuth[iaz-1] - 360.
                else:
                    az_last = azimuth[iaz-1]

                for ir in range(1, nrange):
                    if r < ranges[ir]:
                        break
                GridValue[ix, iy] = interp_ppi(az, r, az_last, azimuth[iaz], ranges[ir-1], ranges[ir], mat_ppi[iaz-1, ir-1], mat_ppi[iaz-1, ir], mat_ppi[iaz, ir-1], mat_ppi[iaz, ir], fillvalue)
    return GridValue

@cython.boundscheck(False)  # Deactivate bounds checking
def get_CAPPI_xy(vol_azimuth, vol_range, cnp.ndarray[cnp.float64_t, ndim=1] fix_elevation, vol_value,
                 double radar_height, cnp.ndarray[cnp.float64_t, ndim=2] GridX,
                 cnp.ndarray[cnp.float64_t, ndim=2] GridY, double level_height,  double fillvalue):
    """
    由雷达体扫数据，插值CAPPI图像
    :param vol_azimuth:存放多个仰角体扫方位角的列表, list, units:degree
    :param vol_range:存放多个仰角体扫距离的列表, list, units:meters
    :param fix_elevation:每个仰角体扫对应的仰角， np.ndarray， 1d
    :param vol_value:存放多个仰角体扫数据的列表, list
    :param radar_height:常量, 雷达距离海平面的高度， units:meters
    :param GridX:要插值的二维格点X, np.ndarray, 2d, units:meters
    :param GridY:要插值的二维格点Y, np.ndarray, 2d, units:meters
    :param level_height:常量，待插值的高度， units:meters
    :param fillvalue:常量，缺测值
    :return:
    """
    cdef int Ne = fix_elevation.shape[0]
    cdef int Nx = GridX.shape[0]
    cdef int Ny = GridX.shape[1]
    cdef cnp.ndarray GridValue = np.full([Nx, Ny], fillvalue, dtype=np.float64)
    cdef int ix, iy, ie, iaz_0, iaz_1, range_0, range_1
    cdef double az, el, r, temp_az, temp_range, az_last_0, az_last_1
    cdef double ER00, ER01, ER10, ER11, IER0, IER1
    for ix in range(Nx):
        for iy in range(Ny):
            az, r, el = cartesian_to_antenna(GridX[ix, iy], GridY[ix, iy], level_height, radar_height)
            if (el <= fix_elevation[-1]) and (el >= fix_elevation[0]):
                for ie in range(Ne):
                    if el < fix_elevation[ie]:
                        break
                for iaz_0, temp_az in enumerate(vol_azimuth[ie-1]):
                    if az < temp_az:
                        break
                else:
                    iaz_0 = 0
                    az = az - 360.

                for range_0, temp_range in enumerate(vol_range[ie-1]):
                    if r < temp_range:
                        break

                if iaz_0 == 0:
                    az_last_0 = vol_azimuth[ie-1][iaz_0-1] - 360.
                else:
                    az_last_0 = vol_azimuth[ie-1][iaz_0-1]

                for iaz_1, temp_az in enumerate(vol_azimuth[ie]):
                    if az < temp_az:
                        break
                else:
                    iaz_1 = 0
                for range_1, temp_range in enumerate(vol_range[ie]):
                    if r < temp_range:
                        break
                if iaz_1 == 0:
                    az_last_1 = vol_azimuth[ie][iaz_1-1] - 360.
                else:
                    az_last_1 = vol_azimuth[ie][iaz_1-1]

                if (r <= vol_range[ie][-1]) and (r <= vol_range[ie-1][-1]):
                    ER00 = interp_azimuth(az, az_last_0, vol_azimuth[ie-1][iaz_0], vol_value[ie-1][iaz_0-1, range_0-1],
                                          vol_value[ie-1][iaz_0, range_0-1], fillvalue)
                    ER01 = interp_azimuth(az, az_last_0, vol_azimuth[ie-1][iaz_0], vol_value[ie-1][iaz_0-1, range_0],
                                          vol_value[ie-1][iaz_0, range_0], fillvalue)
                    ER10 = interp_azimuth(az, az_last_1, vol_azimuth[ie][iaz_1], vol_value[ie][iaz_1-1, range_1-1],
                                          vol_value[ie][iaz_1, range_1-1], fillvalue)
                    ER11 = interp_azimuth(az, az_last_1, vol_azimuth[ie][iaz_1], vol_value[ie][iaz_1-1, range_1],
                                          vol_value[ie][iaz_1, range_1], fillvalue)
                    IER0 = interp_azimuth(r, vol_range[ie-1][range_0-1], vol_range[ie-1][range_0], ER00, ER01, fillvalue)
                    IER1 = interp_azimuth(r, vol_range[ie][range_1-1], vol_range[ie][range_1], ER10, ER11, fillvalue)
                    GridValue[ix, iy] = interp_azimuth(el, fix_elevation[ie-1], fix_elevation[ie], IER0, IER1, fillvalue)
    return GridValue

def interp_azimuth(double az, double az_0, double az_1, double dat_0, double dat_1, double fillvalue):
    """
    在两个方位角或者距离之间进行插值
    """
    if (dat_0 != fillvalue) and (dat_1 != fillvalue):
        return ((az_1 - az)*dat_0 + (az - az_0) * dat_1)/(az_1 - az_0)
    elif dat_0 == fillvalue:
        return dat_1
    else:
        return dat_0

@cython.boundscheck(False)  # Deactivate bounds checking
def get_CR_xy(vol_azimuth, vol_range, cnp.ndarray[cnp.float64_t, ndim=1] fix_elevation, vol_value, double radar_height,
              cnp.ndarray[cnp.float64_t, ndim=2] GridX, cnp.ndarray[cnp.float64_t, ndim=2] GridY, double fillvalue):
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
    cdef int Ne = fix_elevation.shape[0]
    cdef int Nx = GridX.shape[0]
    cdef int Ny = GridY.shape[1]
    cdef cnp.ndarray GridValue = np.zeros([Ne, Nx, Ny], dtype=np.float64)
    cdef int ie
    for ie in range(Ne):
        GridValue[ie,:,:] = ppi_to_grid(vol_azimuth[ie], vol_range[ie], fix_elevation[ie],
                                       vol_value[ie], radar_height, GridX, GridY, fillvalue)

    GridValue = np.where(GridValue == fillvalue, np.nan, GridValue)
    GridValue = np.nanmax(GridValue, axis=0)
    GridValue = np.where(np.isnan(GridValue), fillvalue, GridValue)
    return GridValue