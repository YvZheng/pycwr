from libc.math cimport sin, cos, asin, acos, tan, atan2, sqrt, fabs, NAN
import numpy as np
cimport numpy as cnp
import cython

cdef double PI = 3.141592653589793
cdef double DEG2RAD = PI / 180.0
cdef double RAD2DEG = 180.0 / PI
cdef double DEFAULT_EFFECTIVE_EARTH_RADIUS = 8494666.6666666661


cdef inline double _resolve_effective_earth_radius(object effective_earth_radius):
    cdef double radius
    if effective_earth_radius is None:
        return DEFAULT_EFFECTIVE_EARTH_RADIUS
    radius = float(effective_earth_radius)
    if radius <= 0.0:
        raise ValueError("effective_earth_radius must be positive")
    return radius


cdef inline double _clip_unit(double value) noexcept nogil:
    if value > 1.0:
        return 1.0
    if value < -1.0:
        return -1.0
    return value


cdef inline double _xy_norm(double x, double y) noexcept nogil:
    return sqrt(x * x + y * y)


cdef inline double _azimuth_from_xy(double x, double y) noexcept nogil:
    cdef double az = PI / 2.0 - atan2(y, x)
    if az < 0.0:
        az += 2.0 * PI
    return az * RAD2DEG


cdef inline int _search_right_1d(cnp.float64_t[:] values, double target) noexcept nogil:
    cdef Py_ssize_t lo = 0
    cdef Py_ssize_t hi = values.shape[0]
    cdef Py_ssize_t mid
    while lo < hi:
        mid = (lo + hi) // 2
        if target < values[mid]:
            hi = mid
        else:
            lo = mid + 1
    return <int>lo


cdef inline int _resolve_blind_code(object blind_method):
    cdef str method = "mask" if blind_method is None else str(blind_method).lower()
    if method == "mask":
        return 0
    if method == "nearest_gate":
        return 1
    if method == "lowest_sweep":
        return 2
    if method == "hybrid":
        return 3
    if method == "nearest_valid":
        return 4
    raise ValueError(
        "blind_method must be one of 'mask', 'nearest_gate', 'lowest_sweep', "
        "'hybrid', or 'nearest_valid'."
    )


cdef inline void _antenna_to_cartesian_impl(
    double ranges,
    double azimuth,
    double elevation,
    double h,
    double effective_earth_radius,
    double* x,
    double* y,
    double* z,
) noexcept nogil:
    cdef double theta_a = azimuth * DEG2RAD
    cdef double theta_e = elevation * DEG2RAD
    cdef double sin_e = sin(theta_e)
    cdef double cos_e = cos(theta_e)
    cdef double radar_radius = effective_earth_radius + h
    cdef double arc_arg
    cdef double s

    z[0] = sqrt(ranges * ranges + radar_radius * radar_radius +
                2.0 * ranges * radar_radius * sin_e) - effective_earth_radius
    arc_arg = _clip_unit(ranges * cos_e / (effective_earth_radius + z[0]))
    s = effective_earth_radius * asin(arc_arg)
    x[0] = s * sin(theta_a)
    y[0] = s * cos(theta_a)


cdef inline void _xye_to_antenna_impl(
    double x,
    double y,
    double elevation,
    double h,
    double effective_earth_radius,
    double* azimuth,
    double* ranges,
    double* z,
) noexcept nogil:
    cdef double s = _xy_norm(x, y)
    cdef double phi = s / effective_earth_radius
    cdef double theta_e = elevation * DEG2RAD
    cdef double sin_phi = sin(phi)
    cdef double cos_phi = cos(phi)
    cdef double denom = cos_phi - sin_phi * tan(theta_e)
    cdef double radar_radius = effective_earth_radius + h
    cdef double gate_radius
    cdef double ranges_sq

    azimuth[0] = _azimuth_from_xy(x, y)
    if denom == 0.0:
        ranges[0] = NAN
        z[0] = NAN
        return

    gate_radius = radar_radius / denom
    z[0] = gate_radius - effective_earth_radius
    ranges_sq = radar_radius * radar_radius + gate_radius * gate_radius - 2.0 * radar_radius * gate_radius * cos_phi
    if ranges_sq < 0.0:
        ranges_sq = 0.0
    ranges[0] = sqrt(ranges_sq)


cdef inline void _cartesian_to_antenna_impl(
    double x,
    double y,
    double z,
    double h,
    double effective_earth_radius,
    double* azimuth,
    double* ranges,
    double* elevation,
) noexcept nogil:
    cdef double s = _xy_norm(x, y)
    cdef double phi = s / effective_earth_radius
    cdef double radar_radius = effective_earth_radius + h
    cdef double gate_radius = effective_earth_radius + z
    cdef double ranges_sq = radar_radius * radar_radius + gate_radius * gate_radius - 2.0 * radar_radius * gate_radius * cos(phi)
    cdef double cos_arg

    if ranges_sq < 0.0:
        ranges_sq = 0.0
    ranges[0] = sqrt(ranges_sq)
    if ranges[0] == 0.0:
        elevation[0] = 0.0
    else:
        cos_arg = _clip_unit((radar_radius * radar_radius + ranges_sq - gate_radius * gate_radius) /
                             (2.0 * radar_radius * ranges[0]))
        elevation[0] = (acos(cos_arg) - PI / 2.0) * RAD2DEG
    azimuth[0] = _azimuth_from_xy(x, y)


def antenna_to_cartesian(double ranges, double azimuth, double elevation, double h, effective_earth_radius=None):
    """Convert antenna coordinates to Cartesian coordinates."""
    cdef double x = 0.0
    cdef double y = 0.0
    cdef double z = 0.0
    cdef double radius = _resolve_effective_earth_radius(effective_earth_radius)
    _antenna_to_cartesian_impl(ranges, azimuth, elevation, h, radius, &x, &y, &z)
    return x, y, z

def xye_to_antenna(double x, double y, double elevation, double h, effective_earth_radius=None):
    """Convert Cartesian x/y plus elevation back to antenna coordinates."""
    cdef double azimuth = 0.0
    cdef double ranges = 0.0
    cdef double z = 0.0
    cdef double radius = _resolve_effective_earth_radius(effective_earth_radius)
    _xye_to_antenna_impl(x, y, elevation, h, radius, &azimuth, &ranges, &z)
    return azimuth, ranges, z

def cartesian_to_antenna(double x, double y, double z, double h, effective_earth_radius=None):
    """Convert Cartesian coordinates back to antenna coordinates."""
    cdef double elevation = 0.0
    cdef double ranges = 0.0
    cdef double azimuth = 0.0
    cdef double radius = _resolve_effective_earth_radius(effective_earth_radius)
    _cartesian_to_antenna_impl(x, y, z, h, radius, &azimuth, &ranges, &elevation)
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
    return _azimuth_from_xy(x, y)

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
    cdef double az_span = az_1 - az_0
    cdef double r_span = r_1 - r_0
    if az_span == 0.0 or r_span == 0.0:
        return fillvalue
    if ((mat_00 != fillvalue) and (mat_01 != fillvalue)) and ((mat_10 != fillvalue) and(mat_11 != fillvalue)):
        interped = (mat_00 * (az_1 - az) * (r_1 - r) + mat_10 * (az - az_0) * (r_1 - r) + mat_01 * (az_1 - az) * (r - r_0) + mat_11 * (az - az_0) * (r - r_0)) / r_span / az_span
    elif (mat_00 != fillvalue) and (mat_01 != fillvalue):
        interped = (mat_00 * (r_1 - r) + mat_01 * (r - r_0)) / r_span
    elif ((mat_10 != fillvalue) and (mat_11 != fillvalue)):
        interped = (mat_10 * (r_1 - r) + mat_11 * (r - r_0)) / r_span
    elif ((mat_00 != fillvalue) and (mat_10 != fillvalue)):
        interped = (mat_00 * (az_1 - az) + mat_10 * (az - az_0)) / az_span
    elif ((mat_01 != fillvalue) and (mat_11 != fillvalue)):
        interped = (mat_01 * (az_1 - az) + mat_11 * (az - az_0)) / az_span
    else:
        interped = fillvalue
    return interped


cdef inline double _interp_linear(
    double coord,
    double coord_0,
    double coord_1,
    double dat_0,
    double dat_1,
    double fillvalue,
) noexcept nogil:
    if coord_1 == coord_0:
        return fillvalue
    if dat_0 == fillvalue and dat_1 == fillvalue:
        return fillvalue
    if dat_0 == fillvalue:
        return dat_1
    if dat_1 == fillvalue:
        return dat_0
    return ((coord_1 - coord) * dat_0 + (coord - coord_0) * dat_1) / (coord_1 - coord_0)


cdef inline double _interp_sweep_sample(
    cnp.float64_t[:] azimuth_view,
    cnp.float64_t[:] range_view,
    cnp.float64_t[:, :] value_view,
    double az,
    double r,
    double fillvalue,
    bint use_nearest_gate,
) noexcept nogil:
    cdef int naz = azimuth_view.shape[0]
    cdef int nrange = range_view.shape[0]
    cdef int iaz
    cdef int ir
    cdef double az_value = az
    cdef double az_last = 0.0
    cdef double r_value = r
    cdef double er0
    cdef double er1
    if naz < 2 or nrange < 2:
        return fillvalue
    if r_value > range_view[nrange - 1]:
        return fillvalue
    if r_value < range_view[0]:
        if not use_nearest_gate:
            return fillvalue
        r_value = range_view[0]

    iaz = _search_right_1d(azimuth_view, az_value)
    if iaz >= naz:
        iaz = 0
        az_value -= 360.0
    if iaz == 0:
        az_last = azimuth_view[naz - 1] - 360.0
    else:
        az_last = azimuth_view[iaz - 1]

    ir = _search_right_1d(range_view, r_value)
    if ir <= 0:
        ir = 1
    elif ir >= nrange:
        ir = nrange - 1

    er0 = _interp_linear(az_value, az_last, azimuth_view[iaz], value_view[iaz - 1, ir - 1], value_view[iaz, ir - 1], fillvalue)
    er1 = _interp_linear(az_value, az_last, azimuth_view[iaz], value_view[iaz - 1, ir], value_view[iaz, ir], fillvalue)
    return _interp_linear(r_value, range_view[ir - 1], range_view[ir], er0, er1, fillvalue)

@cython.boundscheck(False)
@cython.wraparound(False)
def ppi_to_grid(cnp.ndarray[cnp.float64_t, ndim=1] azimuth, cnp.ndarray[cnp.float64_t, ndim=1] ranges,
                double elevation, cnp.ndarray[cnp.float64_t, ndim=2] mat_ppi, double radar_height,
                cnp.ndarray[cnp.float64_t, ndim=2] GridX, cnp.ndarray[cnp.float64_t, ndim=2] GridY,
                double fillvalue, effective_earth_radius=None, blind_method="mask"):
    """
    Grid a PPI sweep onto a Cartesian plane.
    :param azimuth: azimuth coordinates for the first ``mat_ppi`` dimension.
    :param ranges: slant-range coordinates for the second ``mat_ppi`` dimension.
    :param elevation: sweep elevation angle.
    :param mat_ppi: PPI field data to interpolate.
    :param radar_height: radar altitude above mean sea level.
    :param GridX: target Cartesian X coordinates.
    :param GridY: target Cartesian Y coordinates.
    :param fillvalue: fill value used in ``mat_ppi``
    :return:
    """
    cdef int naz = azimuth.shape[0]
    cdef int nrange = ranges.shape[0]
    cdef int Nx = GridX.shape[0]
    cdef int Ny = GridY.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2] GridValue = np.empty((Nx, Ny), dtype=np.float64)
    cdef cnp.float64_t[:] azimuth_view = azimuth
    cdef cnp.float64_t[:] range_view = ranges
    cdef cnp.float64_t[:, :] mat_view = mat_ppi
    cdef cnp.float64_t[:, :] grid_x_view = GridX
    cdef cnp.float64_t[:, :] grid_y_view = GridY
    cdef cnp.float64_t[:, :] grid_value_view = GridValue
    cdef double az = 0.0
    cdef double r = 0.0
    cdef double z = 0.0
    cdef double radius = _resolve_effective_earth_radius(effective_earth_radius)
    cdef int blind_code = _resolve_blind_code(blind_method)
    cdef bint use_nearest_gate = blind_code in (1, 3, 4)
    cdef int ix, iy, iaz, ir
    if naz < 2 or nrange < 2:
        GridValue.fill(fillvalue)
        return GridValue
    for ix in range(Nx):
        for iy in range(Ny):
            _xye_to_antenna_impl(grid_x_view[ix, iy], grid_y_view[ix, iy], elevation, radar_height, radius, &az, &r, &z)
            grid_value_view[ix, iy] = _interp_sweep_sample(
                azimuth_view,
                range_view,
                mat_view,
                az,
                r,
                fillvalue,
                use_nearest_gate,
            )
    return GridValue

@cython.boundscheck(False)
@cython.wraparound(False)
def get_CAPPI_xy(vol_azimuth, vol_range, cnp.ndarray[cnp.float64_t, ndim=1] fix_elevation, vol_value,
                 double radar_height, cnp.ndarray[cnp.float64_t, ndim=2] GridX,
                 cnp.ndarray[cnp.float64_t, ndim=2] GridY, double level_height,  double fillvalue,
                 effective_earth_radius=None, blind_method="mask", double beam_width_deg=1.0):
    """
    Interpolate a CAPPI field from the radar volume.
    :param vol_azimuth: list of azimuth arrays, one per sweep.
    :param vol_range: list of range arrays, one per sweep.
    :param fix_elevation: fixed elevation angle for each sweep.
    :param vol_value: list of sweep data arrays.
    :param radar_height: radar altitude above mean sea level.
    :param GridX: target Cartesian X coordinates.
    :param GridY: target Cartesian Y coordinates.
    :param level_height: requested CAPPI height.
    :param fillvalue: fill value.
    :return:
    """
    cdef int Ne = fix_elevation.shape[0]
    cdef int Nx = GridX.shape[0]
    cdef int Ny = GridX.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=2] GridValue = np.full((Nx, Ny), fillvalue, dtype=np.float64)
    cdef cnp.float64_t[:] elevation_view = fix_elevation
    cdef cnp.float64_t[:, :] grid_x_view = GridX
    cdef cnp.float64_t[:, :] grid_y_view = GridY
    cdef cnp.float64_t[:, :] grid_value_view = GridValue
    cdef cnp.float64_t[:] azimuth_0
    cdef cnp.float64_t[:] azimuth_1
    cdef cnp.float64_t[:] range_0_view
    cdef cnp.float64_t[:] range_1_view
    cdef cnp.float64_t[:, :] value_0
    cdef cnp.float64_t[:, :] value_1
    cdef int ix, iy, ie, lower_index, upper_index
    cdef double az = 0.0
    cdef double el = 0.0
    cdef double r = 0.0
    cdef double radius = _resolve_effective_earth_radius(effective_earth_radius)
    cdef double IER0, IER1
    cdef double beam_half = 0.0
    cdef int blind_code = _resolve_blind_code(blind_method)
    cdef bint use_nearest_gate = blind_code in (1, 3, 4)
    cdef bint use_lowest_sweep = blind_code in (2, 3, 4)
    cdef bint use_highest_sweep = blind_code == 4
    if Ne == 0:
        return GridValue
    if beam_width_deg <= 0.0:
        beam_width_deg = 1.0
    beam_half = beam_width_deg * 0.5
    if beam_half < 0.1:
        beam_half = 0.1
    for ix in range(Nx):
        for iy in range(Ny):
            _cartesian_to_antenna_impl(grid_x_view[ix, iy], grid_y_view[ix, iy], level_height, radar_height, radius, &az, &r, &el)
            if Ne == 1:
                if fabs(el - elevation_view[0]) > beam_half:
                    continue
                azimuth_0 = vol_azimuth[0]
                range_0_view = vol_range[0]
                value_0 = vol_value[0]
                grid_value_view[ix, iy] = _interp_sweep_sample(
                    azimuth_0,
                    range_0_view,
                    value_0,
                    az,
                    r,
                    fillvalue,
                    use_nearest_gate,
                )
                continue
            if el < elevation_view[0]:
                if not use_lowest_sweep:
                    continue
                lower_index = 0
                upper_index = 0
            elif el > elevation_view[Ne - 1]:
                if not use_highest_sweep:
                    continue
                lower_index = Ne - 1
                upper_index = lower_index
            else:
                ie = _search_right_1d(elevation_view, el)
                if ie >= Ne:
                    ie = Ne - 1
                lower_index = ie - 1
                upper_index = ie

            azimuth_0 = vol_azimuth[lower_index]
            range_0_view = vol_range[lower_index]
            value_0 = vol_value[lower_index]
            IER0 = _interp_sweep_sample(
                azimuth_0,
                range_0_view,
                value_0,
                az,
                r,
                fillvalue,
                use_nearest_gate,
            )
            if lower_index == upper_index:
                grid_value_view[ix, iy] = IER0
                continue

            azimuth_1 = vol_azimuth[upper_index]
            range_1_view = vol_range[upper_index]
            value_1 = vol_value[upper_index]
            IER1 = _interp_sweep_sample(
                azimuth_1,
                range_1_view,
                value_1,
                az,
                r,
                fillvalue,
                use_nearest_gate,
            )
            grid_value_view[ix, iy] = _interp_linear(
                el,
                elevation_view[lower_index],
                elevation_view[upper_index],
                IER0,
                IER1,
                fillvalue,
            )
    return GridValue

def interp_azimuth(double az, double az_0, double az_1, double dat_0, double dat_1, double fillvalue):
    """Linearly interpolate between two azimuth or range samples."""
    return _interp_linear(az, az_0, az_1, dat_0, dat_1, fillvalue)

@cython.boundscheck(False)
@cython.wraparound(False)
def get_CAPPI_3d(
    vol_azimuth,
    vol_range,
    cnp.ndarray[cnp.float64_t, ndim=1] fix_elevation,
    vol_value,
    double radar_height,
    cnp.ndarray[cnp.float64_t, ndim=2] GridX,
    cnp.ndarray[cnp.float64_t, ndim=2] GridY,
    cnp.ndarray[cnp.float64_t, ndim=1] level_heights,
    double fillvalue,
    effective_earth_radius=None,
    blind_method="mask",
):
    """
    计算单部雷达在多个高度层上的 3D CAPPI 体.
    返回数组 shape 为 (nz, nx, ny)。
    """
    cdef int nlevel = level_heights.shape[0]
    cdef int Nx = GridX.shape[0]
    cdef int Ny = GridX.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=3] GridValue = np.full((nlevel, Nx, Ny), fillvalue, dtype=np.float64)
    cdef int ilevel
    for ilevel in range(nlevel):
        GridValue[ilevel, :, :] = get_CAPPI_xy(
            vol_azimuth,
            vol_range,
            fix_elevation,
            vol_value,
            radar_height,
            GridX,
            GridY,
            level_heights[ilevel],
            fillvalue,
            effective_earth_radius,
            blind_method,
        )
    return GridValue

@cython.boundscheck(False)
@cython.wraparound(False)
def get_mosaic_CAPPI_3d(
    vol_azimuth,
    vol_range,
    fix_elevation,
    vol_value,
    cnp.ndarray[cnp.float64_t, ndim=1] radar_x,
    cnp.ndarray[cnp.float64_t, ndim=1] radar_y,
    cnp.ndarray[cnp.float64_t, ndim=1] radar_height,
    cnp.ndarray[cnp.float64_t, ndim=2] GridX,
    cnp.ndarray[cnp.float64_t, ndim=2] GridY,
    cnp.ndarray[cnp.float64_t, ndim=1] level_heights,
    double fillvalue,
    effective_earth_radius=None,
):
    """
    多部雷达 3D 组网 CAPPI，使用最大值合成。
    输入雷达位置使用统一投影坐标系下的 x/y 偏移，单位为米。
    返回数组 shape 为 (nz, nx, ny)。
    """
    cdef int nradar = radar_height.shape[0]
    cdef int nlevel = level_heights.shape[0]
    cdef int Nx = GridX.shape[0]
    cdef int Ny = GridX.shape[1]
    cdef cnp.ndarray[cnp.float64_t, ndim=3] GridValue = np.full((nlevel, Nx, Ny), fillvalue, dtype=np.float64)
    cdef cnp.ndarray[cnp.float64_t, ndim=3] RadarValue
    cdef cnp.ndarray[cnp.float64_t, ndim=2] LocalX = np.empty((Nx, Ny), dtype=np.float64)
    cdef cnp.ndarray[cnp.float64_t, ndim=2] LocalY = np.empty((Nx, Ny), dtype=np.float64)
    cdef cnp.float64_t[:, :] grid_x_view = GridX
    cdef cnp.float64_t[:, :] grid_y_view = GridY
    cdef cnp.float64_t[:, :] local_x_view = LocalX
    cdef cnp.float64_t[:, :] local_y_view = LocalY
    cdef cnp.float64_t[:, :, :] grid_value_view = GridValue
    cdef cnp.float64_t[:, :, :] radar_value_view
    cdef object radar_radius = effective_earth_radius
    cdef int iradar, ilevel, ix, iy
    if not (len(vol_azimuth) == len(vol_range) == len(fix_elevation) == len(vol_value) == nradar):
        raise ValueError("Radar volume lists and radar position arrays must have the same length.")
    if radar_radius is None:
        radar_radius = [None] * nradar
    elif not isinstance(radar_radius, (list, tuple)):
        radar_radius = [radar_radius] * nradar
    elif len(radar_radius) != nradar:
        raise ValueError("effective_earth_radius must be None, a scalar, or a sequence matching the radar count.")
    for iradar in range(nradar):
        for ix in range(Nx):
            for iy in range(Ny):
                local_x_view[ix, iy] = grid_x_view[ix, iy] - radar_x[iradar]
                local_y_view[ix, iy] = grid_y_view[ix, iy] - radar_y[iradar]
        RadarValue = get_CAPPI_3d(
            vol_azimuth[iradar],
            vol_range[iradar],
            fix_elevation[iradar],
            vol_value[iradar],
            radar_height[iradar],
            LocalX,
            LocalY,
            level_heights,
            fillvalue,
            radar_radius[iradar],
        )
        radar_value_view = RadarValue
        for ilevel in range(nlevel):
            for ix in range(Nx):
                for iy in range(Ny):
                    if radar_value_view[ilevel, ix, iy] == fillvalue:
                        continue
                    if grid_value_view[ilevel, ix, iy] == fillvalue or radar_value_view[ilevel, ix, iy] > grid_value_view[ilevel, ix, iy]:
                        grid_value_view[ilevel, ix, iy] = radar_value_view[ilevel, ix, iy]
    return GridValue

@cython.boundscheck(False)  # Deactivate bounds checking
def get_CR_xy(vol_azimuth, vol_range, cnp.ndarray[cnp.float64_t, ndim=1] fix_elevation, vol_value, double radar_height,
              cnp.ndarray[cnp.float64_t, ndim=2] GridX, cnp.ndarray[cnp.float64_t, ndim=2] GridY, double fillvalue,
              effective_earth_radius=None):
    """
    Compute composite reflectivity from the radar volume.
    :param vol_azimuth: list of azimuth arrays, one per sweep.
    :param vol_range: list of range arrays, one per sweep.
    :param fix_elevation: fixed elevation angle for each sweep.
    :param vol_value: list of sweep data arrays.
    :param radar_height: radar altitude above mean sea level.
    :param GridX: target Cartesian X coordinates.
    :param GridY: target Cartesian Y coordinates.
    :param fillvalue: fill value.
    :return:
    """
    cdef int Ne = fix_elevation.shape[0]
    cdef int Nx = GridX.shape[0]
    cdef int Ny = GridY.shape[1]
    cdef cnp.ndarray GridValue = np.zeros([Ne, Nx, Ny], dtype=np.float64)
    cdef int ie
    for ie in range(Ne):
        GridValue[ie,:,:] = ppi_to_grid(vol_azimuth[ie], vol_range[ie], fix_elevation[ie],
                                       vol_value[ie], radar_height, GridX, GridY, fillvalue, effective_earth_radius)

    GridValue = np.where(GridValue == fillvalue, np.nan, GridValue)
    valid_mask = ~np.isnan(GridValue)
    GridValue = np.max(np.where(valid_mask, GridValue, -np.inf), axis=0)
    GridValue = np.where(valid_mask.any(axis=0), GridValue, fillvalue)
    return GridValue
