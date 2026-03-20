import numpy as np

from .transforms import (
    _azimuth,
    antenna_to_cartesian_cwr,
    cartesian_xy_elevation_to_range_z,
    cartesian_xyz_to_antenna,
)


def _search_right_1d(values, target):
    return int(np.searchsorted(values, target, side="right"))


def _merge_fill_max(base, candidate, fillvalue):
    out = np.array(base, copy=True)
    valid = candidate != fillvalue
    if not np.any(valid):
        return out
    empty = out == fillvalue
    write_mask = valid & empty
    out[write_mask] = candidate[write_mask]
    merge_mask = valid & (~empty)
    out[merge_mask] = np.maximum(out[merge_mask], candidate[merge_mask])
    return out


def _resolve_blind_flags(blind_method):
    blind_method = "mask" if blind_method is None else str(blind_method).lower()
    if blind_method == "mask":
        return False, False, False
    if blind_method == "nearest_gate":
        return True, False, False
    if blind_method == "lowest_sweep":
        return False, True, False
    if blind_method == "hybrid":
        return True, True, False
    if blind_method == "nearest_valid":
        return True, True, True
    raise ValueError(
        "blind_method must be one of 'mask', 'nearest_gate', 'lowest_sweep', "
        "'hybrid', or 'nearest_valid'."
    )


def _interp_sweep_sample(azimuth, ranges, values, target_az, target_range, fillvalue, use_nearest_gate):
    if azimuth.size < 2 or ranges.size < 2:
        return fillvalue
    if target_range > ranges[-1]:
        return fillvalue
    if target_range < ranges[0]:
        if not use_nearest_gate:
            return fillvalue
        target_range = ranges[0]

    iaz = _search_right_1d(azimuth, target_az)
    if iaz >= azimuth.size:
        iaz = 0
        az_value = target_az - 360.0
    else:
        az_value = target_az
    if iaz == 0:
        az_last = azimuth[-1] - 360.0
    else:
        az_last = azimuth[iaz - 1]

    ir = _search_right_1d(ranges, target_range)
    if ir <= 0:
        ir = 1
    elif ir >= ranges.size:
        ir = ranges.size - 1

    er0 = interp_azimuth(
        az_value,
        az_last,
        azimuth[iaz],
        values[iaz - 1, ir - 1],
        values[iaz, ir - 1],
        fillvalue,
    )
    er1 = interp_azimuth(
        az_value,
        az_last,
        azimuth[iaz],
        values[iaz - 1, ir],
        values[iaz, ir],
        fillvalue,
    )
    return interp_azimuth(
        target_range,
        ranges[ir - 1],
        ranges[ir],
        er0,
        er1,
        fillvalue,
    )


antenna_to_cartesian = antenna_to_cartesian_cwr
xye_to_antenna = cartesian_xy_elevation_to_range_z
cartesian_to_antenna = cartesian_xyz_to_antenna
xy_to_azimuth = _azimuth


def interp_ppi(az, r, az_0, az_1, r_0, r_1, mat_00, mat_01, mat_10, mat_11, fillvalue=-999.0):
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
    az_span = az_1 - az_0
    r_span = r_1 - r_0
    if az_span == 0.0 or r_span == 0.0:
        return fillvalue
    if ((mat_00 != fillvalue) and (mat_01 != fillvalue)) and ((mat_10 != fillvalue) and (mat_11 != fillvalue)):
        interped = (
            mat_00 * (az_1 - az) * (r_1 - r)
            + mat_10 * (az - az_0) * (r_1 - r)
            + mat_01 * (az_1 - az) * (r - r_0)
            + mat_11 * (az - az_0) * (r - r_0)
        ) / r_span / az_span
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


def interp_azimuth(az, az_0, az_1, dat_0, dat_1, fillvalue=-999.0):
    """
    在两个方位角或者距离之间进行插值
    """
    if az_1 == az_0:
        return fillvalue
    if (dat_0 == fillvalue) and (dat_1 == fillvalue):
        return fillvalue
    if dat_0 == fillvalue:
        return dat_1
    if dat_1 == fillvalue:
        return dat_0
    return ((az_1 - az) * dat_0 + (az - az_0) * dat_1) / (az_1 - az_0)


def ppi_to_grid(azimuth, ranges, elevation, mat_ppi, radar_height, GridX, GridY, fillvalue=-999.0,
                effective_earth_radius=None, blind_method="mask"):
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
    GridX = np.asarray(GridX, dtype=np.float64)
    GridY = np.asarray(GridY, dtype=np.float64)
    azimuth = np.asarray(azimuth, dtype=np.float64)
    ranges = np.asarray(ranges, dtype=np.float64)
    mat_ppi = np.asarray(mat_ppi, dtype=np.float64)

    GridValue = np.full_like(GridX, fillvalue, dtype=np.float64)
    if azimuth.size < 2 or ranges.size < 2:
        return GridValue
    use_nearest_gate, _, _ = _resolve_blind_flags(blind_method)

    Grid_az, Grid_range, _ = cartesian_xy_elevation_to_range_z(
        GridX,
        GridY,
        elevation,
        radar_height,
        effective_earth_radius=effective_earth_radius,
    )
    Nx, Ny = GridX.shape
    for ix in range(Nx):
        for iy in range(Ny):
            GridValue[ix, iy] = _interp_sweep_sample(
                azimuth,
                ranges,
                mat_ppi,
                Grid_az[ix, iy],
                Grid_range[ix, iy],
                fillvalue,
                use_nearest_gate,
            )
    return GridValue


def get_CR_xy(vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, GridX, GridY, fillvalue=-999.0,
              effective_earth_radius=None):
    """
    计算组合反射率，利用雷达体扫的数据
    """
    fix_elevation = np.asarray(fix_elevation, dtype=np.float64)
    Nx, Ny = GridX.shape
    GridValue = np.zeros([fix_elevation.shape[0], Nx, Ny], dtype=np.float64)
    for ie in range(fix_elevation.shape[0]):
        GridValue[ie] = ppi_to_grid(
            vol_azimuth[ie],
            vol_range[ie],
            fix_elevation[ie],
            vol_value[ie],
            radar_height,
            GridX,
            GridY,
            fillvalue,
            effective_earth_radius=effective_earth_radius,
        )
    GridValue = np.where(GridValue == fillvalue, np.nan, GridValue)
    valid_mask = ~np.isnan(GridValue)
    reduced = np.max(np.where(valid_mask, GridValue, -np.inf), axis=0)
    return np.where(valid_mask.any(axis=0), reduced, fillvalue)


def get_CAPPI_xy(vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, GridX, GridY, level_height,
                 fillvalue=-999.0, effective_earth_radius=None, blind_method="mask", beam_width_deg=1.0):
    """
    由雷达体扫数据，插值CAPPI
    """
    fix_elevation = np.asarray(fix_elevation, dtype=np.float64)
    GridX = np.asarray(GridX, dtype=np.float64)
    GridY = np.asarray(GridY, dtype=np.float64)
    GridValue = np.full_like(GridX, fillvalue, dtype=np.float64)
    if fix_elevation.size == 0:
        return GridValue
    use_nearest_gate, use_lowest_sweep, use_highest_sweep = _resolve_blind_flags(blind_method)

    Grid_az, Grid_range, Grid_el = cartesian_xyz_to_antenna(
        GridX,
        GridY,
        level_height,
        radar_height,
        effective_earth_radius=effective_earth_radius,
    )
    Nx, Ny = GridX.shape
    for ix in range(Nx):
        for iy in range(Ny):
            current_el = Grid_el[ix, iy]
            if fix_elevation.size == 1:
                if np.abs(current_el - fix_elevation[0]) > max(float(beam_width_deg) * 0.5, 0.1):
                    continue
                GridValue[ix, iy] = _interp_sweep_sample(
                    np.asarray(vol_azimuth[0], dtype=np.float64),
                    np.asarray(vol_range[0], dtype=np.float64),
                    np.asarray(vol_value[0], dtype=np.float64),
                    Grid_az[ix, iy],
                    Grid_range[ix, iy],
                    fillvalue,
                    use_nearest_gate,
                )
                continue
            if current_el < fix_elevation[0]:
                if not use_lowest_sweep:
                    continue
                lower_index = 0
                upper_index = 0
            elif current_el > fix_elevation[-1]:
                if not use_highest_sweep:
                    continue
                lower_index = fix_elevation.size - 1
                upper_index = lower_index
            else:
                upper_index = _search_right_1d(fix_elevation, current_el)
                if upper_index >= fix_elevation.size:
                    upper_index = fix_elevation.size - 1
                lower_index = upper_index - 1

            current_az = Grid_az[ix, iy]
            current_range = Grid_range[ix, iy]
            azimuth_0 = np.asarray(vol_azimuth[lower_index], dtype=np.float64)
            range_0_view = np.asarray(vol_range[lower_index], dtype=np.float64)
            value_0 = np.asarray(vol_value[lower_index], dtype=np.float64)
            IER0 = _interp_sweep_sample(
                azimuth_0,
                range_0_view,
                value_0,
                current_az,
                current_range,
                fillvalue,
                use_nearest_gate,
            )
            if lower_index == upper_index:
                GridValue[ix, iy] = IER0
                continue

            azimuth_1 = np.asarray(vol_azimuth[upper_index], dtype=np.float64)
            range_1_view = np.asarray(vol_range[upper_index], dtype=np.float64)
            value_1 = np.asarray(vol_value[upper_index], dtype=np.float64)
            IER1 = _interp_sweep_sample(
                azimuth_1,
                range_1_view,
                value_1,
                current_az,
                current_range,
                fillvalue,
                use_nearest_gate,
            )
            GridValue[ix, iy] = interp_azimuth(
                current_el,
                fix_elevation[lower_index],
                fix_elevation[upper_index],
                IER0,
                IER1,
                fillvalue,
            )

    return GridValue


def get_CAPPI_3d(vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, GridX, GridY, level_heights,
                 fillvalue=-999.0, effective_earth_radius=None, blind_method="mask"):
    """
    计算单部雷达在多个高度层上的 3D CAPPI 体.
    返回 shape 为 (nz, nx, ny) 的数组。
    """
    level_heights = np.asarray(level_heights, dtype=np.float64)
    GridX = np.asarray(GridX, dtype=np.float64)
    GridY = np.asarray(GridY, dtype=np.float64)
    GridValue = np.full((level_heights.size, GridX.shape[0], GridX.shape[1]), fillvalue, dtype=np.float64)
    for iz, level_height in enumerate(level_heights):
        GridValue[iz] = get_CAPPI_xy(
            vol_azimuth,
            vol_range,
            fix_elevation,
            vol_value,
            radar_height,
            GridX,
            GridY,
            float(level_height),
            fillvalue,
            effective_earth_radius=effective_earth_radius,
            blind_method=blind_method,
        )
    return GridValue


def get_mosaic_CAPPI_3d(
    vol_azimuth,
    vol_range,
    fix_elevation,
    vol_value,
    radar_x,
    radar_y,
    radar_height,
    GridX,
    GridY,
    level_heights,
    fillvalue=-999.0,
    effective_earth_radius=None,
):
    """
    多部雷达 3D 组网 CAPPI，使用最大值合成。
    输入雷达位置使用统一投影坐标系下的 x/y 偏移，单位为米。
    返回 shape 为 (nz, nx, ny) 的数组。
    """
    radar_x = np.asarray(radar_x, dtype=np.float64)
    radar_y = np.asarray(radar_y, dtype=np.float64)
    radar_height = np.asarray(radar_height, dtype=np.float64)
    GridX = np.asarray(GridX, dtype=np.float64)
    GridY = np.asarray(GridY, dtype=np.float64)
    level_heights = np.asarray(level_heights, dtype=np.float64)

    nradar = radar_height.size
    if not (len(vol_azimuth) == len(vol_range) == len(fix_elevation) == len(vol_value) == nradar):
        raise ValueError("Radar volume lists and radar position arrays must have the same length.")
    if effective_earth_radius is None or np.isscalar(effective_earth_radius):
        effective_earth_radius = [effective_earth_radius] * nradar
    elif len(effective_earth_radius) != nradar:
        raise ValueError("effective_earth_radius must be None, a scalar, or a sequence matching the radar count.")

    GridValue = np.full((level_heights.size, GridX.shape[0], GridX.shape[1]), fillvalue, dtype=np.float64)
    for iradar in range(nradar):
        radar_grid = get_CAPPI_3d(
            vol_azimuth[iradar],
            vol_range[iradar],
            fix_elevation[iradar],
            vol_value[iradar],
            radar_height[iradar],
            GridX - radar_x[iradar],
            GridY - radar_y[iradar],
            level_heights,
            fillvalue,
            effective_earth_radius=effective_earth_radius[iradar],
        )
        GridValue = _merge_fill_max(GridValue, radar_grid, fillvalue)
    return GridValue
