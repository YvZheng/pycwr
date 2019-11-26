import numpy as np
from scipy import spatial
from ..configure.config import cfg

def get_weight(dist, r, method="barnes"):
    """
    barnes权重函数,给定距离dist和影响半径r,返回权重值
    :param dist: 数据点距离插值点的距离
    :param r: 有效影响半径
    :param method 插值方法
    :return: 该点在插值上的权重
    """
    if method == "barnes":
        weight = np.exp(-4*dist**2/r**2)
    elif method == "cressman":
        weight = (r ** 2 - dist ** 2) / (dist ** 2 + r ** 2)
    else:
        raise Exception("Unidentified method!, must be cressman, barnes")
    return weight

def _get_interp_around_point(point_old, point_new, around_r):
    """
    需要point_new周围r范围以内的所有点
    :param point_old: 原始的点 np.c_[x, y, z] or np.c_[x, y]
    :param point_new: 需要插值到的点 np.c_[x', y', z'] or np.c_[x', y']
    :param around_r: 范围 unit:m
    :return:
    """
    kdtree = spatial.cKDTree(point_old)
    index_nearest = kdtree.query_ball_point(point_new, around_r) ##找到距离around_r以内所有点在point_old的index
    dist = [np.sqrt(np.sum(np.square(point_old[i,...] - itarget), axis=1)) for i, itarget \
            in zip(index_nearest, point_new)]   ##对应每个index_nearest的距离
    return index_nearest, dist

def radar_interp2d(points, values, xi, around_r,  influence_radius=None, method="barnes", fill_value=np.nan):
    """
    Interpolate unstructured D-dimensional data.
    Parameters
    ----------
    points : ndarray of floats, shape (n, D)
        Data point coordinates. Can either be an array of
        shape (n, D), or a tuple of `ndim` arrays.
    values : ndarray of float or complex, shape (n,)
        Data values.
    xi : 2-D ndarray of float or tuple of 1-D array, shape (M, D)
        Points at which to interpolate data.
    around_r: 只取周围around_r点插值
    influence_radius: 插值函数的影响半径
    method : {'barnes', 'cressman'}
        Method of interpolation. One of 'barnes', 'cressman'
    fill_value : float, optional
        Value used to fill in for requested points outside of the
        convex hull of the input points.  If not provided, then the
        default is ``nan``. This option has no effect for the
        'nearest' method.
    Returns
    -------
    ndarray
        Array of interpolated values.
    """
    if influence_radius is None:
        influence_radius = around_r
    grid_shape = xi[0].shape
    target = np.column_stack([xi_grid.ravel() for xi_grid in xi])
    index, distance = _get_interp_around_point(points, target, around_r)
    nrows, _ = target.shape
    grid_vals = np.empty(nrows)
    for i in range(nrows):
        if index[i]:
            weight = get_weight(distance[i], influence_radius, method=method)
            grid_vals[i] = np.dot(values[index[i]], weight)/np.sum(weight)
        else:
            grid_vals[i] = fill_value
    return grid_vals.reshape(grid_shape)

def _get_interp_around_point_var(point_old, point_new, bandwidth=1):
    """
    最小影响半径min_roi设置为200m, 影响半径随着雷达距离变化
    :param point_old:
    :param point_new:
    :param bandwidth: 波束宽度 degree
    :return:
    """
    min_roi = cfg.interp.mroi
    kdtree = spatial.cKDTree(point_old)
    index_nearest = []
    nrows = point_new.shape[0]
    roi = np.empty(nrows)
    for i, itarget in enumerate(point_new):
        roi[i] = min_roi + ((itarget[0] / 1000.) ** 2 + (itarget[1] / 1000.) ** 2) ** 0.5 * bandwidth * cfg.interp.coeff
        index_nearest.append(kdtree.query_ball_point(itarget, roi[i]))  ###idx放的是索引
    distance = [np.sqrt(np.sum(np.square(point_old[i, :] - j), axis=1)) for i, j \
            in zip(index_nearest, point_new)]
    return index_nearest, distance, roi

def radar_interp2d_var(points, values, xi, bandwidth=1, method="barnes", fill_value=np.nan):
    """
    影响半径随着距离雷达中心的距离发生变化
    :param points:要插值出去的点
    :param values:
    :param xi:
    :param bandwidth: 波束宽度
    :param method:
    :param fill_value:
    :return:
    """

    grid_shape = xi[0].shape
    target = np.column_stack([xi_grid.ravel() for xi_grid in xi])
    index, distance, roi = _get_interp_around_point_var(points, target, bandwidth=bandwidth)
    nrows, _ = target.shape
    grid_vals = np.empty(nrows)
    for i in range(nrows):
        if index[i]:
            weight = get_weight(distance[i], roi[i], method=method)
            grid_vals[i] = np.dot(values[index[i]], weight) / np.sum(weight)
        else:
            grid_vals[i] = fill_value
    return grid_vals.reshape(grid_shape)


