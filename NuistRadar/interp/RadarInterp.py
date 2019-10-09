import numpy as np
from scipy import spatial


def weight_barnes(dist, r):
    """
    barnes权重函数,给定距离dist和影响半径r,返回权重值
    :param dist: 数据点距离插值点的距离
    :param r:  Barnes的有效影响半径
    :return: 该点在插值上的权重
    """
    weight = np.exp(-4*dist**2/r**2)
    return weight

def weight_cressman(dist, r):
    """
    cressman权重函数
    :param dist: 数据点距离插值点的距离
    :param r: cressman的有效影响半径
    :return: 对应的权重值
    """
    weight = (r**2-dist**2)/(dist**2+r**2)
    return weight