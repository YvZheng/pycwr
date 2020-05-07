# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 18:47:02 2020
part of code is from https://github.com/CSU-Radarmet/CSU_RadarTools
双偏振雷达水凝物分类算法 (Hydrometeor identification algorithms)
@author: zy
"""
from ..configure.location_config import mbf_path
import xarray as xr
import numpy as np

DEFAULT_WEIGHTS = {'dBZ': 1.5, 'ZDR': 0.8, 'KDP': 1.0, 'CC': 0.8, 'LDR': 0.5, 'T': 0.4}
beta_params = xr.open_dataset(mbf_path)

def hid_beta_function(x, m, a, b):
    """
    :param x: input data(Zh, Zdr, Kdp...)
    :param m: center
    :param a: width
    :param b: slope
    :return:
    """
    if None in x:
        return 0
    else:
        return 1 / (1 + ((x - m) / a) ** (2 * b))

def fhc_HCL(dBZ=None, ZDR=None, KDP=None, CC=None, LDR=None, T=None, method="hybrid",
            band="C", weights=DEFAULT_WEIGHTS):
    """
    Does FHC for warm-season precip. using beta function for Fuzzy Logic
    HCL types:           Species #:
    -------------------------------
    Drizzle                  1
    Rain                     2
    Ice Crystals             3
    Dry Aggregates Snow      4
    Wet Snow                 5
    Vertical Ice             6
    Low-Density Graupel      7
    High-Density Graupel     8
    Hail                     9
    Big Drops                10
    -------------------------------
    cite{Dolan, Brenda , et al. "A Robust C-Band Hydrometeor Identification Algorithm and Application to
     a Long-Term Polarimetric Radar Dataset." Journal of Applied Meteorology and Climatology 52.9(2013):2162-2186.}
    :param weights:
    :param dBZ: Input reflectivity scalar/array
    :param ZDR: Input differential reflectivity scalar/array
    :param KDP: Input specific differential phase scalar/array
    :param CC: Input cross correlation ratio scalar/array
    :param LDR: Input linear depolarization ratio scalar/array
    :param T: Input temperature scalar/array
    :param method: Currently support 'hybrid' or 'linear' methods; hybrid preferred
    :param band: 'X', 'C', or 'S'
    :return: hydrometeor species [1-10]
    """
    if LDR is None and ZDR is None and KDP is None and CC is None:
        print("No Polarimetric variable input!")
        return
    if dBZ is None:
        print("No reflectivity variable input!")
        return
    if LDR is None:
        weights['LDR'] = 0
    if ZDR is None:
        weights['ZDR'] = 0
    if KDP is None:
        weights['KDP'] = 0
    if T is None:
        weights['T'] = 0
    if CC is None:
        weights['CC'] = 0
    Beta_ZDR = hid_beta_function(np.stack([ZDR] * 10, axis=-1),
                                 beta_params.MBF.sel(Band=band, Feature="ZDR", Param="m").values,
                                 beta_params.MBF.sel(Band=band, Feature="ZDR", Param="a").values,
                                 beta_params.MBF.sel(Band=band, Feature="ZDR", Param="b").values)
    Beta_dBZ = hid_beta_function(np.stack([dBZ] * 10, axis=-1),
                                 beta_params.MBF.sel(Band=band, Feature="dBZ", Param="m").values,
                                 beta_params.MBF.sel(Band=band, Feature="dBZ", Param="a").values,
                                 beta_params.MBF.sel(Band=band, Feature="dBZ", Param="b").values)
    Beta_KDP = hid_beta_function(np.stack([KDP] * 10, axis=-1),
                                 beta_params.MBF.sel(Band=band, Feature="KDP", Param="m").values,
                                 beta_params.MBF.sel(Band=band, Feature="KDP", Param="a").values,
                                 beta_params.MBF.sel(Band=band, Feature="KDP", Param="b").values)
    Beta_CC = hid_beta_function(np.stack([CC] * 10, axis=-1),
                                beta_params.MBF.sel(Band=band, Feature="CC", Param="m").values,
                                beta_params.MBF.sel(Band=band, Feature="CC", Param="a").values,
                                beta_params.MBF.sel(Band=band, Feature="CC", Param="b").values)
    Beta_LDR = hid_beta_function(np.stack([LDR] * 10, axis=-1),
                                 beta_params.MBF.sel(Band=band, Feature="LDR", Param="m").values,
                                 beta_params.MBF.sel(Band=band, Feature="LDR", Param="a").values,
                                 beta_params.MBF.sel(Band=band, Feature="LDR", Param="b").values)
    if T is None:
        Beta_T = 1
    else:
        Beta_T = hid_beta_function(np.stack([T] * 10, axis=-1),
                                   beta_params.MBF.sel(Band=band, Feature="T", Param="m").values,
                                   beta_params.MBF.sel(Band=band, Feature="T", Param="a").values,
                                   beta_params.MBF.sel(Band=band, Feature="T", Param="b").values)
    if method == "hybrid":
        HCL = (Beta_LDR*weights['LDR'] + Beta_KDP*weights['KDP'] + Beta_ZDR*weights['ZDR'] + Beta_CC*weights['CC'])/\
              (weights['LDR'] + weights['KDP'] + weights['ZDR'] + weights['CC'])*Beta_T*Beta_dBZ
    elif method=="linear":
        HCL = (Beta_LDR*weights['LDR'] + Beta_KDP*weights['KDP'] + Beta_ZDR*weights['ZDR'] + Beta_CC*weights['CC'] +\
               Beta_T * weights['T'] + Beta_dBZ * weights['dBZ'])/(weights['LDR'] + weights['KDP'] + weights['ZDR'] +\
                                                                   weights['CC'] + weights['T'] + weights['dBZ'])
    else:
        print("No weighting method defined, use hybrid or linear")
        return
    return np.where(np.any(np.isnan(HCL), axis=-1), np.nan, np.argmax(HCL, axis=-1) + 1)

