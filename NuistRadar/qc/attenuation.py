# -*- coding: utf-8 -*-
"""
edited in 2019/11/25
"""

import numpy as np

def correct_attenuation_HB(Ref, a=1.67e-4, b=0.7, gate_length=0.075, thrs=59):
    """
    Hitschfeld1954, fillvalue set np.nan
    :param Ref: radar ref //dbz [naz, nbins]
    :param a: coefficients a
    :param b: coefficients b
    :param gate_length: bin length //m
    :param thrs: thrs //dbz
    :return:
    """
    if hasattr(Ref, "range"):
        gate_length = ((Ref.range[1] - Ref.range[0])/1000.).values ##判断ref有库的属性,或者库长

    Ref = np.where(np.isnan(Ref), -999, Ref)
    pia = np.zeros(Ref.shape)
    ksum = 0.
    # multidimensional version
    # assumes that iteration is only along the last dimension
    # (i.e. range gates) all other dimensions are calculated simultaneously
    # to gain some speed
    for gate in range(Ref.shape[-1] - 1):
        # calculate k in dB/km from k-Z relation
        # c.f. Krämer2008(p. 147)
        k = a * (10.0 ** ((Ref[:, gate] + ksum) / 10.0)) ** b * 2.0 * gate_length
        # k = 10**(log10(a)+0.1*bin*b)
        # dBkn = 10*math.log10(a) + (bin+ksum)*b + 10*math.log10(2*gate_length)
        ksum += k
        pia[:, gate + 1] = ksum
        # stop-criterion, if corrected reflectivity is larger than 59 dBZ
        overflow = (Ref[..., gate + 1] + ksum) > thrs
        if np.any(overflow):
            pia[:, gate + 1][overflow] = np.nan
    pia = np.where(np.isnan(pia), 0, pia)
    Zc = np.where(Ref==-999, np.nan, pia + Ref)
    return Zc, pia

def correct_attenuation(ref, wavelength="C", rscale=0.075):

    '''Scientific paper: Ośródka, K., Szturc, J., and Jurczyk, A., 2012.
    Chain of data quality algorithms for 3-D single-polarization radar
    reflectivity (RADVOL-QC system). Meteorol. Appl. (Early View).
    rscale : radar gates length(km)
    ref :reflectivity factor (dBZ) fillvalue=np.nan'''
    if wavelength == "C":
        a = 0.0044  # Coefficient a in attenuation formula
        b = 1.17  # Coefficient b in attenuation formula
    elif wavelength == "X":
        a = 0.0148
        b = 1.31
    elif wavelength == "S":
        a = 0.0006
        b = 1.00
    else:
        raise ValueError('Unrecognized type')

    if hasattr(ref, "range"):
        rscale = ((ref.range[1] - ref.range[0])/1000.).values ##判断ref有库的属性

    Ra = 200  # Coefficient a in Z-R formula
    Rb = 1.6  # Coefficient b in Z-R formula
    ATT_refl = 4.0  # 低于此值不计算衰减
    a_att = 1.0
    b_att = 5.0  # Quality characterization
    nray, ngate = ref.shape
    Z = np.where(ref > ATT_refl, 10 ** (ref / 10.), 0)  ##lower than refl set 0
    A = rscale * a * (1. / Ra * Z) ** (b / Rb)
    PIA = np.zeros_like(A)
    Asum = np.zeros((nray,))
    for igate in range(1, ngate):
        Asum = Asum + A[:, igate - 1]
        PIA[:, igate] = Asum
    Z_corr = np.where(np.isnan(ref), ref, ref + PIA)
    QI = (b_att - PIA) / (b_att - a_att)
    QI = np.where(PIA < a_att, 1, QI)
    QI = np.where(PIA > b_att, 0, QI)
    return Z_corr, QI

def pia_from_kdp(kdp, dr, gamma=0.08):
    """Retrieving path integrated attenuation from specific differential \
    phase (Kdp).
    The default value of gamma is based on :cite:`Carey2000`.
    Parameters
    ----------
    kdp : :class:`numpy:numpy.ndarray`
       array specific differential phase
       Range dimension must be the last dimension.
    dr : gate length (km)
    gamma : float
       linear coefficient (default value: 0.08) in the relation between
       Kdp phase and specific attenuation (alpha)
    Returns
    -------
    output : :class:`numpy:numpy.ndarray`
        array of same shape as kdp containing the path integrated attenuation
    """
    alpha = gamma * kdp
    return 2 * np.cumsum(alpha, axis=-1) * dr

