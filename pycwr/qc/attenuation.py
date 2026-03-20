# -*- coding: utf-8 -*-
"""
Attenuation correction helpers.
"""

import numpy as np


KDP_ATTENUATION_COEFFICIENTS = {
    "S": {"gamma": 0.02, "beta": 0.002},
    "C": {"gamma": 0.08, "beta": 0.02},
    "X": {"gamma": 0.28, "beta": 0.04},
}


def _as_float_array(data):
    return np.asanyarray(data, dtype=float)


def _resolve_gate_length_km(data, gate_length):
    if hasattr(data, "range"):
        range_values = np.asanyarray(data.range.values, dtype=float)
        if range_values.size > 1:
            return float(np.median(np.diff(range_values)) / 1000.0)
    return float(gate_length)


def resolve_kdp_coefficients(band="C", gamma=None, beta=None):
    """Return attenuation coefficients for a radar wavelength band."""
    band = str(band).upper()
    if band not in KDP_ATTENUATION_COEFFICIENTS:
        raise ValueError("band must be one of 'S', 'C' or 'X'")
    coefficients = KDP_ATTENUATION_COEFFICIENTS[band]
    if gamma is None:
        gamma = coefficients["gamma"]
    if beta is None:
        beta = coefficients["beta"]
    return float(gamma), float(beta)


def _prepare_mask(data, mask=None):
    valid = np.isfinite(data)
    if mask is None:
        return valid
    return valid & np.asarray(mask, dtype=bool)


def correct_attenuation_HB(Ref, a=1.67e-4, b=0.7, gate_length=0.075, thrs=59):
    """
    Hitschfeld-Bordan attenuation correction.

    Parameters
    ----------
    Ref : array-like
        Reflectivity in dBZ. Range dimension must be the last dimension.
    a, b : float
        Coefficients in the k-Z relation.
    gate_length : float
        Gate length in km. If ``Ref`` is an xarray field exposing ``range``,
        the gate spacing is derived from it automatically.
    thrs : float
        Upper reflectivity threshold. Once exceeded, the current valid segment
        is terminated and the remaining gates in that segment are marked NaN.
    """
    reflectivity = _as_float_array(Ref)
    dr = _resolve_gate_length_km(Ref, gate_length)
    flat_ref = reflectivity.reshape(-1, reflectivity.shape[-1])
    flat_zc = np.full_like(flat_ref, np.nan)
    flat_pia = np.full_like(flat_ref, np.nan)

    for iray, row in enumerate(flat_ref):
        running = 0.0
        overflowed = False
        for igate, value in enumerate(row):
            if not np.isfinite(value):
                running = 0.0
                overflowed = False
                continue
            if overflowed:
                continue
            corrected = value + running
            if corrected > thrs:
                overflowed = True
                continue
            flat_pia[iray, igate] = running
            flat_zc[iray, igate] = corrected
            specific_attenuation = a * (10.0 ** (corrected / 10.0)) ** b
            running += 2.0 * dr * specific_attenuation

    return flat_zc.reshape(reflectivity.shape), flat_pia.reshape(reflectivity.shape)


def correct_attenuation(ref, wavelength="C", rscale=0.075):
    """
    Legacy single-polarization attenuation correction from the historical API.

    This function is retained for backward compatibility. New dual-polarization
    workflows should prefer :func:`correct_attenuation_kdp` or
    :func:`pycwr.qc.run_dualpol_qc`.
    """
    if wavelength == "C":
        a = 0.0044
        b = 1.17
    elif wavelength == "X":
        a = 0.0148
        b = 1.31
    elif wavelength == "S":
        a = 0.0006
        b = 1.00
    else:
        raise ValueError("Unrecognized type")

    if hasattr(ref, "range"):
        rscale = float((ref.range[1] - ref.range[0]) / 1000.0)

    Ra = 200
    Rb = 1.6
    ATT_refl = 4.0
    a_att = 1.0
    b_att = 5.0
    reflectivity = _as_float_array(ref)
    nray, ngate = reflectivity.shape
    Z = np.where(reflectivity > ATT_refl, 10 ** (reflectivity / 10.0), 0.0)
    A = rscale * a * (1.0 / Ra * Z) ** (b / Rb)
    PIA = np.zeros_like(A)
    Asum = np.zeros((nray,))
    for igate in range(1, ngate):
        Asum = Asum + A[:, igate - 1]
        PIA[:, igate] = Asum
    Z_corr = np.where(np.isnan(reflectivity), reflectivity, reflectivity + PIA)
    QI = (b_att - PIA) / (b_att - a_att)
    QI = np.where(PIA < a_att, 1, QI)
    QI = np.where(PIA > b_att, 0, QI)
    return Z_corr, QI


def pia_from_kdp(kdp, dr, gamma=0.08, mask=None, clip_negative=True):
    """
    Retrieve path integrated attenuation from specific differential phase.

    Parameters
    ----------
    kdp : array-like
        Specific differential phase (deg/km). Range dimension must be the last
        dimension.
    dr : float
        Gate length in km.
    gamma : float
        Linear coefficient between KDP and specific attenuation.
    mask : array-like, optional
        Boolean mask indicating gates to include in the integration. Invalid
        gates terminate the current contiguous segment.
    clip_negative : bool
        If True, negative KDP is clipped to zero before integration.
    """
    kdp = _as_float_array(kdp)
    gate_length = float(dr)
    valid = _prepare_mask(kdp, mask=mask)
    flat_kdp = kdp.reshape(-1, kdp.shape[-1])
    flat_valid = valid.reshape(-1, valid.shape[-1])
    flat_pia = np.full_like(flat_kdp, np.nan)

    for iray, (row, row_valid) in enumerate(zip(flat_kdp, flat_valid)):
        running = 0.0
        for igate, (value, is_valid) in enumerate(zip(row, row_valid)):
            if not is_valid:
                running = 0.0
                continue
            if clip_negative:
                value = max(value, 0.0)
            running += 2.0 * gamma * value * gate_length
            flat_pia[iray, igate] = running

    return flat_pia.reshape(kdp.shape)


def correct_attenuation_kdp(ref, kdp, dr, zdr=None, gamma=0.08, beta=0.02, mask=None):
    """
    Dual-polarization attenuation correction driven by KDP.

    Parameters
    ----------
    ref : array-like
        Reflectivity in dBZ.
    kdp : array-like
        Specific differential phase in deg/km.
    dr : float
        Gate length in km.
    zdr : array-like, optional
        Differential reflectivity in dB.
    gamma, beta : float
        Linear coefficients converting KDP to specific attenuation for Z and
        ZDR respectively.
    mask : array-like, optional
        Boolean mask of gates to be corrected.
    """
    reflectivity = _as_float_array(ref)
    kdp = _as_float_array(kdp)
    valid_ref = np.isfinite(reflectivity)
    valid_kdp = np.isfinite(kdp)
    correction_mask = valid_ref & valid_kdp
    if mask is not None:
        correction_mask &= np.asarray(mask, dtype=bool)

    pia = pia_from_kdp(kdp, dr=dr, gamma=gamma, mask=correction_mask, clip_negative=True)
    z_corr = np.where(correction_mask, reflectivity + pia, reflectivity)

    if zdr is None:
        zdr_corr = np.full_like(reflectivity, np.nan)
        pia_zdr = np.full_like(reflectivity, np.nan)
    else:
        differential_reflectivity = _as_float_array(zdr)
        valid_zdr = np.isfinite(differential_reflectivity)
        zdr_mask = correction_mask & valid_zdr
        pia_zdr = pia_from_kdp(kdp, dr=dr, gamma=beta, mask=zdr_mask, clip_negative=True)
        zdr_corr = np.where(zdr_mask, differential_reflectivity + pia_zdr, differential_reflectivity)

    return z_corr, pia, zdr_corr, pia_zdr
