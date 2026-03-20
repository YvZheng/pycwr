# -*- coding: utf-8 -*-
"""
Classical polarimetric quality-control primitives.
"""

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy import ndimage

QC_REFERENCE_NOTES = {
    "wu2018": "吴翀;双偏振雷达的资料质量分析,相态识別及组网应用[D];南京信息工程大学;2018年。第2章将晴空回波视为应在QC中分离处理的非气象回波类型。",
    "wdtd_dpqc": "NOAA WDTD dual-pol QC guidance removes non-weather echoes using polarimetric signatures before downstream products are derived.",
    "tang2020": "Tang et al. (2020, JTECH) describe MRMS quality control as a weather / non-weather separation problem with dedicated weak-echo handling.",
}


def _as_float_array(data):
    return np.asanyarray(data, dtype=float)


def _normalize_window(window):
    window = int(window)
    if window < 1:
        raise ValueError("window must be a positive integer")
    if window % 2 == 0:
        window += 1
    return window


def _pad_last_axis(data, pad):
    pad_width = [(0, 0)] * data.ndim
    pad_width[-1] = (pad, pad)
    return np.pad(data, pad_width, mode="constant", constant_values=np.nan)


def _nanmedian_last_axis(data, window):
    window = _normalize_window(window)
    pad = window // 2
    padded = _pad_last_axis(data, pad)
    windows = sliding_window_view(padded, window_shape=window, axis=-1)
    return np.nanmedian(windows, axis=-1)


def _nanmean_last_axis(data, window):
    window = _normalize_window(window)
    pad = window // 2
    padded = _pad_last_axis(data, pad)
    windows = sliding_window_view(padded, window_shape=window, axis=-1)
    return np.nanmean(windows, axis=-1)


def smooth_phidp(phidp, mask=None, median_window=5, fit_window=7, enforce_monotonic=True):
    """
    Smooth PhiDP with median de-spiking and local averaging.
    """
    phase = _as_float_array(phidp)
    valid = np.isfinite(phase)
    if mask is not None:
        valid &= np.asarray(mask, dtype=bool)
    phase = np.where(valid, phase, np.nan)
    phase = _nanmedian_last_axis(phase, median_window)
    phase = _nanmean_last_axis(phase, fit_window)

    if enforce_monotonic:
        flat_phase = phase.reshape(-1, phase.shape[-1])
        for row in flat_phase:
            running = np.nan
            for igate, value in enumerate(row):
                if not np.isfinite(value):
                    running = np.nan
                    continue
                if np.isnan(running) or value >= running:
                    running = value
                else:
                    row[igate] = running
        phase = flat_phase.reshape(phase.shape)

    return phase


def kdp_from_phidp(phidp_smooth, dr, fit_window=7):
    """
    Estimate KDP from a smoothed PhiDP field using a local linear fit.
    """
    phase = _as_float_array(phidp_smooth)
    fit_window = _normalize_window(fit_window)
    pad = fit_window // 2
    padded = _pad_last_axis(phase, pad)
    windows = sliding_window_view(padded, window_shape=fit_window, axis=-1)

    x = np.arange(fit_window, dtype=float) - pad
    valid = np.isfinite(windows)
    count = valid.sum(axis=-1).astype(float)
    x_masked = valid * x
    x_mean = np.divide(
        x_masked.sum(axis=-1),
        count,
        out=np.zeros_like(count),
        where=count > 0,
    )
    y_sum = np.nansum(windows, axis=-1)
    y_mean = np.divide(
        y_sum,
        count,
        out=np.zeros_like(count),
        where=count > 0,
    )
    x_centered = x - x_mean[..., np.newaxis]
    y_centered = windows - y_mean[..., np.newaxis]
    covariance = np.nansum(x_centered * y_centered * valid, axis=-1)
    variance = np.sum((x_centered ** 2) * valid, axis=-1)
    slope = np.divide(
        covariance,
        variance,
        out=np.full_like(covariance, np.nan),
        where=(variance > 0) & (count >= 2),
    )
    return 0.5 * slope / float(dr)


def phidp_texture(phidp, window=7):
    """
    Compute a simple PhiDP texture field from local standard deviation.
    """
    phase = _as_float_array(phidp)
    window = _normalize_window(window)
    pad = window // 2
    padded = _pad_last_axis(phase, pad)
    windows = sliding_window_view(padded, window_shape=window, axis=-1)
    return np.nanstd(windows, axis=-1)


def build_meteo_mask(
    ref,
    rhohv=None,
    phidp_texture=None,
    snr=None,
    min_ref=0.0,
    min_rhohv=0.85,
    max_phidp_texture=20.0,
    min_snr=3.0,
):
    """
    Build a conservative meteorological echo mask from classic thresholds.
    """
    reflectivity = _as_float_array(ref)
    mask = np.isfinite(reflectivity) & (reflectivity >= min_ref)

    if rhohv is not None:
        rhohv = _as_float_array(rhohv)
        mask &= np.isfinite(rhohv) & (rhohv >= min_rhohv)

    if phidp_texture is not None:
        texture = _as_float_array(phidp_texture)
        mask &= np.isfinite(texture) & (texture <= max_phidp_texture)

    if snr is not None:
        snr = _as_float_array(snr)
        mask &= np.isfinite(snr) & (snr >= min_snr)

    return mask


def build_clear_air_mask(
    ref,
    rhohv=None,
    phidp_texture=None,
    snr=None,
    max_ref=15.0,
    max_rhohv=0.97,
    max_phidp_texture=10.0,
    max_snr=20.0,
):
    """
    Build a conservative likely-clear-air echo mask.

    This is not a full hydrometeor classification algorithm. It is a targeted
    weak-echo diagnostic intended to separate likely clear-air / biological
    echoes from precipitation-focused QC workflows.
    """
    reflectivity = _as_float_array(ref)
    mask = np.isfinite(reflectivity) & (reflectivity <= max_ref)
    diagnostic_terms = 0

    if rhohv is not None and max_rhohv is not None:
        rhohv = _as_float_array(rhohv)
        mask &= np.isfinite(rhohv) & (rhohv <= max_rhohv)
        diagnostic_terms += 1

    if phidp_texture is not None and max_phidp_texture is not None:
        texture = _as_float_array(phidp_texture)
        mask &= np.isfinite(texture) & (texture <= max_phidp_texture)
        diagnostic_terms += 1

    if snr is not None and max_snr is not None:
        snr = _as_float_array(snr)
        mask &= np.isfinite(snr) & (snr <= max_snr)
        diagnostic_terms += 1

    if diagnostic_terms == 0:
        return np.zeros_like(reflectivity, dtype=bool)
    return mask


def despeckle_mask(mask, min_size=9, connectivity=8):
    """
    Remove small connected components from a boolean mask.
    """
    binary_mask = np.asarray(mask, dtype=bool)
    if not np.any(binary_mask):
        return binary_mask.copy()

    if connectivity not in (4, 8):
        raise ValueError("connectivity must be 4 or 8")
    structure_order = 1 if connectivity == 4 else binary_mask.ndim
    structure = ndimage.generate_binary_structure(binary_mask.ndim, structure_order)
    labels, num_labels = ndimage.label(binary_mask, structure=structure)
    if num_labels == 0:
        return binary_mask.copy()
    label_sizes = np.bincount(labels.ravel())
    keep = label_sizes >= int(min_size)
    keep[0] = False
    return keep[labels]
