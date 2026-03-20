# -*- coding: utf-8 -*-
"""
Dual-polarization QC pipeline.
"""

import copy

import numpy as np
import xarray as xr

from ..configure.default_config import CINRAD_field_mapping, DEFAULT_METADATA
from .attenuation import correct_attenuation_kdp, resolve_kdp_coefficients
from .polarimetric import (
    QC_REFERENCE_NOTES,
    build_clear_air_mask,
    build_meteo_mask,
    despeckle_mask,
    kdp_from_phidp,
    phidp_texture,
    smooth_phidp,
)


def _empty_like(reference):
    return np.full_like(np.asanyarray(reference, dtype=float), np.nan, dtype=float)


def _field_metadata(field_name):
    if field_name in CINRAD_field_mapping:
        metadata_key = CINRAD_field_mapping[field_name]
        if metadata_key in DEFAULT_METADATA:
            return dict(DEFAULT_METADATA[metadata_key])
    return {}


def _build_data_array(template, values, field_name):
    data_array = xr.DataArray(
        np.asanyarray(values),
        dims=template.dims,
        coords=template.coords,
        name=field_name,
    )
    data_array.attrs = _field_metadata(field_name)
    return data_array


def run_dualpol_qc(
    ref,
    zdr=None,
    phidp=None,
    kdp=None,
    rhohv=None,
    snr=None,
    dr=0.075,
    band="C",
    use_existing_kdp=True,
    clear_air_mode="label",
    clear_air_max_ref=15.0,
    clear_air_max_rhohv=0.97,
    clear_air_max_phidp_texture=10.0,
    clear_air_max_snr=20.0,
):
    """
    Run a classical dual-polarization QC chain on one sweep array set.
    """
    reflectivity = np.asanyarray(ref, dtype=float)
    if reflectivity.ndim < 2:
        raise ValueError("ref must have at least two dimensions with range on the last axis")

    gate_length = float(dr)
    gamma, beta = resolve_kdp_coefficients(band=band)
    base_mask = np.isfinite(reflectivity)

    if phidp is not None:
        phase = np.asanyarray(phidp, dtype=float)
        phidp_smooth = smooth_phidp(phase, mask=base_mask)
        texture = phidp_texture(phidp_smooth)
    else:
        phase = None
        phidp_smooth = _empty_like(reflectivity)
        texture = _empty_like(reflectivity)

    if clear_air_mode not in {"ignore", "label", "mask"}:
        raise ValueError("clear_air_mode must be one of 'ignore', 'label', or 'mask'")

    meteo_mask = build_meteo_mask(
        reflectivity,
        rhohv=rhohv,
        phidp_texture=None if phase is None else texture,
        snr=snr,
    )
    if clear_air_mode == "ignore":
        clear_air_mask = np.zeros_like(reflectivity, dtype=bool)
    else:
        clear_air_mask = build_clear_air_mask(
            reflectivity,
            rhohv=rhohv,
            phidp_texture=None if phase is None else texture,
            snr=snr,
            max_ref=clear_air_max_ref,
            max_rhohv=clear_air_max_rhohv,
            max_phidp_texture=clear_air_max_phidp_texture,
            max_snr=clear_air_max_snr,
        )
    accepted_mask = meteo_mask & base_mask
    if clear_air_mode == "mask":
        accepted_mask &= ~clear_air_mask
    qc_mask = despeckle_mask(accepted_mask, min_size=9, connectivity=8)

    if kdp is not None and use_existing_kdp:
        kdp_used = np.where(qc_mask, np.asanyarray(kdp, dtype=float), np.nan)
    elif phase is not None:
        kdp_used = np.where(qc_mask, kdp_from_phidp(phidp_smooth, dr=gate_length), np.nan)
    else:
        raise ValueError("dual-pol QC requires either kdp or phidp")

    z_corr, pia, zdr_corr, pia_zdr = correct_attenuation_kdp(
        reflectivity,
        kdp_used,
        dr=gate_length,
        zdr=zdr,
        gamma=gamma,
        beta=beta,
        mask=qc_mask,
    )

    z_corr = np.where(qc_mask, z_corr, np.nan)
    pia = np.where(qc_mask, pia, np.nan)
    kdp_used = np.where(qc_mask, kdp_used, np.nan)

    if zdr is None:
        zdr_corr = _empty_like(reflectivity)
        pia_zdr = _empty_like(reflectivity)
    else:
        zdr_corr = np.where(qc_mask, zdr_corr, np.nan)
        pia_zdr = np.where(qc_mask, pia_zdr, np.nan)

    return {
        "ref_corrected": z_corr,
        "zdr_corrected": zdr_corr,
        "pia": pia,
        "pia_zdr": pia_zdr,
        "phidp_smooth": phidp_smooth,
        "kdp_used": kdp_used,
        "phidp_texture": texture,
        "meteo_mask": meteo_mask,
        "clear_air_mask": clear_air_mask,
        "qc_mask": qc_mask,
        "qc_reference_notes": dict(QC_REFERENCE_NOTES),
    }


def apply_dualpol_qc(
    prd,
    sweeps=None,
    inplace=False,
    band="C",
    use_existing_kdp=True,
    clear_air_mode="label",
    clear_air_max_ref=15.0,
    clear_air_max_rhohv=0.97,
    clear_air_max_phidp_texture=10.0,
    clear_air_max_snr=20.0,
):
    """
    Apply the dual-polarization QC chain to selected sweeps of a PRD object.
    """
    target = prd if inplace else copy.deepcopy(prd)
    if sweeps is None:
        sweeps = list(range(target.nsweeps))

    for sweep in sweeps:
        sweep = int(sweep)
        dataset = target.fields[sweep]
        if "dBZ" not in dataset:
            raise KeyError("Sweep %d is missing dBZ" % sweep)
        if "KDP" not in dataset and "PhiDP" not in dataset:
            raise ValueError("Sweep %d must contain KDP or PhiDP for dual-pol QC" % sweep)

        reference_field = dataset["dBZ"]
        range_values = np.asanyarray(dataset["range"].values, dtype=float)
        if range_values.size < 2:
            raise ValueError("Sweep %d must contain at least two range gates" % sweep)
        gate_length = float(np.median(np.diff(range_values)) / 1000.0)

        results = run_dualpol_qc(
            reference_field.values,
            zdr=None if "ZDR" not in dataset else dataset["ZDR"].values,
            phidp=None if "PhiDP" not in dataset else dataset["PhiDP"].values,
            kdp=None if "KDP" not in dataset else dataset["KDP"].values,
            rhohv=None if "CC" not in dataset else dataset["CC"].values,
            snr=None if "SNRH" not in dataset else dataset["SNRH"].values,
            dr=gate_length,
            band=band,
            use_existing_kdp=use_existing_kdp,
            clear_air_mode=clear_air_mode,
            clear_air_max_ref=clear_air_max_ref,
            clear_air_max_rhohv=clear_air_max_rhohv,
            clear_air_max_phidp_texture=clear_air_max_phidp_texture,
            clear_air_max_snr=clear_air_max_snr,
        )

        dataset["Zc"] = _build_data_array(reference_field, results["ref_corrected"], "Zc")
        dataset["PIA"] = _build_data_array(reference_field, results["pia"], "PIA")
        dataset["PhiDP_smooth"] = _build_data_array(reference_field, results["phidp_smooth"], "PhiDP_smooth")
        dataset["PhiDP_texture"] = _build_data_array(reference_field, results["phidp_texture"], "PhiDP_texture")
        dataset["KDPc"] = _build_data_array(reference_field, results["kdp_used"], "KDPc")
        dataset["METEO_MASK"] = _build_data_array(reference_field, results["meteo_mask"], "METEO_MASK")
        dataset["CLEAR_AIR_MASK"] = _build_data_array(reference_field, results["clear_air_mask"], "CLEAR_AIR_MASK")
        dataset["QC_MASK"] = _build_data_array(reference_field, results["qc_mask"], "QC_MASK")

        if "ZDR" in dataset:
            zdr_template = dataset["ZDR"]
        else:
            zdr_template = reference_field
        dataset["ZDRc"] = _build_data_array(zdr_template, results["zdr_corrected"], "ZDRc")
        dataset["PIA_ZDR"] = _build_data_array(zdr_template, results["pia_zdr"], "PIA_ZDR")

    if hasattr(target, "_invalidate_cached_views"):
        target._invalidate_cached_views()
    return target
