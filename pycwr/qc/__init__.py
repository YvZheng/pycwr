from importlib import import_module

__all__ = [
    "apply_dualpol_qc",
    "build_clear_air_mask",
    "build_meteo_mask",
    "correct_attenuation",
    "correct_attenuation_HB",
    "correct_attenuation_kdp",
    "despeckle_mask",
    "kdp_from_phidp",
    "phidp_texture",
    "pia_from_kdp",
    "QC_REFERENCE_NOTES",
    "resolve_kdp_coefficients",
    "run_dualpol_qc",
    "smooth_phidp",
]

_LAZY_SYMBOLS = {
    "correct_attenuation": ("attenuation", "correct_attenuation"),
    "correct_attenuation_HB": ("attenuation", "correct_attenuation_HB"),
    "correct_attenuation_kdp": ("attenuation", "correct_attenuation_kdp"),
    "pia_from_kdp": ("attenuation", "pia_from_kdp"),
    "resolve_kdp_coefficients": ("attenuation", "resolve_kdp_coefficients"),
    "apply_dualpol_qc": ("pipeline", "apply_dualpol_qc"),
    "run_dualpol_qc": ("pipeline", "run_dualpol_qc"),
    "QC_REFERENCE_NOTES": ("polarimetric", "QC_REFERENCE_NOTES"),
    "build_clear_air_mask": ("polarimetric", "build_clear_air_mask"),
    "build_meteo_mask": ("polarimetric", "build_meteo_mask"),
    "despeckle_mask": ("polarimetric", "despeckle_mask"),
    "kdp_from_phidp": ("polarimetric", "kdp_from_phidp"),
    "phidp_texture": ("polarimetric", "phidp_texture"),
    "smooth_phidp": ("polarimetric", "smooth_phidp"),
}


def __getattr__(name):
    symbol = _LAZY_SYMBOLS.get(name)
    if symbol is not None:
        module = import_module(f".{symbol[0]}", __name__)
        value = getattr(module, symbol[1])
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
