from importlib import import_module

__all__ = [
    "WIND_REFERENCE_NOTES",
    "VAD",
    "VVP",
    "vad",
    "vvp",
    "retrieve_vad",
    "retrieve_vvp",
    "retrieve_vwp",
    "retrieve_wind_volume_xy",
    "retrieve_wind_volume_lonlat",
    "select_velocity_field",
    "fhc_HCL",
    "fhc_hcl",
    "hid_beta_function",
    "classify_hydrometeors",
    "classify_sweep_hydrometeors",
    "apply_hydrometeor_classification",
    "interpolate_temperature_profile",
    "available_hydrometeor_classes",
    "hydrometeor_class_name",
    "HYDROMETEOR_CLASSES",
    "HYDROMETEOR_CLASS_IDS",
    "HID_REFERENCE_NOTES",
]


def __getattr__(name):
    if name in {
        "VAD",
        "VVP",
        "vad",
        "vvp",
        "retrieve_vad",
        "retrieve_vvp",
        "retrieve_vwp",
        "retrieve_wind_volume_xy",
        "retrieve_wind_volume_lonlat",
        "select_velocity_field",
        "WIND_REFERENCE_NOTES",
    }:
        module = import_module(".WindField", __name__)
        value = getattr(module, name)
        globals()[name] = value
        return value
    if name in {
        "fhc_HCL",
        "fhc_hcl",
        "hid_beta_function",
        "classify_hydrometeors",
        "classify_sweep_hydrometeors",
        "apply_hydrometeor_classification",
        "interpolate_temperature_profile",
        "available_hydrometeor_classes",
        "hydrometeor_class_name",
        "HYDROMETEOR_CLASSES",
        "HYDROMETEOR_CLASS_IDS",
        "HID_REFERENCE_NOTES",
    }:
        module = import_module(".HID", __name__)
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
