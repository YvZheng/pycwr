# -*- coding: utf-8 -*-
"""
Hydrometeor identification (HID / HCL) retrieval helpers.

This module packages a fuzzy-logic hydrometeor classifier around the
membership beta functions shipped with pycwr. The bundled membership table
supports ten warm-season hydrometeor classes:

1. Drizzle
2. Rain
3. Ice Crystals
4. Dry Aggregates Snow
5. Wet Snow
6. Vertical Ice
7. Low-Density Graupel
8. High-Density Graupel
9. Hail
10. Big Drops

The full classifier can use an external temperature field or sounding profile.
When no thermodynamic profile is available, the module falls back to a
reduced-variable fuzzy scheme by setting the temperature contribution to zero.
"""

import copy
import json

import numpy as np
import xarray as xr

from ..configure.default_config import CINRAD_field_mapping, DEFAULT_METADATA


HYDROMETEOR_CLASSES = (
    "Drizzle",
    "Rain",
    "Ice Crystals",
    "Dry Aggregates Snow",
    "Wet Snow",
    "Vertical Ice",
    "Low-Density Graupel",
    "High-Density Graupel",
    "Hail",
    "Big Drops",
)
HYDROMETEOR_CLASS_IDS = {name: index + 1 for index, name in enumerate(HYDROMETEOR_CLASSES)}
HYDROMETEOR_CLASS_MEANINGS = (
    "drizzle rain ice_crystals dry_aggregates_snow wet_snow vertical_ice "
    "low_density_graupel high_density_graupel hail big_drops"
)
HID_REFERENCE_NOTES = [
    "Dolan, B., S. A. Rutledge, S. Lim, V. Chandrasekar, and M. Thurai, 2013: "
    "A robust C-band hydrometeor identification algorithm and application to a long-term "
    "polarimetric radar dataset. Journal of Applied Meteorology and Climatology, 52, 2162-2186.",
    "Marzano, F. S., D. Scaranari, M. Celano, P. P. Alberoni, G. Vulpiani, and M. Montopoli, 2006: "
    "Hydrometeor classification from dual-polarized weather radar: extending fuzzy logic from "
    "S-band to C-band data. Advances in Geosciences, 7, 109-114.",
]
DEFAULT_WEIGHTS = {"dBZ": 1.5, "ZDR": 0.8, "KDP": 1.0, "CC": 0.8, "LDR": 0.5, "T": 0.4}

_SUPPORTED_BANDS = ("S", "C", "X")
_SUPPORTED_FEATURES = ("dBZ", "ZDR", "KDP", "CC", "LDR", "T")
_POLAR_FEATURES = ("LDR", "KDP", "ZDR", "CC")
_DEFAULT_SWEEP_FIELDS = {
    "dBZ": ("Zc", "dBZ"),
    "ZDR": ("ZDRc", "ZDR"),
    "KDP": ("KDPc", "KDP"),
    "CC": ("CC",),
    "LDR": ("LDR",),
}
_BETA_PARAMETER_TABLE = json.loads(
    r"""{"S": {"dBZ": [[2.0, 29.0, 10.0], [41.5, 15.5, 10.0], [-3.0, 22.0, 20.0], [17.0, 17.0, 15.0], [21.0, 22.0, 10.0], [-3.0, 22.0, 20.0], [37.0, 8.0, 8.0], [49.0, 9.0, 6.0], [58.0, 12.0, 10.0], [57.0, 9.0, 10.0]], "ZDR": [[0.3499999940395355, 0.3499999940395355, 5.0], [2.5999999046325684, 2.799999952316284, 9.0], [3.200000047683716, 2.799999952316284, 10.0], [0.6000000238418579, 0.6000000238418579, 7.0], [1.2999999523162842, 1.2999999523162842, 10.0], [-0.8999999761581421, 0.8999999761581421, 10.0], [0.30000001192092896, 0.800000011920929, 6.0], [1.0, 1.899999976158142, 8.0], [0.14000000059604645, 0.550000011920929, 8.0], [4.0, 1.7000000476837158, 8.0]], "KDP": [[0.009999999776482582, 0.009999999776482582, 2.0], [3.700000047683716, 4.0, 10.0], [0.0430000014603138, 0.0430000014603138, 6.0], [0.03999999910593033, 0.05000000074505806, 1.0], [0.12999999523162842, 0.5, 6.0], [-0.23000000417232513, 0.23000000417232513, 3.0], [0.20000000298023224, 0.5600000023841858, 3.0], [0.550000011920929, 1.1549999713897705, 3.0], [0.20000000298023224, 0.800000011920929, 6.0], [1.600000023841858, 1.5, 6.0]], "CC": [[1.0, 0.014999999664723873, 3.0], [1.0, 0.019999999552965164, 2.0], [1.0, 0.019999999552965164, 3.0], [0.9980000257492065, 0.019999999552965164, 3.0], [0.7799999713897705, 0.20000000298023224, 10.0], [0.9700000286102295, 0.03999999910593033, 3.0], [1.0, 0.009999999776482582, 1.0], [1.0, 0.03999999910593033, 4.0], [0.9599999785423279, 0.10000000149011612, 3.0], [0.9800000190734863, 0.03999999910593033, 3.0]], "LDR": [[-48.380001068115234, 5.5, 10.0], [-25.0, 5.5, 4.0], [-23.329999923706055, 7.980000019073486, 20.0], [-46.45500183105469, 19.674999237060547, 20.0], [-13.604999542236328, 6.445000171661377, 8.0], [-29.68000030517578, 11.899999618530273, 20.0], [-43.18000030517578, 12.920000076293945, 8.0], [-35.89500045776367, 13.204999923706055, 3.0], [-22.2450008392334, 7.764999866485596, 8.0], [-34.08000183105469, 3.059999942779541, 10.0]], "T": [[40.0, 41.0, 50.0], [48.0, 51.0, 30.0], [-50.0, 50.0, 25.0], [-25.0, 26.0, 15.0], [1.0, 3.5, 5.0], [-50.0, 50.0, 25.0], [-50.0, 50.0, 25.0], [-2.5, 20.0, 2.0], [0.0, 100.0, 5.0], [48.0, 51.0, 30.0]]}, "C": {"dBZ": [[1.75, 29.0, 10.0], [39.0, 19.0, 10.0], [-3.0, 22.0, 20.0], [17.0, 18.100000381469727, 10.0], [24.0, 21.299999237060547, 10.0], [-3.0, 22.0, 20.0], [37.0, 9.199999809265137, 8.0], [44.29999923706055, 10.199999809265137, 6.0], [62.29999923706055, 14.300000190734863, 10.0], [57.79999923706055, 8.5, 10.0]], "ZDR": [[0.46000000834465027, 0.46000000834465027, 5.0], [2.299999952316284, 2.200000047683716, 9.0], [2.9000000953674316, 2.700000047683716, 10.0], [1.0, 1.100000023841858, 7.0], [1.2999999523162842, 0.8999999761581421, 10.0], [-0.8999999761581421, 0.8999999761581421, 10.0], [0.8999999761581421, 0.8999999761581421, 6.0], [1.600000023841858, 1.2000000476837158, 3.0], [0.14000000059604645, 0.5600000023841858, 8.0], [4.400000095367432, 1.899999976158142, 8.0]], "KDP": [[0.029999999329447746, 0.029999999329447746, 2.0], [5.5, 5.5, 10.0], [0.07999999821186066, 0.07999999821186066, 6.0], [-0.00800000037997961, 0.30000001192092896, 1.0], [0.25, 0.4300000071525574, 6.0], [-0.75, 0.75, 30.0], [0.10000000149011612, 0.07999999821186066, 3.0], [1.899999976158142, 1.899999976158142, 3.0], [0.6000000238418579, 3.5, 6.0], [3.4000000953674316, 3.299999952316284, 6.0]], "CC": [[1.0, 0.017999999225139618, 3.0], [1.0, 0.02500000037252903, 3.0], [0.9800000190734863, 0.02500000037252903, 3.0], [0.9300000071525574, 0.07000000029802322, 3.0], [0.7400000095367432, 0.25, 10.0], [0.9750000238418579, 0.02199999988079071, 3.0], [1.0, 0.02500000037252903, 1.0], [1.0, 0.03999999910593033, 2.0], [0.9700000286102295, 0.10000000149011612, 3.0], [0.9900000095367432, 0.029999999329447746, 3.0]], "LDR": [[-48.33000183105469, 5.5, 10.0], [-25.0, 5.5, 4.0], [-30.059999465942383, 14.630000114440918, 20.0], [-34.400001525878906, 25.600000381469727, 20.0], [-13.194999694824219, 6.855000019073486, 8.0], [-30.700000762939453, 12.899999618530273, 20.0], [-50.84000015258789, 24.709999084472656, 8.0], [-44.584999084472656, 22.565000534057617, 3.0], [-22.239999771118164, 7.039999961853027, 8.0], [-33.279998779296875, 3.130000114440918, 10.0]], "T": [[40.0, 41.0, 50.0], [48.0, 51.0, 30.0], [-50.0, 50.0, 25.0], [-25.0, 26.0, 15.0], [1.0, 3.5, 5.0], [-50.0, 50.0, 25.0], [-50.0, 50.0, 25.0], [-2.5, 20.0, 2.0], [0.0, 100.0, 5.0], [48.0, 51.0, 30.0]]}, "X": {"dBZ": [[2.0, 29.0, 10.0], [42.0, 17.0, 10.0], [-3.0, 22.0, 20.0], [16.0, 17.0, 15.0], [28.0, 29.0, 10.0], [-3.0, 22.0, 20.0], [36.0, 9.0, 8.0], [43.0, 9.0, 5.0], [59.0, 10.300000190734863, 6.0], [56.0, 9.0, 10.0]], "ZDR": [[0.44999998807907104, 0.44999998807907104, 5.0], [2.8499999046325684, 3.0199999809265137, 9.0], [3.200000047683716, 2.5999999046325684, 10.0], [0.699999988079071, 0.699999988079071, 7.0], [1.600000023841858, 1.2000000476837158, 10.0], [-0.8999999761581421, 1.0, 10.0], [0.30000001192092896, 1.0, 6.0], [0.6000000238418579, 1.5, 2.0], [0.11999999731779099, 0.5, 8.0], [4.300000190734863, 2.0, 8.0]], "KDP": [[0.029999999329447746, 0.029999999329447746, 2.0], [12.75, 13.0, 60.0], [0.15000000596046448, 0.15000000596046448, 6.0], [0.20000000298023224, 0.20000000298023224, 1.0], [-7.0, 8.100000381469727, 6.0], [-0.7900000214576721, 0.7900000214576721, 3.0], [0.699999988079071, 2.0999999046325684, 3.0], [1.5499999523162842, 3.049999952316284, 3.0], [0.5, 1.5, 6.0], [3.5999999046325684, 3.299999952316284, 6.0]], "CC": [[1.0, 0.014999999664723873, 2.0], [1.0, 0.03999999910593033, 2.0], [1.0, 0.029999999329447746, 3.0], [0.9980000257492065, 0.019999999552965164, 3.0], [0.699999988079071, 0.25999999046325684, 10.0], [0.9800000190734863, 0.019999999552965164, 3.0], [1.0, 0.014999999664723873, 1.0], [1.0, 0.03500000014901161, 1.0], [0.9700000286102295, 0.15000000596046448, 3.0], [0.9800000190734863, 0.029999999329447746, 3.0]], "LDR": [[-48.33000183105469, 5.5, 10.0], [-25.090999603271484, 5.5, 4.0], [-23.44499969482422, 8.104999542236328, 20.0], [-45.935001373291016, 19.434999465942383, 20.0], [-13.225000381469727, 6.815000057220459, 8.0], [-31.084999084472656, 13.300000190734863, 20.0], [-42.32500076293945, 13.404999732971191, 8.0], [-34.34000015258789, 13.65999984741211, 3.0], [-22.475000381469727, 6.764999866485596, 8.0], [-32.58000183105469, 0.9399999976158142, 10.0]], "T": [[40.0, 41.0, 50.0], [48.0, 51.0, 30.0], [-50.0, 50.0, 25.0], [-25.0, 26.0, 15.0], [1.0, 3.5, 5.0], [-50.0, 50.0, 25.0], [-50.0, 50.0, 25.0], [-2.5, 20.0, 2.0], [0.0, 100.0, 5.0], [48.0, 51.0, 30.0]]}}"""
)


def _get_beta_params():
    return _BETA_PARAMETER_TABLE


def _normalize_band(band):
    band_key = str(band).strip().upper()
    if band_key not in _SUPPORTED_BANDS:
        raise ValueError("band must be one of %s" % (", ".join(_SUPPORTED_BANDS)))
    return band_key


def _normalize_method(method):
    method_key = str(method).strip().lower()
    if method_key not in {"hybrid", "linear"}:
        raise ValueError("method must be 'hybrid' or 'linear'")
    return method_key


def _normalize_weights(weights):
    merged = dict(DEFAULT_WEIGHTS)
    if weights:
        merged.update(weights)
    normalized = {}
    for feature in _SUPPORTED_FEATURES:
        value = float(merged.get(feature, 0.0))
        normalized[feature] = max(value, 0.0)
    return normalized


def _as_float_array(value, name, shape=None):
    if value is None:
        return None
    array = np.asarray(value, dtype=np.float64)
    if shape is None:
        return array
    if array.shape == shape:
        return array
    try:
        return np.broadcast_to(array, shape).astype(np.float64, copy=False)
    except ValueError as exc:
        raise ValueError("%s must be broadcastable to %s, got %s" % (name, shape, array.shape)) from exc


def _feature_beta(feature_name, values, band):
    params = np.asarray(_get_beta_params()[_normalize_band(band)][feature_name], dtype=np.float64)
    return hid_beta_function(
        values,
        params[:, 0],
        params[:, 1],
        params[:, 2],
    )


def _gate_valid_weight(values, weight):
    if values is None or weight <= 0.0:
        return None
    return np.where(np.isfinite(values), float(weight), 0.0)


def _infer_profile_altitude(profile_height, radar_altitude, height_reference):
    profile_height = np.asarray(profile_height, dtype=np.float64)
    if profile_height.ndim != 1:
        raise ValueError("profile_height must be a 1-D array")
    reference = str(height_reference).strip().lower()
    if reference not in {"asl", "msl", "agl"}:
        raise ValueError("height_reference must be 'asl', 'msl', or 'agl'")
    if reference == "agl":
        return profile_height + float(radar_altitude)
    return profile_height


def _field_metadata(field_name):
    if field_name in CINRAD_field_mapping:
        metadata_key = CINRAD_field_mapping[field_name]
        if metadata_key in DEFAULT_METADATA:
            return dict(DEFAULT_METADATA[metadata_key])
    return {}


def _build_data_array(template, values, name, attrs=None):
    data_array = xr.DataArray(
        np.asanyarray(values),
        dims=template.dims,
        coords=template.coords,
        name=name,
    )
    metadata = _field_metadata(name)
    if attrs:
        metadata.update(attrs)
    data_array.attrs = metadata
    return data_array


def _temperature_field_attrs(source, method):
    return {
        "units": "degC",
        "standard_name": "interpolated_profile",
        "long_name": "Hydrometeor classification temperature field",
        "comment": "Temperature field used by hydrometeor classification (%s)." % source,
        "retrieval_method": method,
        "references": json.dumps(HID_REFERENCE_NOTES, ensure_ascii=False),
    }


def _confidence_field_attrs(output_field, method):
    return {
        "units": "unitless",
        "standard_name": "hydrometeor_classification_confidence",
        "long_name": "Hydrometeor classification confidence",
        "comment": "Maximum fuzzy-membership score associated with %s." % output_field,
        "retrieval_method": method,
        "valid_min": 0.0,
        "valid_max": 1.0,
        "references": json.dumps(HID_REFERENCE_NOTES, ensure_ascii=False),
    }


def _hcl_field_attrs(method, band, used_temperature):
    attrs = _field_metadata("HCL")
    attrs.update(
        {
            "units": "class_id",
            "long_name": "Hydrometeor classification",
            "flag_values": list(range(1, len(HYDROMETEOR_CLASSES) + 1)),
            "flag_meanings": HYDROMETEOR_CLASS_MEANINGS,
            "hydrometeor_classes": json.dumps(list(HYDROMETEOR_CLASSES), ensure_ascii=False),
            "references": json.dumps(HID_REFERENCE_NOTES, ensure_ascii=False),
            "retrieval_method": method,
            "radar_band": _normalize_band(band),
            "uses_temperature": "true" if used_temperature else "false",
            "comment": (
                "Fuzzy-logic hydrometeor classification. "
                "When temperature is unavailable, pycwr falls back to a reduced-variable "
                "polarimetric scheme by omitting the temperature term."
            ),
        }
    )
    return attrs


def _select_dataset_field(dataset, candidates, required=False):
    for name in candidates:
        if name in dataset:
            return dataset[name]
    if required:
        raise KeyError("Missing required field: %s" % "/".join(candidates))
    return None


def _resolve_temperature_spec(spec, sweep, total_sweeps):
    if spec is None:
        return None
    if isinstance(spec, dict):
        return spec.get(int(sweep))
    if isinstance(spec, (list, tuple)) and len(spec) == int(total_sweeps):
        return spec[int(sweep)]
    return spec


def available_hydrometeor_classes():
    """Return the supported hydrometeor class names."""
    return HYDROMETEOR_CLASSES


def hydrometeor_class_name(class_id):
    """Return the class label for one hydrometeor class id."""
    index = int(class_id)
    if index < 1 or index > len(HYDROMETEOR_CLASSES):
        raise ValueError("class_id must be between 1 and %d" % len(HYDROMETEOR_CLASSES))
    return HYDROMETEOR_CLASSES[index - 1]


def interpolate_temperature_profile(
    gate_altitude,
    profile_height,
    profile_temperature,
    radar_altitude=0.0,
    height_reference="asl",
):
    """
    Interpolate a 1-D environmental temperature profile onto radar gate heights.

    Parameters
    ----------
    gate_altitude : array-like
        Gate altitude field in meters above sea level, typically ``dataset['z']``.
    profile_height : 1-D array-like
        Profile heights in meters. Use ``height_reference='agl'`` when these
        heights are relative to the radar.
    profile_temperature : 1-D array-like
        Profile temperatures in degrees Celsius.
    radar_altitude : float, optional
        Radar altitude in meters above sea level.
    height_reference : {"asl", "msl", "agl"}, optional
        Reference for ``profile_height``.

    Returns
    -------
    np.ndarray
        Temperature field with the same shape as ``gate_altitude``.
    """
    gate_altitude = np.asarray(gate_altitude, dtype=np.float64)
    profile_altitude = _infer_profile_altitude(profile_height, radar_altitude, height_reference)
    profile_temperature = np.asarray(profile_temperature, dtype=np.float64)
    if profile_temperature.ndim != 1:
        raise ValueError("profile_temperature must be a 1-D array")
    if profile_altitude.size != profile_temperature.size:
        raise ValueError("profile_height and profile_temperature must have the same length")
    if profile_altitude.size == 0:
        raise ValueError("profile_height cannot be empty")

    valid = np.isfinite(profile_altitude) & np.isfinite(profile_temperature)
    if np.count_nonzero(valid) < 2:
        raise ValueError("At least two finite profile points are required")

    altitude = profile_altitude[valid]
    temperature = profile_temperature[valid]
    order = np.argsort(altitude, kind="mergesort")
    altitude = altitude[order]
    temperature = temperature[order]
    flat = gate_altitude.reshape(-1)
    interpolated = np.interp(flat, altitude, temperature, left=temperature[0], right=temperature[-1])
    return interpolated.reshape(gate_altitude.shape)


def hid_beta_function(x, m, a, b):
    """
    Evaluate the fuzzy beta membership function.

    Parameters
    ----------
    x : array-like
        Input feature values.
    m, a, b : array-like
        Beta-function center, width, and slope parameters for each class.

    Returns
    -------
    np.ndarray
        Membership values with shape ``x.shape + m.shape``.
    """
    values = np.asarray(x, dtype=np.float64)
    center = np.asarray(m, dtype=np.float64)
    width = np.asarray(a, dtype=np.float64)
    slope = np.asarray(b, dtype=np.float64)
    if np.any(width == 0.0):
        raise ValueError("beta-function width cannot be zero")
    expanded = np.expand_dims(values, axis=-1)
    membership = 1.0 / (1.0 + np.power((expanded - center) / width, 2.0 * slope))
    return np.where(np.isfinite(expanded), membership, np.nan)


def classify_hydrometeors(
    dBZ=None,
    ZDR=None,
    KDP=None,
    CC=None,
    LDR=None,
    T=None,
    method="hybrid",
    band="C",
    weights=None,
    return_scores=False,
    return_confidence=False,
):
    """
    Hydrometeor classification from radar variables.

    Parameters
    ----------
    dBZ : array-like
        Reflectivity field. Required.
    ZDR, KDP, CC, LDR : array-like, optional
        Polarimetric feature fields. At least one is required.
    T : array-like, optional
        Temperature field in degrees Celsius. If omitted, the classifier uses
        a no-profile reduced-variable mode.
    method : {"hybrid", "linear"}, optional
        Weighting scheme.
    band : {"S", "C", "X"}, optional
        Radar band used to select the packaged beta parameters.
    weights : dict, optional
        Feature weights overriding :data:`DEFAULT_WEIGHTS`.
    return_scores : bool, optional
        When ``True``, also return the per-class fuzzy score cube.
    return_confidence : bool, optional
        When ``True``, also return the winning-class confidence field.

    Returns
    -------
    np.ndarray or tuple
        Class ids in ``[1, 10]`` with ``NaN`` for unclassified gates. Optional
        score and confidence outputs follow the class field.
    """
    method_key = _normalize_method(method)
    band_key = _normalize_band(band)
    if dBZ is None:
        raise ValueError("dBZ is required for hydrometeor classification")

    reflectivity = _as_float_array(dBZ, "dBZ")
    shape = reflectivity.shape
    features = {
        "dBZ": reflectivity,
        "ZDR": _as_float_array(ZDR, "ZDR", shape=shape),
        "KDP": _as_float_array(KDP, "KDP", shape=shape),
        "CC": _as_float_array(CC, "CC", shape=shape),
        "LDR": _as_float_array(LDR, "LDR", shape=shape),
        "T": _as_float_array(T, "T", shape=shape),
    }
    if not any(features[name] is not None for name in _POLAR_FEATURES):
        raise ValueError("At least one of ZDR, KDP, CC, or LDR is required")

    class_count = len(HYDROMETEOR_CLASSES)
    normalized_weights = _normalize_weights(weights)
    scores = np.zeros(shape + (class_count,), dtype=np.float64)
    denominator = np.zeros(shape + (1,), dtype=np.float64)

    if method_key == "hybrid":
        active_polar = 0
        for feature_name in _POLAR_FEATURES:
            values = features[feature_name]
            weight = normalized_weights[feature_name]
            if values is None or weight <= 0.0:
                continue
            beta = _feature_beta(feature_name, values, band_key)
            gate_weight = _gate_valid_weight(values, weight)
            scores += np.where(gate_weight[..., None] > 0.0, beta * gate_weight[..., None], 0.0)
            denominator += gate_weight[..., None]
            active_polar += 1
        if active_polar == 0:
            raise ValueError("No valid polarimetric feature remains after applying weights")
        scores = np.divide(scores, denominator, out=np.zeros_like(scores), where=denominator > 0.0)
        scores *= np.where(np.isfinite(reflectivity[..., None]), _feature_beta("dBZ", reflectivity, band_key), np.nan)
        if features["T"] is not None and normalized_weights["T"] > 0.0:
            temp_beta = _feature_beta("T", features["T"], band_key)
            scores *= np.where(np.isfinite(features["T"][..., None]), temp_beta, 1.0)
        valid_gate = np.isfinite(reflectivity) & (denominator[..., 0] > 0.0)
    else:
        active_features = 0
        for feature_name in _SUPPORTED_FEATURES:
            values = features[feature_name]
            weight = normalized_weights[feature_name]
            if values is None or weight <= 0.0:
                continue
            beta = _feature_beta(feature_name, values, band_key)
            gate_weight = _gate_valid_weight(values, weight)
            scores += np.where(gate_weight[..., None] > 0.0, beta * gate_weight[..., None], 0.0)
            denominator += gate_weight[..., None]
            active_features += 1
        if active_features == 0:
            raise ValueError("No valid feature remains after applying weights")
        scores = np.divide(scores, denominator, out=np.zeros_like(scores), where=denominator > 0.0)
        valid_gate = np.isfinite(reflectivity) & (denominator[..., 0] > 0.0)

    confidence = np.where(valid_gate, np.max(scores, axis=-1), np.nan)
    classes = np.where(valid_gate, np.argmax(scores, axis=-1) + 1, np.nan).astype(np.float32)

    outputs = [classes]
    if return_scores:
        outputs.append(scores.astype(np.float32))
    if return_confidence:
        outputs.append(confidence.astype(np.float32))
    if len(outputs) == 1:
        return outputs[0]
    return tuple(outputs)


def classify_sweep_hydrometeors(
    sweep_dataset,
    band="C",
    method="hybrid",
    weights=None,
    temperature=None,
    profile_height=None,
    profile_temperature=None,
    radar_altitude=0.0,
    height_reference="asl",
    return_scores=False,
    return_confidence=False,
    field_map=None,
):
    """
    Classify one PRD sweep dataset.

    Parameters
    ----------
    sweep_dataset : xarray.Dataset
        Sweep-level dataset from ``PRD.fields[sweep]``.
    temperature : array-like, optional
        Gate temperature field in degrees Celsius.
    profile_height, profile_temperature : array-like, optional
        1-D sounding profile. When supplied, ``temperature`` must be omitted.
    radar_altitude : float, optional
        Radar altitude in meters above sea level.
    height_reference : {"asl", "msl", "agl"}, optional
        Reference used by ``profile_height``.

    Returns
    -------
    dict
        Always contains ``hcl``. Optional entries are ``confidence``,
        ``scores``, and ``temperature`` when a profile is interpolated.
    """
    if temperature is not None and (profile_height is not None or profile_temperature is not None):
        raise ValueError("Specify either temperature or profile_height/profile_temperature, not both")
    if (profile_height is None) ^ (profile_temperature is None):
        raise ValueError("profile_height and profile_temperature must be provided together")

    mapping = dict(_DEFAULT_SWEEP_FIELDS)
    if field_map:
        mapping.update(field_map)

    reflectivity = _select_dataset_field(sweep_dataset, mapping["dBZ"], required=True)
    kwargs = {"dBZ": reflectivity.values}
    for feature_name in ("ZDR", "KDP", "CC", "LDR"):
        field = _select_dataset_field(sweep_dataset, mapping.get(feature_name, (feature_name,)))
        kwargs[feature_name] = None if field is None else field.values

    used_temperature = None
    if temperature is not None:
        used_temperature = _as_float_array(temperature, "temperature", shape=reflectivity.shape)
    elif profile_height is not None:
        used_temperature = interpolate_temperature_profile(
            np.asarray(sweep_dataset["z"].values, dtype=np.float64),
            profile_height,
            profile_temperature,
            radar_altitude=radar_altitude,
            height_reference=height_reference,
        )
    kwargs["T"] = used_temperature

    outputs = classify_hydrometeors(
        method=method,
        band=band,
        weights=weights,
        return_scores=return_scores,
        return_confidence=return_confidence,
        **kwargs
    )

    if isinstance(outputs, tuple):
        classes = outputs[0]
        cursor = 1
        scores = outputs[cursor] if return_scores else None
        cursor += 1 if return_scores else 0
        confidence = outputs[cursor] if return_confidence else None
    else:
        classes = outputs
        scores = None
        confidence = None

    result = {"hcl": classes.astype(np.float32)}
    if used_temperature is not None:
        result["temperature"] = used_temperature.astype(np.float32)
    if scores is not None:
        result["scores"] = scores.astype(np.float32)
    if confidence is not None:
        result["confidence"] = confidence.astype(np.float32)
    return result


def apply_hydrometeor_classification(
    prd,
    sweeps=None,
    inplace=False,
    band="C",
    method="hybrid",
    weights=None,
    temperature=None,
    profile_height=None,
    profile_temperature=None,
    height_reference="asl",
    output_field="HCL",
    confidence_field=None,
    temperature_field=None,
    field_map=None,
):
    """
    Add hydrometeor classification fields to selected PRD sweeps.

    Parameters
    ----------
    prd : pycwr.core.NRadar.PRD
        Source radar volume.
    sweeps : sequence of int, optional
        Sweeps to process. Defaults to all sweeps.
    inplace : bool, optional
        When ``True``, modify ``prd`` in place.
    temperature : array-like, dict, or sequence, optional
        Gate temperature field. A dict or sequence can provide per-sweep fields.
    profile_height, profile_temperature : 1-D array-like, optional
        Environmental profile used for all sweeps.
    output_field : str, optional
        Name of the class-id field to write. Defaults to ``HCL``.
    confidence_field : str, optional
        Optional field name for the max-membership confidence.
    temperature_field : str, optional
        Optional field name for the temperature field used by the classifier.

    Returns
    -------
    PRD
        The modified PRD object.
    """
    target = prd if inplace else copy.deepcopy(prd)
    if sweeps is None:
        sweeps = list(range(int(target.nsweeps)))

    radar_altitude = float(np.asarray(target.scan_info["altitude"].values))
    for sweep in sweeps:
        sweep = int(sweep)
        dataset = target.fields[sweep]
        reference = _select_dataset_field(dataset, _DEFAULT_SWEEP_FIELDS["dBZ"], required=True)
        sweep_temperature = _resolve_temperature_spec(temperature, sweep, target.nsweeps)
        results = classify_sweep_hydrometeors(
            dataset,
            band=band,
            method=method,
            weights=weights,
            temperature=sweep_temperature,
            profile_height=profile_height,
            profile_temperature=profile_temperature,
            radar_altitude=radar_altitude,
            height_reference=height_reference,
            return_confidence=confidence_field is not None,
            field_map=field_map,
        )
        dataset[output_field] = _build_data_array(
            reference,
            results["hcl"],
            output_field,
            attrs=_hcl_field_attrs(method, band, used_temperature=("temperature" in results)),
        )
        if confidence_field is not None:
            dataset[confidence_field] = _build_data_array(
                reference,
                results["confidence"],
                confidence_field,
                attrs=_confidence_field_attrs(output_field, method),
            )
        if temperature_field is not None and "temperature" in results:
            source = "environmental profile"
            if sweep_temperature is not None:
                source = "explicit temperature field"
            dataset[temperature_field] = _build_data_array(
                reference,
                results["temperature"],
                temperature_field,
                attrs=_temperature_field_attrs(source, method),
            )

    if hasattr(target, "_invalidate_cached_views"):
        target._invalidate_cached_views()
    return target


def fhc_HCL(
    dBZ=None,
    ZDR=None,
    KDP=None,
    CC=None,
    LDR=None,
    T=None,
    method="hybrid",
    band="C",
    weights=DEFAULT_WEIGHTS,
):
    """
    Backward-compatible hydrometeor classification alias.
    """
    return classify_hydrometeors(
        dBZ=dBZ,
        ZDR=ZDR,
        KDP=KDP,
        CC=CC,
        LDR=LDR,
        T=T,
        method=method,
        band=band,
        weights=weights,
    )


fhc_hcl = fhc_HCL
