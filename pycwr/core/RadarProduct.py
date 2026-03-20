import numpy as np


PRODUCT_REFERENCE_NOTES = {
    "VIL": [
        "Greene and Clark (1972), Monthly Weather Review, doi:10.1175/1520-0493(1972)100<0548:VILWNA>2.3.CO;2",
        "NOAA WDTD VIL product guidance: https://vlab.noaa.gov/web/wdtd/-/vertically-integrated-liquid-vil-",
    ],
    "ET": [
        "Lakshmanan et al. (2013), Weather and Forecasting, doi:10.1175/WAF-D-12-00084.1",
        "NOAA WDTD xx dBZ Echo Top guidance: https://vlab.noaa.gov/web/wdtd/-/xx-dbz-echo-top-et-",
        "WSR-88D ROC ICD 2620003Y (Enhanced Echo Tops): https://www.roc.noaa.gov/public-documents/icds/2620003Y.pdf",
    ],
    "MOSAIC": [
        "吴翀;双偏振雷达的资料质量分析,相态识別及组网应用[D];南京信息工程大学;2018年。相关组网讨论见第4章。",
        "Lakshmanan et al. (2006), Weather and Forecasting, doi:10.1175/WAF942.1",
        "NOAA WDTD Radar Quality Index: https://vlab.noaa.gov/web/wdtd/-/radar-quality-index-rqi-",
    ],
}


def dBZ_to_linear(dbz):
    return np.power(10.0, np.asarray(dbz, dtype=np.float64) / 10.0)


def derive_cr(volume, fillvalue=-999.0):
    volume = np.asarray(volume, dtype=np.float64)
    valid = np.isfinite(volume) & (volume != fillvalue)
    reduced = np.max(np.where(valid, volume, -np.inf), axis=0)
    return np.where(np.any(valid, axis=0), reduced, np.nan)


def derive_vil(volume, level_heights, fillvalue=-999.0, min_dbz=18.0, max_dbz_cap=56.0):
    volume = np.asarray(volume, dtype=np.float64)
    level_heights = np.asarray(level_heights, dtype=np.float64)
    valid = np.isfinite(volume) & (volume != fillvalue) & (volume >= float(min_dbz))
    if level_heights.size < 2:
        return np.full(volume.shape[1:], np.nan, dtype=np.float64)
    clipped = np.clip(volume, float(min_dbz), float(max_dbz_cap))
    liquid = 3.44e-6 * np.power(dBZ_to_linear(clipped), 4.0 / 7.0)
    liquid = np.where(valid, liquid, 0.0)
    if hasattr(np, "trapezoid"):
        vil = np.trapezoid(liquid, level_heights, axis=0)
    else:  # pragma: no cover
        vil = np.trapz(liquid, level_heights, axis=0)
    return np.where(np.any(valid, axis=0), vil, np.nan)


def derive_et(volume, level_heights, fillvalue=-999.0, threshold_dbz=18.0, return_topped=False):
    volume = np.asarray(volume, dtype=np.float64)
    level_heights = np.asarray(level_heights, dtype=np.float64)
    if level_heights.size == 0:
        et = np.full(volume.shape[1:], np.nan, dtype=np.float64)
        topped = np.zeros(volume.shape[1:], dtype=np.uint8)
        return (et, topped) if return_topped else et
    valid = np.isfinite(volume) & (volume != fillvalue)
    above = valid & (volume >= float(threshold_dbz))
    has_echo = np.any(above, axis=0)
    et = np.full(volume.shape[1:], np.nan, dtype=np.float64)
    topped = np.zeros(volume.shape[1:], dtype=np.uint8)
    if not np.any(has_echo):
        return (et, topped) if return_topped else et

    top_index = volume.shape[0] - 1 - np.argmax(above[::-1], axis=0)
    flat_et = et.ravel()
    flat_topped = topped.ravel()
    flat_top = top_index.ravel()
    flat_has = has_echo.ravel()
    flat_volume = volume.reshape(volume.shape[0], -1)
    for i in np.where(flat_has)[0]:
        top = int(flat_top[i])
        if top >= level_heights.size - 1:
            flat_et[i] = float(level_heights[top])
            flat_topped[i] = 1
            continue
        z0 = float(level_heights[top])
        z1 = float(level_heights[top + 1])
        v0 = float(flat_volume[top, i])
        v1 = float(flat_volume[top + 1, i])
        if np.isfinite(v1) and v1 < threshold_dbz and v0 > threshold_dbz and v1 != v0:
            weight = (float(threshold_dbz) - v0) / (v1 - v0)
            flat_et[i] = z0 + weight * (z1 - z0)
        else:
            flat_et[i] = z0
    return (et, topped) if return_topped else et
