# -*- coding: utf-8 -*-
"""Web viewer color palette definitions based on the provided XML scheme."""

from dataclasses import dataclass

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, PowerNorm

from ..configure.default_config import HYDROMETEOR_CLASS_COLORS, HYDROMETEOR_CLASS_SHORT_LABELS_ZH
from ..draw._plot_core import ColorbarOptions, PlotStyle, default_colorbar_label, resolve_field_style


SPECIAL_COLORS = {
    "Text1": "#ffffff",
    "Text2": "#ffffff",
    "Axis": "#c0c0c0",
    "Ring": "#ffa07a",
    "Symbol": "#ffffff",
    "SymbolUT": "#0000ff",
    "Background1": "#000000",
    "Background2": "#000000",
    "Wind": "#ffffff",
    "WindUT": "#0000ff",
    "RF": "#7c007c",
    "ND": "#7c007c",
}


@dataclass(frozen=True)
class FieldPalette:
    levels: tuple
    labels: tuple
    colors: tuple

    def boundaries(self):
        levels = np.asarray(self.levels, dtype=float)
        if levels.size == 1:
            return np.array([levels[0] - 0.5, levels[0] + 0.5], dtype=float)
        step = levels[-1] - levels[-2]
        if step == 0:
            step = 1.0
        return np.concatenate([levels, [levels[-1] + step]])

    def discrete_cmap(self, field_name):
        return ListedColormap(list(self.colors), name="pycwr_web_%s" % field_name)

    def continuous_cmap(self, field_name):
        return LinearSegmentedColormap.from_list("pycwr_web_linear_%s" % field_name, list(self.colors))


PALETTES = {
    "dBZ": FieldPalette(
        levels=(-5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65),
        labels=("-5", "0", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", "65"),
        colors=(
            "#00aca4", "#c0c0fe", "#7a72ee", "#1e26d0", "#a6fca8",
            "#00ea00", "#10921a", "#fcf464", "#c8c802", "#8c8c00",
            "#feacac", "#fe6454", "#ee0230", "#d48efe", "#aa24fa",
        ),
    ),
    "dBT": FieldPalette(
        levels=(-5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65),
        labels=("-5", "0", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", "65"),
        colors=(
            "#00aca4", "#c0c0fe", "#7a72ee", "#1e26d0", "#a6fca8",
            "#00ea00", "#10921a", "#fcf464", "#c8c802", "#8c8c00",
            "#feacac", "#fe6454", "#ee0230", "#d48efe", "#aa24fa",
        ),
    ),
    "V": FieldPalette(
        levels=(-30, -27, -20, -15, -10, -5, -1, 0, 1, 5, 10, 15, 20, 27),
        labels=("-30", "-27", "-20", "-15", "-10", "-5", "-1", "0", "1", "5", "10", "15", "20", "27"),
        colors=(
            "#7ee0fe", "#00e0fe", "#00b0b0", "#00fe00", "#00c400", "#008000",
            "#fefefe", "#fcfcfc", "#fe0000", "#fe5858", "#feb0b0", "#fe7c00",
            "#fed200", "#fefe00",
        ),
    ),
    "W": FieldPalette(
        levels=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13),
        labels=("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"),
        colors=(
            "#e0e0e0", "#7ce0e0", "#00e0e0", "#00b0b0", "#00fefe", "#00c400",
            "#008000", "#fefe00", "#fed200", "#fe7c00", "#feb0b0", "#af5858",
            "#fe0000", "#e60000",
        ),
    ),
    "SNRH": FieldPalette(
        levels=(-5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65),
        labels=("-5", "0", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", "65"),
        colors=(
            "#00aca4", "#c0c0fe", "#7a72ee", "#1e26d0", "#a6fca8",
            "#00ea00", "#10921a", "#fcf464", "#c8c802", "#8c8c00",
            "#feacac", "#fe6454", "#ee0230", "#d48efe", "#aa24fa",
        ),
    ),
    "SNRV": FieldPalette(
        levels=(-5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65),
        labels=("-5", "0", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", "65"),
        colors=(
            "#00aca4", "#c0c0fe", "#7a72ee", "#1e26d0", "#a6fca8",
            "#00ea00", "#10921a", "#fcf464", "#c8c802", "#8c8c00",
            "#feacac", "#fe6454", "#ee0230", "#d48efe", "#aa24fa",
        ),
    ),
    "SQI": FieldPalette(
        levels=(0, 0.06, 0.13, 0.19, 0.25, 0.31, 0.38, 0.44, 0.5, 0.56, 0.62, 0.69, 0.75, 0.81, 0.87, 0.94),
        labels=("0", "0.06", "0.13", "0.19", "0.25", "0.31", "0.38", "0.44", "0.5", "0.56", "0.62", "0.69", "0.75", "0.81", "0.87", "0.94"),
        colors=(
            "#003cff", "#00efef", "#00babf", "#00837d", "#008938", "#00b729",
            "#00da0d", "#00ff00", "#ffff3b", "#fff000", "#ffc600", "#ffa500",
            "#ff7200", "#ff1f00", "#c10000", "#d400aa",
        ),
    ),
    "RR": FieldPalette(
        levels=(0.05, 0.1, 0.2, 0.5, 0.8, 1, 2, 5, 8, 10, 20, 50, 80, 100, 200, 500),
        labels=("0.05", "0.1", "0.2", "0.5", "0.8", "1", "2", "5", "8", "10", "20", "50", "80", "100", "200", "500"),
        colors=(
            "#003cff", "#00efef", "#00babf", "#00837d", "#008938", "#00b729",
            "#00da0d", "#00ff00", "#ffff3b", "#fff000", "#ffc600", "#ffa500",
            "#ff7200", "#ff1f00", "#c10000", "#d400aa",
        ),
    ),
    "ZDR": FieldPalette(
        levels=(-4, -3, -2, -1, 0, 0.2, 0.5, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5),
        labels=("-4", "-3", "-2", "-1", "0", "0.2", "0.5", "0.8", "1", "1.5", "2", "2.5", "3", "3.5", "4", "5"),
        colors=(
            "#464646", "#6e6e6e", "#969696", "#c8c8c8", "#dcf0dc", "#00c027",
            "#00e80a", "#24ff24", "#ffff1e", "#ffe600", "#ffbc00", "#ff9800",
            "#ff5e00", "#f20f00", "#bb003a", "#ff00ff",
        ),
    ),
    "LDR": FieldPalette(
        levels=(-30, -28, -26, -24, -22, -20, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0),
        labels=("-30", "-28", "-26", "-24", "-22", "-20", "-18", "-16", "-14", "-12", "-10", "-8", "-6", "-4", "-2", "0"),
        colors=(
            "#003cff", "#00efef", "#00babf", "#00837d", "#008938", "#00b729",
            "#00da0d", "#00ff00", "#ffff3b", "#fff000", "#ffc600", "#ffa500",
            "#ff7200", "#ff1f00", "#c10000", "#d400aa",
        ),
    ),
    "CC": FieldPalette(
        levels=(0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.01),
        labels=("0", "0.1", "0.3", "0.5", "0.6", "0.7", "0.8", "0.85", "0.9", "0.92", "0.94", "0.95", "0.96", "0.97", "0.98", "0.99", "1.01"),
        colors=(
            "#003cff", "#00efef", "#00babf", "#00837d", "#008938", "#00b729",
            "#00da0d", "#00ff00", "#ffff3b", "#fff000", "#ffc600", "#ffa500",
            "#ff7200", "#ff1f00", "#c10000", "#d400aa", "#7c007c",
        ),
    ),
    "PhiDP": FieldPalette(
        levels=(0, 24, 48, 72, 96, 120, 144, 168, 192, 216, 240, 264, 288, 312, 336, 360),
        labels=("0", "24", "48", "72", "96", "120", "144", "168", "192", "216", "240", "264", "288", "312", "336", "360"),
        colors=(
            "#003cff", "#00efef", "#00babf", "#00837d", "#008938", "#00b729",
            "#00da0d", "#00ff00", "#ffff3b", "#fff000", "#ffc600", "#ffa500",
            "#ff7200", "#ff1f00", "#c10000", "#d400aa",
        ),
    ),
    "KDP": FieldPalette(
        levels=(-0.8, -0.4, -0.2, -0.1, 0.1, 0.15, 0.22, 0.33, 0.5, 0.75, 1.1, 1.7, 2.4, 3.1, 7, 20),
        labels=("-0.8", "-0.4", "-0.2", "-0.1", "0.1", "0.15", "0.22", "0.33", "0.5", "0.75", "1.1", "1.7", "2.4", "3.1", "7", "20"),
        colors=(
            "#00ffff", "#00efef", "#00a8ac", "#b4b4b4", "#b4b4b4", "#00c027",
            "#00e80a", "#24ff24", "#ffff1e", "#ffe600", "#ffbc00", "#ff9800",
            "#ff5e00", "#f20f00", "#bb003a", "#ff00ff",
        ),
    ),
    "HCL": FieldPalette(
        levels=tuple(range(1, len(HYDROMETEOR_CLASS_SHORT_LABELS_ZH) + 1)),
        labels=HYDROMETEOR_CLASS_SHORT_LABELS_ZH,
        colors=HYDROMETEOR_CLASS_COLORS,
    ),
}


def get_palette(field_name):
    return PALETTES.get(field_name)


def _finite_values(field_data):
    values = np.asarray(getattr(field_data, "values", field_data), dtype=float)
    finite = values[np.isfinite(values)]
    return finite


def _bounded_range(raw_range, lower, upper):
    vmin = max(float(raw_range[0]), float(lower))
    vmax = min(float(raw_range[1]), float(upper))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
        return float(lower), float(upper)
    return vmin, vmax


def _percentile_range(field_data, low, high, fallback):
    finite = _finite_values(field_data)
    if finite.size == 0:
        return fallback
    vmin = float(np.nanpercentile(finite, low))
    vmax = float(np.nanpercentile(finite, high))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
        return fallback
    return vmin, vmax


def _preset_range(field_key, field_data, palette, preset):
    palette_range = None
    if palette is not None:
        palette_range = (float(palette.levels[0]), float(palette.levels[-1]))
    if preset == "weak_echo":
        if field_key in ("dBZ", "dBT", "SNRH", "SNRV"):
            return _bounded_range((-5.0, 30.0), *(palette_range or (-5.0, 30.0)))
        if field_key == "RR":
            return _bounded_range((0.05, 50.0), *(palette_range or (0.05, 50.0)))
    detail_range = _percentile_range(field_data, 1.0, 99.0, palette_range or (0.0, 1.0))
    if palette_range is None:
        return detail_range
    return _bounded_range(detail_range, *palette_range)


def build_web_style(
    radar,
    sweep,
    field_key,
    field_data,
    *,
    continuous=False,
    preset="standard",
    clabel=None,
    colorbar_visible=True,
    colorbar_orientation="vertical",
):
    palette = get_palette(field_key)
    preset = (preset or "standard").strip().lower()
    if preset not in ("standard", "fine_detail", "weak_echo"):
        preset = "standard"
    if preset != "standard":
        value_range = _preset_range(field_key, field_data, palette, preset)
        gamma = 0.82 if preset == "fine_detail" else 0.72
        colorbar = ColorbarOptions(
            visible=colorbar_visible,
            orientation=colorbar_orientation,
            label=clabel or default_colorbar_label(field_key),
            ticks=list(palette.levels) if palette is not None else None,
            ticklabels=list(palette.labels) if palette is not None else None,
        )
        if palette is not None:
            return PlotStyle(
                cmap=palette.continuous_cmap(field_key),
                value_range=value_range,
                norm=PowerNorm(gamma=gamma, vmin=value_range[0], vmax=value_range[1], clip=True),
                continuous=True,
                colorbar=colorbar,
            )
        fallback = resolve_field_style(
            radar,
            sweep,
            field_key,
            field_data,
            value_range=value_range,
            continuous=True,
            colorbar_visible=colorbar_visible,
            colorbar_orientation=colorbar_orientation,
            colorbar_label=clabel,
        )
        fallback.norm = PowerNorm(gamma=gamma, vmin=value_range[0], vmax=value_range[1], clip=True)
        return fallback
    if palette is None:
        return resolve_field_style(
            radar,
            sweep,
            field_key,
            field_data,
            continuous=continuous,
            colorbar_visible=colorbar_visible,
            colorbar_orientation=colorbar_orientation,
            colorbar_label=clabel,
        )
    if continuous:
        return PlotStyle(
            cmap=palette.continuous_cmap(field_key),
            value_range=(float(palette.levels[0]), float(palette.levels[-1])),
            continuous=True,
            colorbar=ColorbarOptions(
                visible=colorbar_visible,
                orientation=colorbar_orientation,
                label=clabel or default_colorbar_label(field_key),
                ticks=list(palette.levels),
                ticklabels=list(palette.labels),
            ),
        )
    return PlotStyle(
        cmap=palette.discrete_cmap(field_key),
        levels=palette.boundaries(),
        continuous=False,
        colorbar=ColorbarOptions(
            visible=colorbar_visible,
            orientation=colorbar_orientation,
            label=clabel or default_colorbar_label(field_key),
            ticks=list(palette.levels),
            ticklabels=list(palette.labels),
        ),
    )
