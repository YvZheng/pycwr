from importlib import import_module

_LAZY_MODULES = {
    "colormap": "colormap",
    "RadarPlot": "RadarPlot",
    "SingleRadarPlot": "SingleRadarPlot",
    "SingleRadarPlotMap": "SingleRadarPlotMap",
    "VerticalSectionPlot": "VerticalSectionPlot",
    "easy": "easy",
}
_LAZY_SYMBOLS = {
    "Graph": ("RadarPlot", "Graph"),
    "GraphMap": ("RadarPlot", "GraphMap"),
    "RadarDisplay": ("RadarPlot", "RadarDisplay"),
    "RadarMapDisplay": ("RadarPlot", "RadarMapDisplay"),
    "EasyPlotResult": ("easy", "EasyPlotResult"),
    "EasyRadarPlotter": ("easy", "EasyRadarPlotter"),
    "plot": ("easy", "plot"),
    "plot_ppi": ("easy", "plot_ppi"),
    "plot_ppi_map": ("easy", "plot_ppi_map"),
    "plot_rhi": ("easy", "plot_rhi"),
    "plot_section": ("easy", "plot_section"),
    "plot_section_lonlat": ("easy", "plot_section_lonlat"),
    "plot_vvp": ("easy", "plot_vvp"),
    "plot_wind_profile": ("easy", "plot_wind_profile"),
}

__all__ = [
    *sorted(_LAZY_MODULES),
    *sorted(_LAZY_SYMBOLS),
]


def _ensure_colormap_loaded():
    if "colormap" not in globals():
        globals()["colormap"] = import_module(".colormap", __name__)
    return globals()["colormap"]


def __getattr__(name):
    module_name = _LAZY_MODULES.get(name)
    if module_name is not None:
        if module_name != "colormap":
            _ensure_colormap_loaded()
        module = import_module(f".{module_name}", __name__)
        globals()[name] = module
        return module
    symbol = _LAZY_SYMBOLS.get(name)
    if symbol is not None:
        _ensure_colormap_loaded()
        module = import_module(f".{symbol[0]}", __name__)
        value = getattr(module, symbol[1])
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
