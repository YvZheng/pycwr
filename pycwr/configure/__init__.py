from importlib import import_module

_LAZY_MODULES = {
    "default_config",
    "location_config",
    "pyart_config",
    "pyart_default_config",
    "pyart_lazydict",
}

__all__ = sorted(_LAZY_MODULES)


def __getattr__(name):
    if name in _LAZY_MODULES:
        module = import_module(f".{name}", __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
