"""Top-level package for pycwr.

This module avoids importing optional or heavy subpackages at import time
so that ``import pycwr`` remains lightweight and warning-free.
"""

from importlib import import_module

__version__ = "1.0.3"

_LAZY_MODULES = {"core", "io", "interp", "qc", "retrieve", "configure", "draw"}

__all__ = ["__version__", *sorted(_LAZY_MODULES)]

def __getattr__(name):
    if name in _LAZY_MODULES:
        module = import_module(f".{name}", __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
