"""Top-level package for pycwr.

This module avoids importing optional subpackages that require heavy
third-party dependencies (e.g., cartopy) so that ``import pycwr`` works
in minimal environments. Optional modules are imported lazily when
available.
"""

from . import core, io, interp, retrieve, qc

__all__ = ["core", "io", "interp", "qc", "retrieve"]

try:  # pragma: no cover - optional dependency may be missing
    from . import configure  # noqa: F401
    __all__.append("configure")
except Exception:  # ImportError or other runtime error
    configure = None  # type: ignore

try:  # pragma: no cover - optional dependency may be missing
    from . import draw  # noqa: F401
    __all__.append("draw")
except Exception:
    draw = None  # type: ignore
