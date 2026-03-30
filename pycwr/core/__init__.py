import importlib
import sys
import warnings

from . import NRadar, PyartRadar, transforms, RadarGrid, interop, RadarProduct

_RADARGRIDC_ABI_WARNING_FRAGMENT = "numpy.ndarray size changed"


def _load_radargrid_backend(import_module=None):
    """Load the optional C backend, falling back safely on ABI issues."""
    importer = import_module or importlib.import_module
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            module = importer(".RadarGridC", __name__)
        except ImportError:
            return RadarGrid
        except Exception as exc:  # pragma: no cover - defensive fallback for platform-specific import failures.
            warnings.warn(
                "Falling back to pycwr.core.RadarGrid because RadarGridC could not be imported: %s"
                % exc,
                RuntimeWarning,
                stacklevel=2,
            )
            return RadarGrid

    abi_warning = next(
        (
            warning_message
            for warning_message in caught
            if _RADARGRIDC_ABI_WARNING_FRAGMENT in str(warning_message.message)
        ),
        None,
    )
    if abi_warning is not None:
        sys.modules.pop(f"{__name__}.RadarGridC", None)
        warnings.warn(
            "Falling back to pycwr.core.RadarGrid because RadarGridC appears incompatible with the current NumPy/Python ABI.",
            RuntimeWarning,
            stacklevel=2,
        )
        return RadarGrid

    for warning_message in caught:
        warnings.warn(
            str(warning_message.message),
            warning_message.category,
            stacklevel=2,
        )
    return module


RadarGridC = _load_radargrid_backend()

__all__ = ["NRadar", "PyartRadar", "transforms", "RadarGrid", "RadarGridC", "RadarProduct", "interop"]
