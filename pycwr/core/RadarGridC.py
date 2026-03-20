"""
Pure-Python mirror of the ``RadarGridC`` Cython module.

When the compiled extension is unavailable, importing ``pycwr.core.RadarGridC``
will fall back to this source module and preserve the same public API.
"""

from .RadarGrid import (
    get_CAPPI_3d,
    get_CAPPI_xy,
    get_CR_xy,
    get_mosaic_CAPPI_3d,
    interp_azimuth,
    interp_ppi,
    ppi_to_grid,
)
from .transforms import (
    _azimuth as xy_to_azimuth,
    antenna_to_cartesian_cwr as antenna_to_cartesian,
    cartesian_xy_elevation_to_range_z as xye_to_antenna,
    cartesian_xyz_to_antenna as cartesian_to_antenna,
)


__all__ = [
    "antenna_to_cartesian",
    "xye_to_antenna",
    "cartesian_to_antenna",
    "xy_to_azimuth",
    "interp_ppi",
    "ppi_to_grid",
    "get_CAPPI_xy",
    "get_CAPPI_3d",
    "get_mosaic_CAPPI_3d",
    "interp_azimuth",
    "get_CR_xy",
]
