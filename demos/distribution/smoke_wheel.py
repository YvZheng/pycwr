#!/usr/bin/env python3
"""Check the installed binary backend, interpolation, plotting, and NetCDF IO."""
from __future__ import annotations

import argparse
import importlib
import importlib.machinery
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
import tempfile


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="pycwr-wheel-smoke-") as temp:
        temp_path = Path(temp)
        os.environ.setdefault("MPLCONFIGDIR", str(temp_path / "matplotlib"))
        os.environ.setdefault("XDG_CACHE_HOME", str(temp_path / "cache"))
        import numpy as np
        import xarray as xr
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pycwr
        from pycwr import core
        from pycwr.core import RadarGrid
        from pycwr.io import read_auto
        from pycwr.draw import plot_ppi
        from pycwr.qc import apply_dualpol_qc

        assert sys.version_info[:2] == (3, 11), sys.version
        assert importlib.metadata.version("pycwr") == "1.0.8"
        backend = importlib.import_module("pycwr.core.RadarGridC")
        backend_path = str(backend.__file__)
        assert any(backend_path.endswith(ext) for ext in importlib.machinery.EXTENSION_SUFFIXES), backend_path
        assert core.RadarGridC is backend, "pycwr selected its Python fallback instead of the binary backend"
        assert Path(pycwr.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()), pycwr.__file__
        assert all(callable(func) for func in (read_auto, plot_ppi, apply_dualpol_qc))

        # A known bilinear interpolation result catches callable-but-broken extensions.
        np.testing.assert_allclose(backend.interp_ppi(5, 1500, 0, 10, 1000, 2000, 0, 10, 20, 30, -999), 15)
        az = np.arange(0, 360, 5, dtype=np.float64)
        ranges = np.arange(250, 5250, 250, dtype=np.float64)
        values = np.ascontiguousarray(20 + 2 * np.cos(np.deg2rad(az[:, None])) + ranges[None, :] / 1000)
        grid_x, grid_y = np.meshgrid(np.linspace(-2000, 2000, 9), np.linspace(-2000, 2000, 9), indexing="ij")
        call_args = (az, ranges, 0.5, values, 100.0, grid_x, grid_y, -999.0)
        actual = backend.ppi_to_grid(*call_args)
        expected = RadarGrid.ppi_to_grid(*call_args)
        np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-7)
        assert np.count_nonzero(actual != -999.0) >= 72
        assert actual[4, 4] == -999.0, "The unsampled radar origin must remain masked"

        dataset = xr.Dataset({"reflectivity": (("x", "y"), actual)})
        dataset.to_netcdf(temp_path / "grid.nc", engine="netcdf4")
        with xr.open_dataset(temp_path / "grid.nc", engine="netcdf4") as restored:
            np.testing.assert_array_equal(restored.reflectivity.values, actual)
        figure, axis = plt.subplots()
        axis.pcolormesh(grid_x, grid_y, np.ma.masked_equal(actual, -999.0), shading="auto")
        figure.savefig(temp_path / "grid.png")
        plt.close(figure)
        assert (temp_path / "grid.png").stat().st_size > 1000
        result = {"status": "passed", "python": sys.version, "platform": platform.platform(),
                  "machine": platform.machine(), "compiled_extension": backend_path,
                  "numpy": np.__version__, "pycwr": pycwr.__version__,
                  "checks": ["installed binary extension", "no runtime fallback", "known interpolation",
                             "C/Python grid agreement", "blind-zone mask", "NetCDF roundtrip", "PNG output"]}
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
