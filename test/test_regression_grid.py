"""Grid-product consistency checks that also serve as low-level examples."""

import unittest
import warnings

import numpy as np


class RadarGridConsistencyTests(unittest.TestCase):
    def test_load_radargrid_backend_falls_back_on_numpy_abi_warning(self):
        from pycwr import core

        def fake_import(_name, _package=None):
            warnings.warn(
                "numpy.ndarray size changed, may indicate binary incompatibility. Expected 16 from C header, got 96 from PyObject",
                RuntimeWarning,
            )
            return object()

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            backend = core._load_radargrid_backend(import_module=fake_import)

        self.assertIs(backend, core.RadarGrid)
        self.assertTrue(any("Falling back to pycwr.core.RadarGrid" in str(item.message) for item in caught))

    @staticmethod
    def _synthetic_volume(sweep_values):
        azimuth = np.array([0.0, 90.0, 180.0, 270.0], dtype=np.float64)
        ranges = np.array([1000.0, 2000.0], dtype=np.float64)
        vol_azimuth = [azimuth.copy(), azimuth.copy()]
        vol_range = [ranges.copy(), ranges.copy()]
        fix_elevation = np.array([1.0, 2.0], dtype=np.float64)
        vol_value = [
            np.full((azimuth.size, ranges.size), sweep_values[0], dtype=np.float64),
            np.full((azimuth.size, ranges.size), sweep_values[1], dtype=np.float64),
        ]
        return vol_azimuth, vol_range, fix_elevation, vol_value

    def _target_grid(self, radar_height):
        from pycwr.core.transforms import antenna_to_cartesian_cwr

        x, y, z = antenna_to_cartesian_cwr(
            np.array([1500.0], dtype=np.float64),
            np.array([0.0], dtype=np.float64),
            np.array([1.5], dtype=np.float64),
            radar_height,
        )
        return (
            np.array([[x[0]]], dtype=np.float64),
            np.array([[y[0]]], dtype=np.float64),
            np.array([z[0]], dtype=np.float64),
        )

    def test_get_cappi_3d_python_and_cython_match(self):
        from pycwr.core import RadarGrid

        try:
            from pycwr.core import RadarGridC
        except ImportError:
            self.skipTest("RadarGridC extension is not available")

        radar_height = 100.0
        vol_azimuth, vol_range, fix_elevation, vol_value = self._synthetic_volume((10.0, 20.0))
        grid_x, grid_y, levels = self._target_grid(radar_height)

        py_value = RadarGrid.get_CAPPI_3d(
            vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, grid_x, grid_y, levels, -999.0
        )
        cy_value = RadarGridC.get_CAPPI_3d(
            vol_azimuth, vol_range, fix_elevation, vol_value, radar_height, grid_x, grid_y, levels, -999.0
        )

        self.assertTrue(np.allclose(py_value, cy_value, rtol=1e-10, atol=1e-6))
        self.assertTrue(np.isclose(py_value[0, 0, 0], 15.0, rtol=0.0, atol=1e-6))

    def test_get_mosaic_cappi_3d_python_and_cython_match(self):
        from pycwr.core import RadarGrid

        try:
            from pycwr.core import RadarGridC
        except ImportError:
            self.skipTest("RadarGridC extension is not available")

        radar_height = np.array([100.0, 100.0], dtype=np.float64)
        radar_x = np.array([0.0, 0.0], dtype=np.float64)
        radar_y = np.array([0.0, 0.0], dtype=np.float64)
        grid_x, grid_y, levels = self._target_grid(radar_height[0])

        vol_azimuth = []
        vol_range = []
        fix_elevation = []
        vol_value = []
        for sweep_values in ((10.0, 20.0), (30.0, 40.0)):
            radar_vol_azimuth, radar_vol_range, radar_fix_elevation, radar_vol_value = self._synthetic_volume(sweep_values)
            vol_azimuth.append(radar_vol_azimuth)
            vol_range.append(radar_vol_range)
            fix_elevation.append(radar_fix_elevation)
            vol_value.append(radar_vol_value)

        py_value = RadarGrid.get_mosaic_CAPPI_3d(
            vol_azimuth, vol_range, fix_elevation, vol_value, radar_x, radar_y, radar_height, grid_x, grid_y, levels, -999.0
        )
        cy_value = RadarGridC.get_mosaic_CAPPI_3d(
            vol_azimuth, vol_range, fix_elevation, vol_value, radar_x, radar_y, radar_height, grid_x, grid_y, levels, -999.0
        )

        self.assertTrue(np.allclose(py_value, cy_value, rtol=1e-10, atol=1e-6))
        self.assertTrue(np.isclose(py_value[0, 0, 0], 35.0, rtol=0.0, atol=1e-6))

    def test_get_cappi_xy_single_sweep_python_and_cython_match(self):
        from pycwr.core import RadarGrid
        from pycwr.core.transforms import antenna_to_cartesian_cwr

        try:
            from pycwr.core import RadarGridC
        except ImportError:
            self.skipTest("RadarGridC extension is not available")

        radar_height = 100.0
        azimuth = np.array([0.0, 90.0, 180.0, 270.0], dtype=np.float64)
        ranges = np.array([1000.0, 2000.0], dtype=np.float64)
        fix_elevation = np.array([1.0], dtype=np.float64)
        vol_azimuth = [azimuth.copy()]
        vol_range = [ranges.copy()]
        vol_value = [np.full((azimuth.size, ranges.size), 10.0, dtype=np.float64)]
        x, y, _ = antenna_to_cartesian_cwr(
            np.array([1500.0], dtype=np.float64),
            np.array([0.0], dtype=np.float64),
            np.array([1.0], dtype=np.float64),
            radar_height,
        )
        grid_x = np.array([[x[0]]], dtype=np.float64)
        grid_y = np.array([[y[0]]], dtype=np.float64)

        py_value = RadarGrid.get_CAPPI_xy(
            vol_azimuth,
            vol_range,
            fix_elevation,
            vol_value,
            radar_height,
            grid_x,
            grid_y,
            126.0,
            -999.0,
        )
        cy_value = RadarGridC.get_CAPPI_xy(
            vol_azimuth,
            vol_range,
            fix_elevation,
            vol_value,
            radar_height,
            grid_x,
            grid_y,
            126.0,
            -999.0,
        )

        self.assertTrue(np.allclose(py_value, cy_value, rtol=1e-10, atol=1e-6))
        self.assertTrue(np.isclose(py_value[0, 0], 10.0, rtol=0.0, atol=1e-6))


if __name__ == "__main__":
    unittest.main()
