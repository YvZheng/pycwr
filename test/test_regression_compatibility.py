import unittest
import importlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


class CompatibilitySmokeTests(unittest.TestCase):
    @staticmethod
    def _sample_path():
        return (
            Path(__file__).resolve().parents[2]
            / "Z9046"
            / "Z_RADR_I_Z9046_20260317065928_O_DOR_SAD_CAP_FMT.bin.bz2"
        )

    def test_import_top_level_package(self):
        import pycwr

        self.assertIn("core", pycwr.__all__)
        self.assertIn("draw", pycwr.__all__)

    def test_import_config_package(self):
        import pycwr.configure as configure

        self.assertIn("location_config", configure.__all__)
        self.assertIn("pyart_config", configure.__all__)
        self.assertTrue(hasattr(configure.location_config, "radar_info"))
        self.assertGreater(len(configure.location_config.radar_info.index), 0)

    def test_draw_package_exports_lazy_symbols(self):
        import pycwr.draw as draw

        self.assertIn("Graph", draw.__all__)
        self.assertEqual(draw.Graph.__name__, "Graph")

    def test_optional_web_viewer_package_does_not_break_base_import(self):
        import pycwr.GraphicalInterface as gui

        self.assertTrue(hasattr(gui, "__all__"))
        self.assertIn("create_app", gui.__all__)
        self.assertEqual(sorted(gui.__all__), ["create_app", "launch", "web_app", "web_colors"])

    def test_plot_xy_smoke(self):
        from pycwr.draw.RadarPlot import plot_xy

        x = np.array([[0.0, 1000.0], [0.0, 1000.0]])
        y = np.array([[0.0, 0.0], [1000.0, 1000.0]])
        data = np.array([[1.0, 2.0], [3.0, 4.0]])

        fig, ax = plt.subplots()
        try:
            mesh = plot_xy(ax, x, y, data, cbar=False)
            self.assertIsNotNone(mesh)
        finally:
            plt.close(fig)

    def test_sample_read_regression(self):
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        self.assertEqual(str(prd.scan_info.scan_type.values), "ppi")
        self.assertEqual(int(prd.nsweeps), 9)
        self.assertEqual(int(prd.nrays), 3264)
        self.assertEqual(prd.fields[0]["dBZ"].shape, (361, 920))
        self.assertAlmostEqual(float(prd.fields[0]["dBZ"].min(skipna=True)), -22.5)
        self.assertAlmostEqual(float(prd.fields[0]["dBZ"].max(skipna=True)), 65.5)

    def test_extended_reflectivity_sidecar(self):
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        self.assertIn("dBZ", prd.extended_fields)
        self.assertTrue(prd.has_extended_field(0, "dBZ"))
        self.assertTrue(prd.has_extended_field(1, "dBZ"))
        self.assertFalse(prd.has_extended_field(2, "dBZ"))

        default_field = prd.get_sweep_field(0, "dBZ")
        native = prd.get_native_sweep_field(0, "dBZ")
        aligned = prd.fields[0]["dBZ"]
        self.assertEqual(default_field.shape, native.shape)
        self.assertEqual(native.shape[0], aligned.shape[0])
        self.assertGreater(native.shape[1], aligned.shape[1])
        self.assertEqual(native.dtype, np.float32)
        self.assertGreater(float(native.range.values[-1]), float(prd.fields[0].range.values[-1]))

        prd.ordered_az(inplace=True)
        native_sorted = prd.get_native_sweep_field(0, "dBZ")
        self.assertEqual(native_sorted.dims[0], "azimuth")
        self.assertEqual(native_sorted.shape, native.shape)

    def test_wsr98d_special_resolution_codec(self):
        from pycwr.io.WSR98DFile import _decode_wsr98d_resolution, _encode_wsr98d_resolution

        self.assertEqual(_decode_wsr98d_resolution(125), 125.0)
        self.assertEqual(_decode_wsr98d_resolution(39018), 62.5)
        self.assertEqual(_encode_wsr98d_resolution(125.0), 125)
        self.assertEqual(_encode_wsr98d_resolution(62.5), 39018)

    def test_native_range_mode_products_and_plot(self):
        from pycwr.draw.RadarPlot import Graph
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        x = np.arange(-50_000.0, 50_001.0, 5_000.0)
        y = np.arange(-50_000.0, 50_001.0, 5_000.0)
        product_levels = np.array([1000.0, 2000.0, 3000.0], dtype=np.float64)
        prd.add_product_CR_xy(x, y)
        prd.add_product_CAPPI_xy(x, y, 3000.0)
        prd.add_product_VIL_xy(x, y, product_levels)
        prd.add_product_ET_xy(x, y, product_levels)
        self.assertIn("CR_native", prd.product)
        self.assertIn("CAPPI_3000_native", prd.product)
        self.assertIn("VIL_native", prd.product)
        self.assertIn("ET_native", prd.product)
        self.assertIn("ET_TOPPED_native", prd.product)

        fig, ax = plt.subplots()
        try:
            graph = Graph(prd)
            mesh = graph.plot_ppi(ax, 0, "dBZ", range_mode="native", cbar=False)
            self.assertIsNotNone(mesh)
        finally:
            plt.close(fig)

    def test_native_range_mode_sections_and_network(self):
        from pycwr.core.NRadar import grid_3d_network_xy
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        native = prd.get_native_sweep_field(0, "dBZ")
        aligned = prd.fields[0]["dBZ"]

        rhi_range, rhi_height, rhi_field = prd.get_RHI_data(
            float(native.azimuth.values[0]),
            field_name="dBZ",
            range_mode="native",
        )
        self.assertEqual(rhi_field[0].shape[1], native.shape[1])
        self.assertGreater(rhi_field[0].shape[1], aligned.shape[1])
        self.assertEqual(rhi_range[0].shape, rhi_height[0].shape)

        vcs_range, vcs_height, vcs_field = prd.get_vcs_data(
            (-10_000.0, 0.0),
            (10_000.0, 0.0),
            "dBZ",
            range_mode="native",
        )
        self.assertEqual(len(vcs_range), len(vcs_field))
        self.assertGreaterEqual(len(vcs_range), 1)

        x = np.arange(-20_000.0, 20_001.0, 10_000.0)
        y = np.arange(-20_000.0, 20_001.0, 10_000.0)
        network = grid_3d_network_xy(
            [prd],
            x,
            y,
            np.array([3000.0]),
            range_mode="native",
        )
        self.assertIn("network_3d", network)
        self.assertEqual(network["network_3d"].attrs["range_mode"], "native")

    def test_prd_pyart_export_remains_backward_compatible(self):
        from pycwr.core.PyartRadar import Radar as InternalRadar
        from pycwr.io import read_auto
        import pyart

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        self.assertIsNone(prd.PyartRadar)
        radar = prd.ToPyartRadar()
        self.assertIsInstance(radar, pyart.core.Radar)
        self.assertIn("reflectivity", radar.fields)
        self.assertEqual(radar.ngates, radar.range["data"].size)
        self.assertEqual(radar.nrays, 3264)

        legacy_radar = prd.to_pyart_radar(use_external=False, force_rebuild=True)
        self.assertIsInstance(legacy_radar, InternalRadar)

        native_radar = prd.to_pyart_radar(field_names=["dBZ"], use_external=False, force_rebuild=True)
        self.assertIsInstance(native_radar, InternalRadar)
        self.assertIn("reflectivity", native_radar.fields)
        self.assertGreater(native_radar.range["data"].size, radar.range["data"].size)
        self.assertEqual(native_radar.fields["reflectivity"]["data"].shape[0], radar.fields["reflectivity"]["data"].shape[0])

    def test_prd_xradar_exports_are_explicit(self):
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        sweeps = prd.to_xradar_sweeps(field_names=["dBZ"])
        self.assertIn("sweep_0", sweeps)
        self.assertIn("DBZ", sweeps["sweep_0"])
        self.assertEqual(sweeps["sweep_0"].attrs["range_mode"], "native")
        self.assertGreater(sweeps["sweep_0"].sizes["range"], prd.fields[0].sizes["range"])
        self.assertIn("sweep_mode", sweeps["sweep_0"])
        self.assertIn("follow_mode", sweeps["sweep_0"])
        self.assertIn("prt_mode", sweeps["sweep_0"])
        self.assertIn("sweep_fixed_angle", sweeps["sweep_0"])

        tree = prd.to_xradar(strict=True, field_names=["dBZ"])
        self.assertIn("sweep_0", tree.children)
        self.assertEqual(tree["sweep_0"].dataset.attrs["range_mode"], "native")
        self.assertEqual(tree.dataset.attrs["Conventions"], "Cf/Radial")
        self.assertIn("sweep_group_name", tree.dataset)
        self.assertIn("volume_number", tree.dataset)
        self.assertIn("instrument_type", tree.dataset)

    def test_real_pyart_interop(self):
        pyart = importlib.import_module("pyart")
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        radar = prd.to_pyart_radar(strict=True, force_rebuild=True)
        self.assertIsInstance(radar, pyart.core.Radar)
        self.assertIn("reflectivity", radar.fields)
        self.assertEqual(radar.nrays, 3264)
        self.assertEqual(radar.nsweeps, 9)
        self.assertEqual(radar.get_field(0, "reflectivity").shape[0], 361)
        self.assertGreater(radar.get_nyquist_vel(0), 0.0)

        display = pyart.graph.RadarDisplay(radar)
        self.assertEqual(type(display).__name__, "RadarDisplay")

        gatefilter = pyart.filters.GateFilter(radar)
        gatefilter.exclude_invalid("reflectivity")
        self.assertEqual(gatefilter.gate_excluded.shape, (radar.nrays, radar.ngates))

        grid = pyart.map.grid_from_radars(
            (radar,),
            grid_shape=(1, 21, 21),
            grid_limits=((1000.0, 1000.0), (-50_000.0, 50_000.0), (-50_000.0, 50_000.0)),
            fields=["reflectivity"],
        )
        self.assertIn("reflectivity", grid.fields)
        self.assertEqual(grid.fields["reflectivity"]["data"].shape, (1, 21, 21))

    def test_reader_level_pyart_export_delegates_to_prd(self):
        pyart = importlib.import_module("pyart")
        from pycwr.io import WSR98DFile

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        reader = WSR98DFile.WSR98D2NRadar(WSR98DFile.WSR98DBaseData(str(sample)))
        radar = reader.to_pyart_radar(strict=True)
        self.assertIsInstance(radar, pyart.core.Radar)
        self.assertIn("reflectivity", radar.fields)
        legacy_radar = reader.to_pyart_radar(use_external=False, force_rebuild=True)
        self.assertEqual(type(legacy_radar).__name__, "Radar")

    def test_ordered_az_cache_and_sorting(self):
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        ordered = prd.ordered_az()
        self.assertIs(ordered, prd.ordered_az())
        self.assertGreater(len(ordered.fields), 0)
        for sweep in ordered.fields:
            azimuth = sweep.azimuth.values
            self.assertTrue(np.all(azimuth[:-1] <= azimuth[1:]))


if __name__ == "__main__":
    unittest.main()
