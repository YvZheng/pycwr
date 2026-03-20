"""Public-API sample tests that double as copyable usage examples."""

import tempfile
import unittest
import os
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")


class PublicApiSampleTests(unittest.TestCase):
    @staticmethod
    def _repo_root():
        return Path(__file__).resolve().parents[2]

    @classmethod
    def _sample_z9046(cls):
        return (
            cls._repo_root()
            / "Z9046"
            / "Z_RADR_I_Z9046_20260317065928_O_DOR_SAD_CAP_FMT.bin.bz2"
        )

    @classmethod
    def _sample_z9047(cls):
        return (
            cls._repo_root()
            / "Z9047"
            / "Z_RADR_I_Z9047_20260317065837_O_DOR_SAD_CAP_FMT.bin.bz2"
        )

    def test_prd_summary_and_available_fields(self):
        from pycwr.io import read_auto

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        summary = prd.summary()

        self.assertEqual(summary["sitename"], prd.sitename)
        self.assertEqual(summary["scan_type"], "ppi")
        self.assertEqual(summary["nsweeps"], int(prd.nsweeps))
        self.assertIn("dBZ", summary["fields"])
        self.assertIn("dBZ", prd.available_fields(sweep=0))
        self.assertEqual(len(summary["sweeps"]), int(prd.nsweeps))
        self.assertGreater(summary["sweeps"][0]["aligned_max_range_m"], 0.0)

    def test_format_specific_reader_example_matches_read_auto(self):
        from pycwr.io import read_WSR98D, read_auto

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        auto_prd = read_auto(str(sample))
        direct_prd = read_WSR98D(str(sample))

        self.assertEqual(int(auto_prd.nsweeps), int(direct_prd.nsweeps))
        self.assertEqual(int(auto_prd.nrays), int(direct_prd.nrays))
        np.testing.assert_allclose(
            np.asarray(auto_prd.fields[0]["dBZ"].values, dtype=np.float64),
            np.asarray(direct_prd.fields[0]["dBZ"].values, dtype=np.float64),
            rtol=0.0,
            atol=1e-6,
            equal_nan=True,
        )

    def test_native_reflectivity_example(self):
        from pycwr.io import read_auto

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        default_field = prd.get_sweep_field(0, "dBZ")
        aligned = prd.get_sweep_field(0, "dBZ", range_mode="aligned")
        native = prd.get_sweep_field(0, "dBZ", range_mode="native")

        self.assertTrue(prd.has_extended_field(0, "dBZ"))
        self.assertEqual(default_field.shape, native.shape)
        self.assertEqual(aligned.shape[0], native.shape[0])
        self.assertGreater(native.shape[1], aligned.shape[1])
        self.assertGreater(float(native.range.values[-1]), float(aligned.range.values[-1]))

    def test_product_generation_example(self):
        from pycwr.io import read_auto

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        x = np.arange(-20_000.0, 20_001.0, 10_000.0)
        y = np.arange(-20_000.0, 20_001.0, 10_000.0)
        product_levels = np.array([1000.0, 2000.0, 3000.0], dtype=np.float64)

        prd.add_product_CR_xy(x, y)
        prd.add_product_CAPPI_xy(x, y, 3000.0)
        prd.add_product_CAPPI_3d_xy(x, y, np.array([1000.0, 3000.0], dtype=np.float64))
        prd.add_product_VIL_xy(x, y, product_levels)
        prd.add_product_ET_xy(x, y, product_levels)

        self.assertIn("CR_native", prd.product)
        self.assertIn("CAPPI_3000_native", prd.product)
        self.assertIn("CAPPI_3D_native", prd.product)
        self.assertIn("VIL_native", prd.product)
        self.assertIn("ET_native", prd.product)
        self.assertIn("ET_TOPPED_native", prd.product)

    def test_export_examples_without_external_pyart_dependency(self):
        from pycwr.io import read_auto

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        legacy_radar = prd.to_pyart_radar(use_external=False, field_names=["dBZ"], force_rebuild=True)
        sweeps = prd.to_xradar_sweeps(field_names=["dBZ"], force_rebuild=True)

        self.assertIn("reflectivity", legacy_radar.fields)
        self.assertGreater(
            legacy_radar.range["data"].size,
            prd.fields[0].sizes["range"],
        )
        self.assertIn("sweep_0", sweeps)
        self.assertIn("DBZ", sweeps["sweep_0"].data_vars)
        self.assertEqual(sweeps["sweep_0"].attrs["range_mode"], "native")

    def test_sorted_sweep_field_cache_on_real_sample(self):
        from pycwr.io import read_auto

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        first = prd.get_sweep_field(0, "dBZ", sort_by_azimuth=True)
        second = prd.get_sweep_field(0, "dBZ", sort_by_azimuth=True)

        self.assertIs(first, second)
        self.assertEqual(first.dims[0], "azimuth")
        self.assertTrue((first.azimuth.values[:-1] <= first.azimuth.values[1:]).all())

    def test_real_sample_regression_fingerprint(self):
        from pycwr.io import read_auto

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        field = np.asarray(prd.fields[0]["dBZ"].values, dtype=np.float64)

        self.assertEqual(int(prd.nsweeps), 9)
        self.assertEqual(int(prd.nrays), 3264)
        self.assertEqual(field.shape, (361, 920))
        self.assertEqual(int(np.isfinite(field).sum()), 191144)
        self.assertAlmostEqual(float(np.nanmean(field)), 15.491404386221905, places=6)
        self.assertAlmostEqual(float(np.nanstd(field)), 7.316992112108141, places=6)
        np.testing.assert_allclose(
            np.asarray(prd.scan_info["fixed_angle"].values, dtype=np.float64),
            np.array([0.5, 1.5, 2.4000001, 3.4000001, 4.30000019, 6.0, 9.89999962, 14.60000038, 19.5]),
            rtol=0.0,
            atol=1e-6,
        )

    def test_extract_rhi_matches_sorted_field_interpolation(self):
        from pycwr.io import read_auto

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        target_azimuth = 0.0
        section = prd.extract_rhi(azimuth=target_azimuth, field_name="dBZ")

        expected_rows = []
        expected_ranges = []
        expected_elevations = []
        for sweep in np.asarray(prd._fixed_angle_sort_index, dtype=np.int32):
            sweep_field = prd.get_sweep_field(sweep, "dBZ", sort_by_azimuth=True)
            sweep_ranges = np.asarray(sweep_field.range.values, dtype=np.float64)
            azimuth_row = np.full(sweep_ranges.shape, target_azimuth, dtype=np.float64)
            expected_rows.append(prd._interpolate_ppi_strip(sweep_field, azimuth_row, sweep_ranges))
            expected_ranges.append(sweep_ranges)
            expected_elevations.append(
                np.full(sweep_ranges.shape, float(prd.scan_info.fixed_angle.values[sweep]), dtype=np.float64)
            )

        np.testing.assert_allclose(
            np.asarray(section["dBZ"].values, dtype=np.float64),
            prd._pad_section_rows(expected_rows),
            rtol=0.0,
            atol=1e-6,
            equal_nan=True,
        )
        np.testing.assert_allclose(
            np.asarray(section["range"].values, dtype=np.float64),
            prd._pad_section_rows(expected_ranges),
            rtol=0.0,
            atol=1e-6,
            equal_nan=True,
        )
        np.testing.assert_allclose(
            np.asarray(section["elevation"].values, dtype=np.float64),
            prd._pad_section_rows(expected_elevations),
            rtol=0.0,
            atol=1e-6,
            equal_nan=True,
        )

    def test_easy_plot_quickstart_smoke(self):
        from pycwr.draw import plot, plot_ppi, plot_rhi, plot_section
        from pycwr.io import read_auto

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        with tempfile.TemporaryDirectory() as tmpdir:
            ppi_path = Path(tmpdir) / "ppi.png"
            section_path = Path(tmpdir) / "section.png"
            rhi_path = Path(tmpdir) / "rhi.png"

            ppi = plot_ppi(prd, field="dBZ", sweep=0, save=ppi_path)
            section = plot_section(prd, start=(-30, 0), end=(30, 0), field="dBZ", save=section_path)
            rhi = plot_rhi(prd, azimuth=0.0, field="dBZ", save=rhi_path)
            quicklook = plot(prd, kind="ppi", field="dBZ")

            self.assertIsNotNone(ppi.artist)
            self.assertIsNotNone(section.artist)
            self.assertIsNotNone(rhi.artist)
            self.assertIsNotNone(quicklook.artist)
            self.assertTrue(ppi_path.exists())
            self.assertTrue(section_path.exists())
            self.assertTrue(rhi_path.exists())

    def test_single_radar_wind_retrieval_example(self):
        from pycwr.io import read_auto
        from pycwr.retrieve import retrieve_vad, retrieve_vvp

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        vad = retrieve_vad(prd, sweeps=0, max_range_km=40.0, gate_step=4)
        vvp = retrieve_vvp(prd, sweep=0, max_range_km=20.0, az_num=91, bin_num=5, azimuth_step=12, range_step=4)
        profile = prd.retrieve_vwp(sweeps=[0, 1, 2], max_range_km=40.0, gate_step=4, height_step=500.0)
        stored = prd.add_product_VWP(sweeps=[0, 1, 2], max_range_km=40.0, gate_step=4, height_step=500.0)

        self.assertEqual(vad.attrs["method"], "VAD")
        self.assertEqual(vvp.attrs["method"], "VVP")
        self.assertEqual(profile.attrs["method"], "VWP")
        self.assertEqual(stored.attrs["method"], "VWP")
        self.assertGreater(int(np.isfinite(vad["u"].values).sum()), 0)
        self.assertGreater(int(np.isfinite(vvp["u"].values).sum()), 0)
        self.assertGreater(int(np.isfinite(profile["u"].values).sum()), 0)
        self.assertIn("VWP_u", prd.product)
        self.assertIn("z_vwp", prd.product.coords)

    def test_web_viewer_metadata_on_real_sample(self):
        try:
            from pycwr.GraphicalInterface import create_app
        except ImportError as exc:  # pragma: no cover
            self.skipTest(f"Flask viewer dependencies are unavailable: {exc}")

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        app = create_app(allowed_roots=[sample.parent], auth_token="sample-token")
        client = app.test_client()
        response = client.get(
            "/api/metadata",
            query_string={"path": str(sample), "token": "sample-token"},
        )

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertTrue(payload["ok"])
        self.assertEqual(payload["name"], sample.name)
        self.assertIn("dBZ", payload["fields"])
        self.assertEqual(payload["scan_type"], "ppi")

    def test_real_network_workflow_small_grid(self):
        from pycwr.interp import run_radar_network_3d

        if os.environ.get("PYCWR_RUN_SLOW_SAMPLE_TESTS") != "1":
            self.skipTest("set PYCWR_RUN_SLOW_SAMPLE_TESTS=1 to run slow real-sample network tests")

        sample_a = self._sample_z9046()
        sample_b = self._sample_z9047()
        if not sample_a.exists() or not sample_b.exists():
            self.skipTest("network sample radar files are not available in this workspace")

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "network.nc"
            dataset = run_radar_network_3d(
                target_time="2026-03-17T07:00:00",
                radar_dirs=[str(sample_a.parent), str(sample_b.parent)],
                lon_min=118.5,
                lon_max=118.5,
                lat_min=32.0,
                lat_max=32.0,
                lon_res_deg=0.05,
                lat_res_deg=0.05,
                level_heights=[1000.0],
                field_names=["dBZ"],
                output_path=str(output_path),
                parallel=False,
                use_qc=False,
            )

            self.assertIn("dBZ", dataset.data_vars)
            self.assertEqual(dataset["dBZ"].dims, ("z", "lat", "lon"))
            self.assertEqual(dataset.attrs["radar_count"], 2)
            self.assertTrue(output_path.exists())


if __name__ == "__main__":
    unittest.main()
