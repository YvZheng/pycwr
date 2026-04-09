"""Public-API sample tests that double as copyable usage examples."""

import tempfile
import unittest
import os
import importlib
import subprocess
import sys
import warnings
from pathlib import Path
from unittest import mock

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

    @classmethod
    def _sample_z9280(cls):
        return cls._repo_root() / "Z_RADR_I_Z9280_20250402034600_O_DOR-CUT_SAD_CAP_5_1_FMT.bin.bz2"

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

    def test_lonlat_product_generation_example(self):
        from pycwr.io import read_auto

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        lon0 = float(prd.scan_info["longitude"].values)
        lat0 = float(prd.scan_info["latitude"].values)
        lon = np.array([lon0 - 0.05, lon0, lon0 + 0.05], dtype=np.float64)
        lat = np.array([lat0 - 0.05, lat0, lat0 + 0.05], dtype=np.float64)
        product_levels = np.array([1000.0, 2000.0, 3000.0], dtype=np.float64)

        prd.add_product_CR_lonlat(lon, lat)
        prd.add_product_CAPPI_lonlat(lon, lat, 3000.0)
        prd.add_product_VIL_lonlat(lon, lat, product_levels)
        prd.add_product_ET_lonlat(lon, lat, product_levels)

        self.assertIn("CR_geo_native", prd.product)
        self.assertIn("CAPPI_geo_3000_native", prd.product)
        self.assertIn("VIL_geo_native", prd.product)
        self.assertIn("ET_geo_native", prd.product)
        self.assertIn("ET_TOPPED_geo_native", prd.product)
        self.assertEqual(prd.product["CR_geo_native"].dims, ("lon_cr", "lat_cr"))
        self.assertEqual(prd.product["CAPPI_geo_3000_native"].dims, ("lon_cappi_3000", "lat_cappi_3000"))
        self.assertEqual(prd.product["VIL_geo_native"].dims, ("lon_vil", "lat_vil"))
        self.assertEqual(prd.product["ET_geo_native"].dims, ("lon_et", "lat_et"))
        self.assertEqual(prd.product["ET_TOPPED_geo_native"].dims, ("lon_et", "lat_et"))

    def test_export_examples_without_external_pyart_dependency(self):
        from pycwr.io import read_auto

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        legacy_radar = prd.to_pyart_radar(use_external=False, field_names=["dBZ"], force_rebuild=True)

        self.assertIn("reflectivity", legacy_radar.fields)
        self.assertGreater(
            legacy_radar.range["data"].size,
            prd.fields[0].sizes["range"],
        )
        try:
            importlib.import_module("xradar")
        except ImportError:
            return

        sweeps = prd.to_xradar_sweeps(field_names=["dBZ"], force_rebuild=True)
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

    def test_z9280_cut_sample_is_readable(self):
        from pycwr.io import read_auto

        sample = self._sample_z9280()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))

        self.assertEqual(prd.scan_info.scan_type.item(), "ppi")
        self.assertEqual(int(prd.nsweeps), 1)
        self.assertEqual(int(prd.nrays), 366)
        self.assertEqual(
            sorted(prd.available_fields()),
            ["CC", "KDP", "PhiDP", "SNRH", "ZDR", "dBT", "dBZ"],
        )
        summary = prd.sweep_summary()
        self.assertEqual(len(summary), 1)
        self.assertAlmostEqual(float(summary[0]["fixed_angle"]), 0.5, places=6)

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

    def test_direct_graphmap_import_auto_loads_pycwr_colormap(self):
        try:
            import cartopy  # noqa: F401
        except ImportError as exc:  # pragma: no cover
            self.skipTest(f"map plotting dependencies are unavailable: {exc}")

        sample = self._sample_z9280()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        script = """
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cartopy import crs as ccrs
from pycwr.io import read_auto
from pycwr.draw.RadarPlot import GraphMap

sample = r\"\"\"%s\"\"\"
prd = read_auto(sample)
fig = plt.figure()
try:
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    GraphMap(prd, ccrs.PlateCarree()).plot_ppi_map(ax, 0, "dBZ")
finally:
    plt.close(fig)
""" % sample
        result = subprocess.run(
            [sys.executable, "-c", script],
            cwd=str(Path(__file__).resolve().parents[1]),
            env=os.environ.copy(),
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr or result.stdout)

    def test_missing_velocity_field_reports_field_error_on_z9280(self):
        try:
            import matplotlib.pyplot as plt
            from cartopy import crs as ccrs
        except ImportError as exc:  # pragma: no cover
            self.skipTest(f"map plotting dependencies are unavailable: {exc}")

        from pycwr.draw.RadarPlot import GraphMap
        from pycwr.io import read_auto

        sample = self._sample_z9280()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        fig = plt.figure()
        try:
            ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
            with self.assertRaisesRegex(KeyError, "field V not found"):
                GraphMap(prd, ccrs.PlateCarree()).plot_ppi_map(ax, 0, "V")
        finally:
            plt.close(fig)

    def test_plot_ppi_default_cmap_strategy_smoke(self):
        from pycwr.draw import plot_ppi
        from pycwr.io import read_auto

        sample = self._sample_z9280()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "z9280_ppi.png"
            result = plot_ppi(prd, field="dBZ", sweep=0, save=output)
            self.assertTrue(output.exists())

        self.assertIsNotNone(result.artist)

    def test_plot_ppi_map_default_cmap_strategy_smoke(self):
        try:
            import cartopy  # noqa: F401
        except ImportError as exc:  # pragma: no cover
            self.skipTest(f"map plotting dependencies are unavailable: {exc}")

        from pycwr.draw import plot_ppi_map
        from pycwr.io import read_auto

        sample = self._sample_z9280()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "z9280_ppi_map.png"
            result = plot_ppi_map(prd, field="dBZ", sweep=0, save=output)
            self.assertTrue(output.exists())

        self.assertIsNotNone(result.artist)

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

    def test_apply_map_axes_expands_degenerate_extent_without_cartopy_warning(self):
        try:
            import matplotlib.pyplot as plt
            from cartopy import crs as ccrs
        except ImportError as exc:  # pragma: no cover
            self.skipTest(f"map plotting dependencies are unavailable: {exc}")

        from pycwr.draw._plot_core import MapOptions, apply_map_axes

        fig = plt.figure()
        try:
            ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                apply_map_axes(
                    ax,
                    (120.0, 120.0, 30.0, 30.0),
                    MapOptions(data_crs=ccrs.PlateCarree()),
                )
            messages = "\n".join(str(item.message) for item in caught)
            self.assertNotIn("Attempting to set identical low and high xlims", messages)
            self.assertNotIn("Attempting to set identical low and high ylims", messages)
        finally:
            plt.close(fig)

    def test_single_radar_wind_retrieval_example(self):
        from pycwr.io import read_auto
        from pycwr.retrieve import retrieve_vad, retrieve_vvp, retrieve_wind_volume_xy

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        vad = retrieve_vad(prd, sweeps=0, max_range_km=40.0, gate_step=4)
        vvp = retrieve_vvp(prd, sweep=0, max_range_km=20.0, az_num=91, bin_num=5, azimuth_step=12, range_step=4)
        profile = prd.retrieve_vwp(sweeps=[0, 1, 2], max_range_km=40.0, gate_step=4, height_step=500.0)
        stored = prd.add_product_VWP(sweeps=[0, 1, 2], max_range_km=40.0, gate_step=4, height_step=500.0)
        wind = retrieve_wind_volume_xy(
            prd,
            XRange=np.array([0.0], dtype=np.float64),
            YRange=np.array([0.0], dtype=np.float64),
            level_heights=np.array([500.0, 1000.0], dtype=np.float64),
            max_range_km=20.0,
            az_num=31,
            bin_num=5,
            azimuth_step=24,
            range_step=8,
            horizontal_radius_m=8_000.0,
            horizontal_min_neighbors=3,
            vertical_tolerance_m=400.0,
        )
        stored_wind = prd.add_product_WIND_VOLUME_xy(
            XRange=np.array([0.0], dtype=np.float64),
            YRange=np.array([0.0], dtype=np.float64),
            level_heights=np.array([500.0, 1000.0], dtype=np.float64),
            max_range_km=20.0,
            az_num=31,
            bin_num=5,
            azimuth_step=24,
            range_step=8,
            horizontal_radius_m=8_000.0,
            horizontal_min_neighbors=3,
            vertical_tolerance_m=400.0,
        )

        self.assertEqual(vad.attrs["method"], "VAD")
        self.assertEqual(vvp.attrs["method"], "VVP")
        self.assertEqual(profile.attrs["method"], "VWP")
        self.assertEqual(stored.attrs["method"], "VWP")
        self.assertEqual(wind.attrs["method"], "single_radar_horizontal_wind_volume")
        self.assertEqual(stored_wind.attrs["method"], "single_radar_horizontal_wind_volume")
        self.assertEqual(wind.attrs["selection_mode"], "auto")
        self.assertGreater(int(np.isfinite(vad["u"].values).sum()), 0)
        self.assertGreater(int(np.isfinite(vvp["u"].values).sum()), 0)
        self.assertGreater(int(np.isfinite(profile["u"].values).sum()), 0)
        self.assertGreater(int(np.isfinite(wind["u"].values).sum()), 0)
        self.assertIn("VWP_u", prd.product)
        self.assertIn("z_vwp", prd.product.coords)
        self.assertIn("WIND_VOLUME_u", prd.product)
        self.assertIn("WIND_VOLUME_quality_flag", prd.product)
        self.assertIn("z_wind", prd.product.coords)

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

    def test_web_viewer_launch_smoke(self):
        try:
            from pycwr.GraphicalInterface import launch
        except ImportError as exc:  # pragma: no cover
            self.skipTest(f"Flask viewer dependencies are unavailable: {exc}")

        fake_app = mock.Mock()
        with mock.patch("pycwr.GraphicalInterface.web_app.create_app", return_value=fake_app):
            launch(host="127.0.0.1", port=8765, open_browser=False)

        fake_app.run.assert_called_once_with(
            host="127.0.0.1",
            port=8765,
            debug=False,
            use_reloader=False,
        )

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
