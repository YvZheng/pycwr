"""Network interpolation examples and regression tests."""

import tempfile
import unittest
from pathlib import Path
from unittest import mock
import json

import numpy as np
import xarray as xr


class RadarInterpWorkflowTests(unittest.TestCase):
    def test_parse_radar_time_from_filename(self):
        from pycwr.interp.RadarInterp import parse_radar_time_from_filename

        parsed = parse_radar_time_from_filename(
            "Z_RADR_I_Z9046_20260317065928_O_DOR_SAD_CAP_FMT.bin.bz2"
        )
        self.assertEqual(parsed.isoformat(), "2026-03-17T06:59:28")

    def test_select_radar_files_uses_nearest_time_within_tolerance(self):
        from pycwr.interp.RadarInterp import select_radar_files

        with tempfile.TemporaryDirectory() as tmpdir:
            radar_a = Path(tmpdir) / "Z9001"
            radar_b = Path(tmpdir) / "Z9002"
            radar_a.mkdir()
            radar_b.mkdir()
            (radar_a / "A_20260317065800_test.bin").write_text("", encoding="ascii")
            (radar_a / "A_20260317070100_test.bin").write_text("", encoding="ascii")
            (radar_b / "B_20260317070030_test.bin").write_text("", encoding="ascii")
            (radar_b / "B_20260317073000_test.bin").write_text("", encoding="ascii")

            selected = select_radar_files(
                [str(radar_a), str(radar_b)],
                "2026-03-17T07:00:00",
                tolerance_minutes=5,
            )

        self.assertEqual(len(selected), 2)
        selected_paths = {Path(item["path"]).name for item in selected}
        self.assertEqual(
            selected_paths,
            {"A_20260317070100_test.bin", "B_20260317070030_test.bin"},
        )

    def test_build_latlon_grid_returns_lat_lon_mesh(self):
        from pycwr.interp.RadarInterp import build_latlon_grid

        lon, lat, grid_lon, grid_lat = build_latlon_grid(120.0, 120.1, 30.0, 30.1, 0.1, 0.1)

        self.assertTrue(np.allclose(lon, np.array([120.0, 120.1])))
        self.assertTrue(np.allclose(lat, np.array([30.0, 30.1])))
        self.assertEqual(grid_lon.shape, (2, 2))
        self.assertEqual(grid_lat.shape, (2, 2))
        self.assertTrue(np.allclose(grid_lon[0], lon))
        self.assertTrue(np.allclose(grid_lat[:, 0], lat))

    def test_compose_network_volume_uses_exponential_weighting(self):
        from pycwr.interp.RadarInterp import compose_network_volume

        grid_x = np.array([[0.0]], dtype=np.float64)
        grid_y = np.array([[0.0]], dtype=np.float64)
        radar_sites = np.array([[-1000.0, 0.0], [1000.0, 0.0]], dtype=np.float64)
        radar_volumes = [
            np.array([[[10.0]]], dtype=np.float64),
            np.array([[[30.0]]], dtype=np.float64),
        ]

        composite, valid_count = compose_network_volume(
            radar_volumes,
            radar_sites,
            grid_x,
            grid_y,
            fillvalue=-999.0,
            influence_radius_m=300000.0,
        )

        self.assertTrue(np.isclose(composite[0, 0, 0], 20.0, rtol=0.0, atol=1e-6))
        self.assertEqual(int(valid_count[0, 0, 0]), 2)

    def test_grid_single_radar_masks_values_outside_source_range(self):
        from pycwr.interp.RadarInterp import grid_single_radar_to_latlon_3d

        class FakePRD(object):
            def get_vol_data(self, field_name="dBZ", fillvalue=-999.0, range_mode="aligned", max_range_km=None):
                self.vol = (
                    [np.array([0.0], dtype=np.float64)],
                    [np.array([1000.0], dtype=np.float64)],
                    np.array([1.0], dtype=np.float64),
                    [np.array([[10.0]], dtype=np.float64)],
                    100.0,
                    120.0,
                    30.0,
                )

        with mock.patch(
            "pycwr.interp.RadarInterp.get_CAPPI_3d",
            return_value=np.array([[[9999.0]]], dtype=np.float64),
        ):
            grid_value = grid_single_radar_to_latlon_3d(
                FakePRD(),
                "dBZ",
                np.array([[0.0]], dtype=np.float64),
                np.array([[0.0]], dtype=np.float64),
                np.array([1000.0], dtype=np.float64),
                0.0,
                0.0,
                fillvalue=-999.0,
                range_mode="native",
                max_range_km=460.0,
            )

        self.assertEqual(float(grid_value[0, 0, 0]), -999.0)

    def test_derive_et_reports_topped_flag(self):
        from pycwr.interp.RadarInterp import _derive_et

        volume = np.array(
            [
                [[10.0, 20.0]],
                [[30.0, 25.0]],
                [[40.0, 19.0]],
            ],
            dtype=np.float64,
        )
        heights = np.array([1000.0, 2000.0, 3000.0], dtype=np.float64)

        et, topped = _derive_et(
            volume,
            heights,
            fillvalue=-999.0,
            threshold_dbz=18.0,
            return_topped=True,
        )

        self.assertTrue(np.isclose(et[0, 0], 3000.0))
        self.assertEqual(int(topped[0, 0]), 1)
        self.assertTrue(np.isclose(et[0, 1], 3000.0))
        self.assertEqual(int(topped[0, 1]), 1)

    def test_run_radar_network_3d_supports_multiple_fields_and_netcdf(self):
        from pycwr.interp.RadarInterp import radar_network_3d_to_netcdf, run_radar_network_3d

        class FakePRD(object):
            def __init__(self, radar_id):
                self.radar_id = radar_id
                self.effective_earth_radius = None
                self.scan_info = xr.Dataset(
                    data_vars={
                        "longitude": xr.DataArray(120.0),
                        "latitude": xr.DataArray(30.0),
                    }
                )

            def get_vol_data(self, field_name="dBZ", fillvalue=-999.0, range_mode="aligned", max_range_km=None):
                self.vol = (
                    [
                        np.array([0.0, 90.0], dtype=np.float64),
                        np.array([0.0, 90.0], dtype=np.float64),
                    ],
                    [
                        np.array([1000.0, 2000.0], dtype=np.float64),
                        np.array([1000.0, 2000.0], dtype=np.float64),
                    ],
                    np.array([1.0, 2.0], dtype=np.float64),
                    [
                        np.full((2, 2), 10.0, dtype=np.float64),
                        np.full((2, 2), 12.0, dtype=np.float64),
                    ],
                    100.0,
                    120.0,
                    30.0,
                )

        field_lookup = {
            ("Z9001", "dBZ"): 10.0,
            ("Z9002", "dBZ"): 30.0,
            ("Z9001", "ZDR"): 1.0,
            ("Z9002", "ZDR"): 3.0,
        }

        def fake_read_auto(path, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
            return FakePRD(Path(path).parent.name)

        def fake_get_radar_info(path):
            return 30.0, 120.0, 100.0, 2.8

        def fake_grid_single(prd, field_name, grid_x, grid_y, level_heights, radar_x, radar_y, fillvalue=-999.0,
                             effective_earth_radius=None, range_mode="aligned", max_range_km=None,
                             blind_method="mask", return_metadata=False):
            value = field_lookup[(prd.radar_id, field_name)]
            grid = np.full((len(level_heights),) + grid_x.shape, value, dtype=np.float64)
            if return_metadata:
                return grid, {"actual_max_range_km": 460.0, "range_mode": range_mode}
            return grid

        with tempfile.TemporaryDirectory() as tmpdir:
            radar_a = Path(tmpdir) / "Z9001"
            radar_b = Path(tmpdir) / "Z9002"
            radar_a.mkdir()
            radar_b.mkdir()
            (radar_a / "A_20260317070000_test.bin").write_text("", encoding="ascii")
            (radar_b / "B_20260317070000_test.bin").write_text("", encoding="ascii")
            output_path = Path(tmpdir) / "network.nc"

            with mock.patch("pycwr.interp.RadarInterp.read_auto", side_effect=fake_read_auto):
                with mock.patch("pycwr.interp.RadarInterp.get_radar_info", side_effect=fake_get_radar_info):
                    with mock.patch(
                        "pycwr.interp.RadarInterp.grid_single_radar_to_latlon_3d",
                        side_effect=fake_grid_single,
                    ):
                        dataset = run_radar_network_3d(
                            target_time="2026-03-17T07:00:00",
                            radar_dirs=[str(radar_a), str(radar_b)],
                            lon_min=120.0,
                            lon_max=120.0,
                            lat_min=30.0,
                            lat_max=30.0,
                            lon_res_deg=0.01,
                            lat_res_deg=0.01,
                            level_heights=np.array([1000.0, 3000.0], dtype=np.float64),
                            field_names=["dBZ", "ZDR"],
                            parallel=False,
                            composite_method="exp_weighted",
                        )
                        written = radar_network_3d_to_netcdf(dataset, output_path)

            self.assertEqual(Path(written), output_path)
            self.assertIn("dBZ", dataset.data_vars)
            self.assertIn("ZDR", dataset.data_vars)
            self.assertEqual(dataset["dBZ"].dims, ("z", "lat", "lon"))
            self.assertTrue(np.allclose(dataset["dBZ"].values, 20.0))
            self.assertTrue(np.allclose(dataset["ZDR"].values, 2.0))
            self.assertEqual(dataset.attrs["composite_method"], "exp_weighted")
            self.assertEqual(json.loads(dataset.attrs["field_range_mode"])["dBZ"], "native")
            self.assertEqual(json.loads(dataset["dBZ"].attrs["actual_max_range_km_by_radar"])["Z9001"], 460.0)

            reopened = xr.load_dataset(output_path)
            try:
                self.assertEqual(tuple(reopened["dBZ"].dims), ("z", "lat", "lon"))
                self.assertTrue(np.allclose(reopened["ZDR"].values, 2.0))
            finally:
                reopened.close()

    def test_run_radar_network_3d_skips_missing_field_on_one_radar(self):
        from pycwr.interp.RadarInterp import run_radar_network_3d

        class FakePRD(object):
            def __init__(self, radar_id):
                self.radar_id = radar_id
                self.effective_earth_radius = None
                self.scan_info = xr.Dataset(
                    data_vars={
                        "longitude": xr.DataArray(120.0),
                        "latitude": xr.DataArray(30.0),
                    }
                )

            def get_vol_data(self, field_name="dBZ", fillvalue=-999.0, range_mode="aligned", max_range_km=None):
                self.vol = (
                    [
                        np.array([0.0, 90.0], dtype=np.float64),
                        np.array([0.0, 90.0], dtype=np.float64),
                    ],
                    [
                        np.array([1000.0, 2000.0], dtype=np.float64),
                        np.array([1000.0, 2000.0], dtype=np.float64),
                    ],
                    np.array([1.0, 2.0], dtype=np.float64),
                    [
                        np.full((2, 2), 10.0, dtype=np.float64),
                        np.full((2, 2), 12.0, dtype=np.float64),
                    ],
                    100.0,
                    120.0,
                    30.0,
                )

        def fake_read_auto(path, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
            return FakePRD(Path(path).parent.name)

        def fake_get_radar_info(path):
            return 30.0, 120.0, 100.0, 2.8

        def fake_grid_single(prd, field_name, grid_x, grid_y, level_heights, radar_x, radar_y, fillvalue=-999.0,
                             effective_earth_radius=None, range_mode="aligned", max_range_km=None,
                             blind_method="mask", return_metadata=False):
            if prd.radar_id == "Z9002" and field_name == "ZDR":
                raise KeyError(field_name)
            value = 20.0 if field_name == "dBZ" else 2.0
            grid = np.full((len(level_heights),) + grid_x.shape, value, dtype=np.float64)
            if return_metadata:
                return grid, {"actual_max_range_km": 460.0, "range_mode": range_mode}
            return grid

        with tempfile.TemporaryDirectory() as tmpdir:
            radar_a = Path(tmpdir) / "Z9001"
            radar_b = Path(tmpdir) / "Z9002"
            radar_a.mkdir()
            radar_b.mkdir()
            (radar_a / "A_20260317070000_test.bin").write_text("", encoding="ascii")
            (radar_b / "B_20260317070000_test.bin").write_text("", encoding="ascii")

            with mock.patch("pycwr.interp.RadarInterp.read_auto", side_effect=fake_read_auto):
                with mock.patch("pycwr.interp.RadarInterp.get_radar_info", side_effect=fake_get_radar_info):
                    with mock.patch(
                        "pycwr.interp.RadarInterp.grid_single_radar_to_latlon_3d",
                        side_effect=fake_grid_single,
                    ):
                        dataset = run_radar_network_3d(
                            target_time="2026-03-17T07:00:00",
                            radar_dirs=[str(radar_a), str(radar_b)],
                            lon_min=120.0,
                            lon_max=120.0,
                            lat_min=30.0,
                            lat_max=30.0,
                            lon_res_deg=0.01,
                            lat_res_deg=0.01,
                            level_heights=np.array([1000.0], dtype=np.float64),
                            field_names=["dBZ", "ZDR"],
                            parallel=False,
                            composite_method="exp_weighted",
                        )

        self.assertTrue(np.allclose(dataset["dBZ"].values, 20.0))
        self.assertTrue(np.allclose(dataset["ZDR"].values, 2.0))
        self.assertIsNone(json.loads(dataset["ZDR"].attrs["actual_max_range_km_by_radar"])["Z9002"])

    def test_run_radar_network_3d_uses_station_lookup_without_preread(self):
        from pycwr.interp.RadarInterp import run_radar_network_3d

        class FakePRD(object):
            def __init__(self, radar_id, lon, lat):
                self.radar_id = radar_id
                self.effective_earth_radius = None
                self.scan_info = xr.Dataset(
                    data_vars={
                        "longitude": xr.DataArray(lon),
                        "latitude": xr.DataArray(lat),
                        "altitude": xr.DataArray(100.0),
                    }
                )

            def get_vol_data(self, field_name="dBZ", fillvalue=-999.0, range_mode="aligned", max_range_km=None):
                self.vol = (
                    [np.array([0.0, 90.0], dtype=np.float64)],
                    [np.array([1000.0, 2000.0], dtype=np.float64)],
                    np.array([1.0], dtype=np.float64),
                    [np.full((2, 2), 10.0, dtype=np.float64)],
                    100.0,
                    float(self.scan_info["longitude"].values),
                    float(self.scan_info["latitude"].values),
                )

        read_calls = []

        def fake_read_auto(path, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
            radar_id = Path(path).parent.name
            read_calls.append((radar_id, station_lon, station_lat, station_alt))
            if radar_id == "Z9001":
                return FakePRD(radar_id, 120.0, 30.0)
            return FakePRD(radar_id, 121.0, 30.0)

        def fake_grid_single(prd, field_name, grid_x, grid_y, level_heights, radar_x, radar_y, fillvalue=-999.0,
                             effective_earth_radius=None, range_mode="aligned", max_range_km=None,
                             blind_method="mask", return_metadata=False):
            value = 10.0 if prd.radar_id == "Z9001" else 30.0
            grid = np.full((len(level_heights),) + grid_x.shape, value, dtype=np.float64)
            if return_metadata:
                return grid, {"actual_max_range_km": 230.0, "range_mode": range_mode}
            return grid

        station_lookup_calls = []

        def fake_get_radar_info(path):
            radar_id = Path(path).parent.name
            station_lookup_calls.append(radar_id)
            if radar_id == "Z9001":
                return 30.0, 120.0, 100.0, 2.8
            return 30.0, 121.0, 100.0, 2.8

        with tempfile.TemporaryDirectory() as tmpdir:
            radar_a = Path(tmpdir) / "Z9001"
            radar_b = Path(tmpdir) / "Z9002"
            radar_a.mkdir()
            radar_b.mkdir()
            (radar_a / "A_20260317070000_test.bin").write_text("", encoding="ascii")
            (radar_b / "B_20260317070000_test.bin").write_text("", encoding="ascii")

            with mock.patch("pycwr.interp.RadarInterp.read_auto", side_effect=fake_read_auto):
                with mock.patch("pycwr.interp.RadarInterp.get_radar_info", side_effect=fake_get_radar_info):
                    with mock.patch(
                        "pycwr.interp.RadarInterp.grid_single_radar_to_latlon_3d",
                        side_effect=fake_grid_single,
                    ):
                        with mock.patch(
                            "pycwr.interp.RadarInterp.grid_single_radar_cr_to_latlon",
                            return_value=np.zeros((1, 2), dtype=np.float64),
                        ):
                            with mock.patch(
                                "pycwr.interp.RadarInterp._compute_observability_mask",
                                return_value=np.ones((1, 1, 2), dtype=bool),
                            ):
                                dataset = run_radar_network_3d(
                                    target_time="2026-03-17T07:00:00",
                                    radar_dirs=[str(radar_a), str(radar_b)],
                                    lon_min=120.0,
                                    lon_max=121.0,
                                    lat_min=30.0,
                                    lat_max=30.0,
                                    lon_res_deg=1.0,
                                    lat_res_deg=1.0,
                                    level_heights=np.array([1000.0], dtype=np.float64),
                                    field_names=["dBZ"],
                                    output_products=[],
                                    plot_overview=False,
                                    plot_height_levels=[],
                                    parallel=False,
                                    composite_method="nearest",
                                )

        self.assertEqual(station_lookup_calls, ["Z9001", "Z9002"])
        self.assertEqual(read_calls, [("Z9001", 120.0, 30.0, 100.0), ("Z9002", 121.0, 30.0, 100.0)])
        self.assertTrue(np.allclose(dataset["dBZ"].values[0, 0], np.array([10.0, 30.0])))

    def test_run_radar_network_3d_can_apply_qc_before_mosaic(self):
        from pycwr.interp.RadarInterp import run_radar_network_3d

        class FakePRD(object):
            def __init__(self, radar_id):
                self.radar_id = radar_id
                self.nsweeps = 1
                self.effective_earth_radius = None
                self.scan_info = xr.Dataset(
                    data_vars={
                        "longitude": xr.DataArray(120.0),
                        "latitude": xr.DataArray(30.0),
                    }
                )
                self.fields = [
                    xr.Dataset(
                        data_vars={
                            "dBZ": xr.DataArray(np.ones((2, 2)), dims=("time", "range")),
                            "PhiDP": xr.DataArray(np.ones((2, 2)), dims=("time", "range")),
                        }
                    )
                ]

            def apply_dualpol_qc(
                self,
                sweeps=None,
                inplace=False,
                band="C",
                use_existing_kdp=True,
                clear_air_mode="label",
                clear_air_max_ref=15.0,
                clear_air_max_rhohv=0.97,
                clear_air_max_phidp_texture=10.0,
                clear_air_max_snr=20.0,
            ):
                self.fields[0]["Zc"] = xr.DataArray(np.ones((2, 2)), dims=("time", "range"))
                self.fields[0]["QC_MASK"] = xr.DataArray(np.ones((2, 2)), dims=("time", "range"))
                self.fields[0]["CLEAR_AIR_MASK"] = xr.DataArray(np.zeros((2, 2)), dims=("time", "range"))
                return self

            def get_vol_data(self, field_name="dBZ", fillvalue=-999.0, range_mode="aligned", max_range_km=None):
                self.vol = (
                    [
                        np.array([0.0, 90.0], dtype=np.float64),
                        np.array([0.0, 90.0], dtype=np.float64),
                    ],
                    [
                        np.array([1000.0, 2000.0], dtype=np.float64),
                        np.array([1000.0, 2000.0], dtype=np.float64),
                    ],
                    np.array([1.0, 2.0], dtype=np.float64),
                    [
                        np.full((2, 2), 5.0, dtype=np.float64),
                        np.full((2, 2), 6.0, dtype=np.float64),
                    ],
                    100.0,
                    120.0,
                    30.0,
                )

        def fake_read_auto(path, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
            return FakePRD(Path(path).parent.name)

        def fake_get_radar_info(path):
            return 30.0, 120.0, 100.0, 2.8

        observed_fields = []

        def fake_grid_single(prd, field_name, grid_x, grid_y, level_heights, radar_x, radar_y, fillvalue=-999.0,
                             effective_earth_radius=None, range_mode="aligned", max_range_km=None,
                             blind_method="mask", return_metadata=False):
            observed_fields.append(field_name)
            grid = np.full((len(level_heights),) + grid_x.shape, 5.0, dtype=np.float64)
            if return_metadata:
                return grid, {"actual_max_range_km": 460.0, "range_mode": range_mode}
            return grid

        with tempfile.TemporaryDirectory() as tmpdir:
            radar_a = Path(tmpdir) / "Z9001"
            radar_a.mkdir()
            (radar_a / "A_20260317070000_test.bin").write_text("", encoding="ascii")

            with mock.patch("pycwr.interp.RadarInterp.read_auto", side_effect=fake_read_auto):
                with mock.patch("pycwr.interp.RadarInterp.get_radar_info", side_effect=fake_get_radar_info):
                    with mock.patch(
                        "pycwr.interp.RadarInterp.grid_single_radar_to_latlon_3d",
                        side_effect=fake_grid_single,
                    ):
                        dataset = run_radar_network_3d(
                            target_time="2026-03-17T07:00:00",
                            radar_dirs=[str(radar_a)],
                            lon_min=120.0,
                            lon_max=120.0,
                            lat_min=30.0,
                            lat_max=30.0,
                            lon_res_deg=0.01,
                            lat_res_deg=0.01,
                            level_heights=np.array([1000.0], dtype=np.float64),
                            field_names=["dBZ"],
                            parallel=False,
                            use_qc=True,
                        )

        self.assertEqual(observed_fields, ["Zc", "QC_MASK"])
        self.assertTrue(dataset.attrs["use_qc"])
        self.assertTrue(json.loads(dataset.attrs["qc_applied_by_radar"])["Z9001"]["applied"])
        self.assertEqual(json.loads(dataset["dBZ"].attrs["source_field_name_by_radar"])["Z9001"], "Zc")
        self.assertEqual(dataset["dBZ"].attrs["qc_applied"], "true")

    def test_hybrid_blind_method_fills_near_radar_blind_zone(self):
        from pycwr.core.RadarGrid import get_CAPPI_3d

        vol_azimuth = [
            np.array([0.0, 90.0, 180.0, 270.0], dtype=np.float64),
            np.array([0.0, 90.0, 180.0, 270.0], dtype=np.float64),
        ]
        vol_range = [
            np.array([1000.0, 2000.0], dtype=np.float64),
            np.array([1000.0, 2000.0], dtype=np.float64),
        ]
        fix_elevation = np.array([1.0, 2.0], dtype=np.float64)
        vol_value = [
            np.full((4, 2), 10.0, dtype=np.float64),
            np.full((4, 2), 20.0, dtype=np.float64),
        ]
        grid_x = np.array([[500.0]], dtype=np.float64)
        grid_y = np.array([[0.0]], dtype=np.float64)
        level_heights = np.array([0.0], dtype=np.float64)

        masked = get_CAPPI_3d(
            vol_azimuth,
            vol_range,
            fix_elevation,
            vol_value,
            0.0,
            grid_x,
            grid_y,
            level_heights,
            fillvalue=-999.0,
            blind_method="mask",
        )
        hybrid = get_CAPPI_3d(
            vol_azimuth,
            vol_range,
            fix_elevation,
            vol_value,
            0.0,
            grid_x,
            grid_y,
            level_heights,
            fillvalue=-999.0,
            blind_method="hybrid",
        )

        self.assertEqual(float(masked[0, 0, 0]), -999.0)
        self.assertEqual(float(hybrid[0, 0, 0]), 10.0)

    def test_standard_product_derivations(self):
        from pycwr.interp.RadarInterp import _derive_cr, _derive_et, _derive_vil

        volume = np.array(
            [
                [[10.0]],
                [[30.0]],
                [[20.0]],
            ],
            dtype=np.float64,
        )
        heights = np.array([1000.0, 2000.0, 3000.0], dtype=np.float64)

        cr = _derive_cr(volume, fillvalue=-999.0)
        et = _derive_et(volume, heights, fillvalue=-999.0, threshold_dbz=18.0)
        vil = _derive_vil(volume, heights, fillvalue=-999.0, min_dbz=18.0, max_dbz_cap=56.0)

        self.assertEqual(float(cr[0, 0]), 30.0)
        self.assertTrue(2000.0 <= float(et[0, 0]) <= 3000.0)
        self.assertTrue(float(vil[0, 0]) > 0.0)

    def test_run_radar_network_3d_supports_external_config_file(self):
        from pycwr.interp.RadarInterp import run_radar_network_3d

        class FakePRD(object):
            def __init__(self, radar_id):
                self.radar_id = radar_id
                self.effective_earth_radius = None
                self.scan_info = xr.Dataset(
                    data_vars={
                        "longitude": xr.DataArray(120.0),
                        "latitude": xr.DataArray(30.0),
                    }
                )

            def get_vol_data(self, field_name="dBZ", fillvalue=-999.0, range_mode="aligned", max_range_km=None):
                self.vol = (
                    [
                        np.array([0.0, 90.0], dtype=np.float64),
                        np.array([0.0, 90.0], dtype=np.float64),
                    ],
                    [
                        np.array([1000.0, 2000.0], dtype=np.float64),
                        np.array([1000.0, 2000.0], dtype=np.float64),
                    ],
                    np.array([1.0, 2.0], dtype=np.float64),
                    [
                        np.full((2, 2), 7.0, dtype=np.float64),
                        np.full((2, 2), 8.0, dtype=np.float64),
                    ],
                    100.0,
                    120.0,
                    30.0,
                )

        observed_blind_methods = []

        def fake_read_auto(path, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
            return FakePRD(Path(path).parent.name)

        def fake_get_radar_info(path):
            return 30.0, 120.0, 100.0, 2.8

        def fake_grid_single(prd, field_name, grid_x, grid_y, level_heights, radar_x, radar_y, fillvalue=-999.0,
                             effective_earth_radius=None, range_mode="aligned", max_range_km=None,
                             blind_method="mask", return_metadata=False):
            observed_blind_methods.append(blind_method)
            grid = np.full((len(level_heights),) + grid_x.shape, 7.0, dtype=np.float64)
            if return_metadata:
                return grid, {"actual_max_range_km": 230.0, "range_mode": range_mode}
            return grid

        with tempfile.TemporaryDirectory() as tmpdir:
            radar_a = Path(tmpdir) / "Z9001"
            radar_a.mkdir()
            (radar_a / "A_20260317070000_test.bin").write_text("", encoding="ascii")
            config_path = Path(tmpdir) / "network.json"
            config_path.write_text(
                json.dumps(
                    {
                        "network": {
                            "radar_dirs": [str(radar_a)],
                            "field_names": ["dBZ"],
                            "lon_min": 120.0,
                            "lon_max": 120.0,
                            "lat_min": 30.0,
                            "lat_max": 30.0,
                            "lon_res_deg": 0.01,
                            "lat_res_deg": 0.01,
                            "level_heights": [1000.0],
                            "blind_method": "hybrid",
                            "composite_method": "exp_weighted",
                            "parallel": False,
                            "output_dir": tmpdir
                        }
                    }
                ),
                encoding="utf-8",
            )

            with mock.patch("pycwr.interp.RadarInterp.read_auto", side_effect=fake_read_auto):
                with mock.patch("pycwr.interp.RadarInterp.get_radar_info", side_effect=fake_get_radar_info):
                    with mock.patch(
                        "pycwr.interp.RadarInterp.grid_single_radar_to_latlon_3d",
                        side_effect=fake_grid_single,
                    ):
                        dataset = run_radar_network_3d(
                            target_time="2026-03-17T07:00:00",
                            config_path=str(config_path),
                        )

        self.assertEqual(observed_blind_methods, ["hybrid"])
        self.assertEqual(dataset.attrs["blind_method"], "hybrid")
        self.assertEqual(dataset.attrs["network_config_path"], str(config_path))
        self.assertTrue(np.allclose(dataset["dBZ"].values, 7.0))

    def test_reference_plot_resolution_guard_rejects_coarse_grid(self):
        from pycwr.interp.RadarInterp import _validate_reference_plot_resolution

        with self.assertRaises(ValueError):
            _validate_reference_plot_resolution(0.05, 0.01)

    def test_build_reference_cr_colormap_matches_reference_bins(self):
        from pycwr.interp.RadarInterp import _REFERENCE_CR_BOUNDS, _REFERENCE_CR_COLORS, _build_reference_cr_colormap

        cmap, norm = _build_reference_cr_colormap()

        self.assertEqual(len(_REFERENCE_CR_COLORS), 14)
        self.assertTrue(np.allclose(_REFERENCE_CR_BOUNDS, np.arange(-5.0, 70.0, 5.0)))
        self.assertEqual(cmap.N, 14)
        self.assertEqual(int(norm.boundaries[0]), -5)
        self.assertEqual(int(norm.boundaries[-1]), 65)

    def test_plot_reference_cr_field_outputs_reference_canvas(self):
        from datetime import datetime
        from PIL import Image
        from pycwr.interp.RadarInterp import _plot_reference_cr_field
        from pycwr.draw._plot_core import ccrs

        if ccrs is None:
            self.skipTest("cartopy is not available")

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "reference.png"
            written = _plot_reference_cr_field(
                np.array([120.0, 121.0], dtype=np.float64),
                np.array([30.0, 31.0], dtype=np.float64),
                np.array([[10.0, 20.0], [30.0, np.nan]], dtype=np.float64),
                datetime.fromisoformat("2026-03-17T07:00:00"),
                output_path,
            )

            with Image.open(written) as img:
                self.assertEqual(img.size, (920, 790))

    def test_compose_network_volume_quality_weighted_prefers_higher_quality(self):
        from pycwr.interp.RadarInterp import compose_network_volume

        grid_x = np.array([[0.0]], dtype=np.float64)
        grid_y = np.array([[0.0]], dtype=np.float64)
        radar_sites = np.array([[-1000.0, 0.0], [1000.0, 0.0]], dtype=np.float64)
        radar_volumes = [
            np.array([[[10.0]]], dtype=np.float64),
            np.array([[[30.0]]], dtype=np.float64),
        ]
        quality_weights = [
            np.array([[[1.0]]], dtype=np.float64),
            np.array([[[0.25]]], dtype=np.float64),
        ]

        composite, valid_count = compose_network_volume(
            radar_volumes,
            radar_sites,
            grid_x,
            grid_y,
            fillvalue=-999.0,
            method="quality_weighted",
            influence_radius_m=300000.0,
            quality_weights=quality_weights,
        )

        self.assertTrue(np.isclose(composite[0, 0, 0], 14.0, rtol=0.0, atol=1e-6))
        self.assertEqual(int(valid_count[0, 0, 0]), 2)

    def test_compute_observability_mask_supports_single_sweep(self):
        from pycwr.interp.RadarInterp import _compute_observability_mask

        mask = _compute_observability_mask(
            [np.array([1000.0, 2000.0], dtype=np.float64)],
            np.array([1.0], dtype=np.float64),
            0.0,
            np.array([[1500.0]], dtype=np.float64),
            np.array([[0.0]], dtype=np.float64),
            np.array([26.0], dtype=np.float64),
            0.0,
            0.0,
            beam_widths=np.array([2.0], dtype=np.float64),
        )

        self.assertTrue(bool(mask[0, 0, 0]))

    def test_run_radar_network_3d_rejects_velocity_field(self):
        from pycwr.interp.RadarInterp import run_radar_network_3d

        with self.assertRaises(ValueError):
            run_radar_network_3d(
                target_time="2026-03-17T07:00:00",
                radar_dirs=["/tmp/does-not-matter"],
                lon_min=120.0,
                lon_max=120.0,
                lat_min=30.0,
                lat_max=30.0,
                lon_res_deg=0.01,
                lat_res_deg=0.01,
                level_heights=np.array([1000.0], dtype=np.float64),
                field_names=["V"],
                parallel=False,
            )

    def test_run_radar_network_3d_marks_qc_fallback_as_false(self):
        from pycwr.interp.RadarInterp import run_radar_network_3d

        class FakePRD(object):
            def __init__(self, radar_id):
                self.radar_id = radar_id
                self.effective_earth_radius = None
                self.nsweeps = 1
                self.scan_info = xr.Dataset(
                    data_vars={
                        "longitude": xr.DataArray(120.0),
                        "latitude": xr.DataArray(30.0),
                    }
                )
                self.fields = [
                    xr.Dataset(
                        data_vars={
                            "dBZ": xr.DataArray(np.ones((2, 2)), dims=("time", "range")),
                        }
                    )
                ]

            def get_vol_data(self, field_name="dBZ", fillvalue=-999.0, range_mode="aligned", max_range_km=None):
                self.vol = (
                    [np.array([0.0, 90.0], dtype=np.float64)],
                    [np.array([1000.0, 2000.0], dtype=np.float64)],
                    np.array([1.0], dtype=np.float64),
                    [np.full((2, 2), 5.0, dtype=np.float64)],
                    100.0,
                    120.0,
                    30.0,
                )

        def fake_read_auto(path, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
            return FakePRD(Path(path).parent.name)

        def fake_get_radar_info(path):
            return 30.0, 120.0, 100.0, 2.8

        with tempfile.TemporaryDirectory() as tmpdir:
            radar_a = Path(tmpdir) / "Z9001"
            radar_a.mkdir()
            (radar_a / "A_20260317070000_test.bin").write_text("", encoding="ascii")

            with mock.patch("pycwr.interp.RadarInterp.read_auto", side_effect=fake_read_auto):
                with mock.patch("pycwr.interp.RadarInterp.get_radar_info", side_effect=fake_get_radar_info):
                    dataset = run_radar_network_3d(
                        target_time="2026-03-17T07:00:00",
                        radar_dirs=[str(radar_a)],
                        lon_min=120.0,
                        lon_max=120.0,
                        lat_min=30.0,
                        lat_max=30.0,
                        lon_res_deg=0.01,
                        lat_res_deg=0.01,
                        level_heights=np.array([1000.0], dtype=np.float64),
                        field_names=["dBZ"],
                        parallel=False,
                        use_qc=True,
                    )

        self.assertEqual(dataset["dBZ"].attrs["qc_applied"], "false")
        self.assertFalse(json.loads(dataset.attrs["qc_applied_by_radar"])["Z9001"]["applied"])

    def test_run_radar_network_3d_uses_product_level_heights_for_vil_et(self):
        from pycwr.interp.RadarInterp import run_radar_network_3d

        class FakePRD(object):
            def __init__(self, radar_id):
                self.radar_id = radar_id
                self.effective_earth_radius = None
                self.scan_info = xr.Dataset(
                    data_vars={
                        "longitude": xr.DataArray(120.0),
                        "latitude": xr.DataArray(30.0),
                    }
                )

            def get_vol_data(self, field_name="dBZ", fillvalue=-999.0, range_mode="aligned", max_range_km=None):
                self.vol = (
                    [
                        np.array([0.0, 90.0], dtype=np.float64),
                        np.array([0.0, 90.0], dtype=np.float64),
                    ],
                    [
                        np.array([1000.0, 2000.0], dtype=np.float64),
                        np.array([1000.0, 2000.0], dtype=np.float64),
                    ],
                    np.array([1.0, 2.0], dtype=np.float64),
                    [
                        np.full((2, 2), 10.0, dtype=np.float64),
                        np.full((2, 2), 12.0, dtype=np.float64),
                    ],
                    100.0,
                    120.0,
                    30.0,
                )

        def fake_read_auto(path, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
            return FakePRD(Path(path).parent.name)

        def fake_get_radar_info(path):
            return 30.0, 120.0, 100.0, 2.8

        def fake_grid_single(prd, field_name, grid_x, grid_y, level_heights, radar_x, radar_y, fillvalue=-999.0,
                             effective_earth_radius=None, range_mode="aligned", max_range_km=None,
                             blind_method="mask", return_metadata=False):
            grid = np.full((len(level_heights),) + grid_x.shape, 7.0, dtype=np.float64)
            if return_metadata:
                return grid, {"actual_max_range_km": 230.0, "range_mode": range_mode}
            return grid

        observed_vil_heights = []
        observed_et_heights = []

        with tempfile.TemporaryDirectory() as tmpdir:
            radar_a = Path(tmpdir) / "Z9001"
            radar_a.mkdir()
            (radar_a / "A_20260317070000_test.bin").write_text("", encoding="ascii")

            with mock.patch("pycwr.interp.RadarInterp.read_auto", side_effect=fake_read_auto):
                with mock.patch("pycwr.interp.RadarInterp.get_radar_info", side_effect=fake_get_radar_info):
                    with mock.patch(
                        "pycwr.interp.RadarInterp.grid_single_radar_to_latlon_3d",
                        side_effect=fake_grid_single,
                    ):
                        with mock.patch(
                            "pycwr.interp.RadarInterp._derive_vil",
                            side_effect=lambda volume, heights, **kwargs: observed_vil_heights.append(
                                tuple(np.asarray(heights, dtype=np.float64))
                            ) or np.zeros(volume.shape[1:], dtype=np.float64),
                        ):
                            with mock.patch(
                                "pycwr.interp.RadarInterp._derive_et",
                                side_effect=lambda volume, heights, **kwargs: (
                                    observed_et_heights.append(tuple(np.asarray(heights, dtype=np.float64)))
                                    or (
                                        np.zeros(volume.shape[1:], dtype=np.float64),
                                        np.zeros(volume.shape[1:], dtype=np.uint8),
                                    )
                                ),
                            ):
                                run_radar_network_3d(
                                    target_time="2026-03-17T07:00:00",
                                    radar_dirs=[str(radar_a)],
                                    lon_min=120.0,
                                    lon_max=120.0,
                                    lat_min=30.0,
                                    lat_max=30.0,
                                    lon_res_deg=0.01,
                                    lat_res_deg=0.01,
                                    level_heights=np.array([1000.0, 3000.0], dtype=np.float64),
                                    product_level_heights=np.array([500.0, 1000.0, 1500.0], dtype=np.float64),
                                    field_names=["dBZ"],
                                    parallel=False,
                                )

        self.assertEqual(observed_vil_heights, [(500.0, 1000.0, 1500.0)])
        self.assertEqual(observed_et_heights, [(500.0, 1000.0, 1500.0)])

    def test_run_radar_network_3d_writes_et_topped_and_vil_reference_metadata(self):
        from pycwr.interp.RadarInterp import run_radar_network_3d

        class FakePRD(object):
            def __init__(self, radar_id):
                self.radar_id = radar_id
                self.effective_earth_radius = None
                self.scan_info = xr.Dataset(
                    data_vars={
                        "longitude": xr.DataArray(120.0),
                        "latitude": xr.DataArray(30.0),
                    }
                )

            def get_vol_data(self, field_name="dBZ", fillvalue=-999.0, range_mode="aligned", max_range_km=None):
                self.vol = (
                    [np.array([0.0], dtype=np.float64)],
                    [np.array([1000.0, 2000.0], dtype=np.float64)],
                    np.array([1.0], dtype=np.float64),
                    [np.full((1, 2), 35.0, dtype=np.float64)],
                    100.0,
                    120.0,
                    30.0,
                )

        def fake_read_auto(path, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
            return FakePRD(Path(path).parent.name)

        def fake_get_radar_info(path):
            return 30.0, 120.0, 100.0, 2.8

        def fake_grid_single(prd, field_name, grid_x, grid_y, level_heights, radar_x, radar_y, fillvalue=-999.0,
                             effective_earth_radius=None, range_mode="aligned", max_range_km=None,
                             blind_method="mask", return_metadata=False):
            values = np.linspace(30.0, 40.0, len(level_heights), dtype=np.float64)[:, np.newaxis, np.newaxis]
            grid = np.broadcast_to(values, (len(level_heights),) + grid_x.shape).copy()
            if return_metadata:
                return grid, {"actual_max_range_km": 230.0, "range_mode": range_mode}
            return grid

        with tempfile.TemporaryDirectory() as tmpdir:
            radar_a = Path(tmpdir) / "Z9001"
            radar_a.mkdir()
            (radar_a / "A_20260317070000_test.bin").write_text("", encoding="ascii")

            with mock.patch("pycwr.interp.RadarInterp.read_auto", side_effect=fake_read_auto):
                with mock.patch("pycwr.interp.RadarInterp.get_radar_info", side_effect=fake_get_radar_info):
                    with mock.patch(
                        "pycwr.interp.RadarInterp.grid_single_radar_to_latlon_3d",
                        side_effect=fake_grid_single,
                    ):
                        dataset = run_radar_network_3d(
                            target_time="2026-03-17T07:00:00",
                            radar_dirs=[str(radar_a)],
                            lon_min=120.0,
                            lon_max=120.0,
                            lat_min=30.0,
                            lat_max=30.0,
                            lon_res_deg=0.01,
                            lat_res_deg=0.01,
                            level_heights=np.array([1000.0, 2000.0], dtype=np.float64),
                            product_level_heights=np.array([1000.0, 2000.0, 3000.0], dtype=np.float64),
                            field_names=["dBZ"],
                            output_products=["VIL", "ET"],
                            composite_method="exp_weighted",
                            parallel=False,
                        )

        self.assertIn("ET_TOPPED", dataset.data_vars)
        self.assertEqual(int(dataset["ET_TOPPED"].values[0, 0]), 1)
        self.assertIn("Greene and Clark", dataset["VIL"].attrs["implementation_note"])
        self.assertIn("Greene and Clark", json.loads(dataset["VIL"].attrs["references"])[0])
        self.assertIn("Lakshmanan", json.loads(dataset["ET"].attrs["references"])[0])

    def test_run_radar_network_3d_uses_reference_plot_for_cr(self):
        from pycwr.interp.RadarInterp import run_radar_network_3d

        class FakePRD(object):
            def __init__(self, radar_id):
                self.radar_id = radar_id
                self.effective_earth_radius = None
                self.scan_info = xr.Dataset(
                    data_vars={
                        "longitude": xr.DataArray(120.0),
                        "latitude": xr.DataArray(30.0),
                    }
                )

            def get_vol_data(self, field_name="dBZ", fillvalue=-999.0, range_mode="aligned", max_range_km=None):
                self.vol = (
                    [
                        np.array([0.0, 90.0], dtype=np.float64),
                        np.array([0.0, 90.0], dtype=np.float64),
                    ],
                    [
                        np.array([1000.0, 2000.0], dtype=np.float64),
                        np.array([1000.0, 2000.0], dtype=np.float64),
                    ],
                    np.array([1.0, 2.0], dtype=np.float64),
                    [
                        np.full((2, 2), 10.0, dtype=np.float64),
                        np.full((2, 2), 12.0, dtype=np.float64),
                    ],
                    100.0,
                    120.0,
                    30.0,
                )

        def fake_read_auto(path, station_lon=None, station_lat=None, station_alt=None, effective_earth_radius=None):
            return FakePRD(Path(path).parent.name)

        def fake_get_radar_info(path):
            return 30.0, 120.0, 100.0, 2.8

        def fake_grid_single(prd, field_name, grid_x, grid_y, level_heights, radar_x, radar_y, fillvalue=-999.0,
                             effective_earth_radius=None, range_mode="aligned", max_range_km=None,
                             blind_method="mask", return_metadata=False):
            grid = np.full((len(level_heights),) + grid_x.shape, 7.0, dtype=np.float64)
            if return_metadata:
                return grid, {"actual_max_range_km": 230.0, "range_mode": range_mode}
            return grid

        with tempfile.TemporaryDirectory() as tmpdir:
            radar_a = Path(tmpdir) / "Z9001"
            radar_a.mkdir()
            (radar_a / "A_20260317070000_test.bin").write_text("", encoding="ascii")

            with mock.patch("pycwr.interp.RadarInterp.read_auto", side_effect=fake_read_auto):
                with mock.patch("pycwr.interp.RadarInterp.get_radar_info", side_effect=fake_get_radar_info):
                    with mock.patch(
                        "pycwr.interp.RadarInterp.grid_single_radar_to_latlon_3d",
                        side_effect=fake_grid_single,
                    ):
                        with mock.patch(
                            "pycwr.interp.RadarInterp._plot_reference_cr_field",
                            return_value="/tmp/reference.png",
                        ) as reference_plot:
                            with mock.patch(
                                "pycwr.interp.RadarInterp._plot_2d_field",
                                return_value="/tmp/simple.png",
                            ) as simple_plot:
                                dataset = run_radar_network_3d(
                                    target_time="2026-03-17T07:00:00",
                                    radar_dirs=[str(radar_a)],
                                    lon_min=120.0,
                                    lon_max=120.0,
                                    lat_min=30.0,
                                    lat_max=30.0,
                                    lon_res_deg=0.01,
                                    lat_res_deg=0.01,
                                    level_heights=np.array([1000.0], dtype=np.float64),
                                    field_names=["dBZ"],
                                    output_products=["CR"],
                                    plot_overview=False,
                                    plot_height_levels=[],
                                    plot_output_dir=tmpdir,
                                    parallel=False,
                                    plot_style="reference",
                                )

        self.assertTrue(reference_plot.called)
        self.assertFalse(simple_plot.called)
        self.assertIn("CR", json.loads(dataset.attrs["plot_files"]))


if __name__ == "__main__":
    unittest.main()
