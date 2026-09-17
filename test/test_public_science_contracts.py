"""Synthetic contracts for public science workflows; no sample files required."""

import unittest

import numpy as np


def make_radar(reflectivity=30.0, native=False):
    from pycwr.core.NRadar import PRD

    angles = np.array([1.0, 4.0, 8.0])
    azimuth_one = np.roll(np.arange(0.0, 360.0, 10.0), 7)
    ranges = np.arange(1.0, 13.0) * 1000.0
    azimuth = np.tile(azimuth_one, angles.size)
    elevation = np.repeat(angles, azimuth_one.size)
    time = np.datetime64("2026-09-17T00:00:00") + np.arange(azimuth.size).astype("timedelta64[s]")
    shape = (azimuth.size, ranges.size)
    velocity = (8.0 * np.sin(np.deg2rad(azimuth)) - 4.0 * np.cos(np.deg2rad(azimuth))) * np.cos(np.deg2rad(elevation))
    fields = {
        "dBZ": np.full(shape, reflectivity), "V": np.repeat(velocity[:, None], ranges.size, axis=1),
        "ZDR": np.full(shape, 0.5), "KDP": np.full(shape, 0.5), "CC": np.full(shape, 0.99),
        "PhiDP": np.broadcast_to(np.arange(ranges.size), shape).copy(),
    }
    extended = {}
    if native:
        native_ranges = np.arange(1.0, 17.0) * 1000.0
        extended["dBZ"] = {}
        for sweep, angle in enumerate(angles):
            extended["dBZ"][sweep] = {
                "range": native_ranges, "azimuth": azimuth_one,
                "elevation": np.full(azimuth_one.size, angle),
                "time": time[sweep * azimuth_one.size:(sweep + 1) * azimuth_one.size],
                "data": np.full((azimuth_one.size, native_ranges.size), reflectivity),
            }
    return PRD(
        fields=fields, scan_type="ppi", time=time, range=ranges, azimuth=azimuth, elevation=elevation,
        latitude=30.0, longitude=120.0, altitude=100.0,
        sweep_start_ray_index=np.arange(angles.size) * azimuth_one.size,
        sweep_end_ray_index=(np.arange(angles.size) + 1) * azimuth_one.size - 1,
        fixed_angle=angles, bins_per_sweep=np.full(angles.size, ranges.size),
        nyquist_velocity=np.full(angles.size, 25.0), frequency=5.6,
        unambiguous_range=np.full(angles.size, 150000.0), nrays=azimuth.size,
        nsweeps=angles.size, sitename="CONTRACT", extended_fields=extended,
    )


class PublicRadarProductContracts(unittest.TestCase):
    def test_native_summary_and_sorted_cache_without_sample_file(self):
        radar = make_radar(native=True)
        summary = radar.summary()
        self.assertEqual(summary["nsweeps"], 3)
        self.assertEqual(summary["sweeps"][0]["native_max_range_m_by_field"]["dBZ"], 16000.0)
        native = radar.get_native_sweep_field(0, "dBZ")
        aligned = radar.get_sweep_field(0, "dBZ", range_mode="aligned")
        self.assertEqual(native.shape, (36, 16))
        self.assertEqual(aligned.shape, (36, 12))
        first = radar.get_sweep_field(0, "dBZ", sort_by_azimuth=True)
        self.assertIs(first, radar.get_sweep_field(0, "dBZ", sort_by_azimuth=True))
        self.assertTrue(np.all(np.diff(first.azimuth) > 0.0))
        self.assertTrue(radar.has_extended_field(0, "dBZ"))
        self.assertEqual(radar.resolve_field_name("reflectivity"), "dBZ")
        with self.assertRaises(ValueError):
            radar.get_sweep_field(0, "dBZ", range_mode="invalid")
        with self.assertRaises(KeyError):
            radar.get_sweep_field(0, "missing")

    def test_xy_product_family_values_coverage_and_metadata(self):
        radar = make_radar()
        axis = np.array([-30000.0, -5000.0, 5000.0, 30000.0])
        levels = np.array([300.0, 500.0, 700.0])
        radar.add_product_CR_xy(axis, axis)
        radar.add_product_CAPPI_xy(axis, axis, 500.0)
        radar.add_product_CAPPI_3d_xy(axis, axis, levels)
        radar.add_product_VIL_xy(axis, axis, levels)
        radar.add_product_ET_xy(axis, axis, levels)
        for name in ("CR_native", "CAPPI_500_native"):
            field = radar.product[name]
            self.assertEqual(field.shape, (4, 4))
            np.testing.assert_allclose(field.values[1:3, 1:3], 30.0)
            self.assertTrue(np.isnan(field.values[0]).all())
            self.assertEqual(field.attrs["units"], "dBZ")
        np.testing.assert_allclose(radar.product.CAPPI_3D_native.values[:, 1:3, 1:3], 30.0)
        np.testing.assert_allclose(radar.product.ET_native.values[1:3, 1:3], levels[-1])
        np.testing.assert_array_equal(radar.product.ET_TOPPED_native.values[1:3, 1:3], 1)
        expected_vil = 3.44e-6 * (10.0 ** 3.0) ** (4.0 / 7.0) * (levels[-1] - levels[0])
        np.testing.assert_allclose(radar.product.VIL_native.values[1:3, 1:3], expected_vil)

    def test_lonlat_product_family_uses_requested_axes(self):
        radar = make_radar()
        lon = np.array([119.96, 120.04])
        lat = np.array([29.96, 30.0, 30.04])
        levels = np.array([300.0, 500.0])
        radar.add_product_CR_lonlat(lon, lat)
        radar.add_product_CAPPI_lonlat(lon, lat, 400.0)
        radar.add_product_VIL_lonlat(lon, lat, levels)
        radar.add_product_ET_lonlat(lon, lat, levels)
        for name in ("CR_geo_native", "CAPPI_geo_400_native", "VIL_geo_native", "ET_geo_native", "ET_TOPPED_geo_native"):
            field = radar.product[name]
            self.assertEqual(field.shape, (2, 3))
            np.testing.assert_array_equal(field.coords[field.dims[0]], lon)
            np.testing.assert_array_equal(field.coords[field.dims[1]], lat)
            self.assertTrue(np.isfinite(field.values).all(), name)
        np.testing.assert_allclose(radar.product.CR_geo_native, 30.0)

    def test_multiradar_public_grid_uses_maximum_and_rejects_empty_network(self):
        from pycwr.core.NRadar import grid_3d_network_xy

        grid = grid_3d_network_xy([make_radar(20.0), make_radar(35.0)], [5000.0], [0.0], [400.0])
        self.assertEqual(grid.network_3d.dims, ("z", "x", "y"))
        np.testing.assert_allclose(grid.network_3d, 35.0)
        with self.assertRaises(ValueError):
            grid_3d_network_xy([], [5000.0], [0.0], [400.0])


class ClassificationAndQCContracts(unittest.TestCase):
    def test_classification_copy_inplace_and_cache_invalidation(self):
        radar = make_radar()
        radar.summary()
        classified = radar.classify_hydrometeors(
            sweeps=[1], profile_height=[0.0, 2000.0], profile_temperature=[20.0, 0.0],
            confidence_field="HCL_CONF", temperature_field="HCL_T",
        )
        self.assertNotIn("HCL", radar.available_fields())
        self.assertNotIn("HCL", classified.fields[0])
        self.assertIn("HCL", classified.summary()["fields"])
        for name in ("HCL", "HCL_CONF", "HCL_T"):
            self.assertEqual(classified.fields[1][name].dims, classified.fields[1].dBZ.dims)
            self.assertTrue(np.isfinite(classified.fields[1][name]).all())
        returned = radar.add_hydrometeor_classification(sweeps=[0])
        self.assertIs(returned, radar)
        self.assertIn("HCL", radar.summary()["fields"])

    def test_classification_failure_paths_are_explicit(self):
        from pycwr.retrieve import classify_hydrometeors, interpolate_temperature_profile

        for kwargs in ({"band": "invalid"}, {"method": "invalid"}, {"ZDR": np.zeros((2, 3))}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                classify_hydrometeors(np.ones((2, 2)), **kwargs)
        with self.assertRaises(ValueError):
            classify_hydrometeors(np.ones((2, 2)))
        with self.assertRaises(ValueError):
            interpolate_temperature_profile([100.0], [0.0], [20.0])

    def test_qc_selected_sweep_copy_and_missing_phase_failure(self):
        radar = make_radar()
        radar.summary()
        corrected = radar.apply_dualpol_qc(sweeps=[1])
        self.assertNotIn("Zc", radar.available_fields())
        self.assertNotIn("Zc", corrected.fields[0])
        self.assertIn("Zc", corrected.summary()["fields"])
        self.assertTrue((corrected.fields[1].Zc >= corrected.fields[1].dBZ).all())
        radar.fields[0] = radar.fields[0].drop_vars(["KDP", "PhiDP"])
        with self.assertRaises(ValueError):
            radar.apply_dualpol_qc(sweeps=[0])

    def test_qc_rejects_invalid_gate_spacing(self):
        from pycwr.qc import correct_attenuation, correct_attenuation_HB, kdp_from_phidp, pia_from_kdp

        data = np.ones((2, 4))
        for spacing in (0.0, -0.1, np.nan, np.inf):
            for func, keyword in ((correct_attenuation, "rscale"), (correct_attenuation_HB, "gate_length"),
                                  (kdp_from_phidp, "dr"), (pia_from_kdp, "dr")):
                with self.subTest(func=func.__name__, spacing=spacing), self.assertRaises(ValueError):
                    func(data, **{keyword: spacing})

    def test_legacy_attenuation_bands_and_single_gate_xarray(self):
        import xarray as xr
        from pycwr.qc import correct_attenuation

        for band in ("S", "C", "X"):
            corrected, quality = correct_attenuation(np.full((2, 4), 30.0), wavelength=band)
            self.assertTrue((corrected >= 30.0).all())
            self.assertTrue(((quality >= 0.0) & (quality <= 1.0)).all())
        single = xr.DataArray([[30.0]], dims=("time", "range"), coords={"range": [1000.0]})
        corrected, quality = correct_attenuation(single)
        np.testing.assert_array_equal(corrected, [[30.0]])
        np.testing.assert_array_equal(quality, [[1.0]])
        with self.assertRaises(ValueError):
            correct_attenuation(single, wavelength="invalid")

    def test_legacy_hid_alias_matches_public_classifier(self):
        from pycwr.retrieve.HID import classify_hydrometeors, fhc_HCL, fhc_hcl

        kwargs = {"dBZ": np.full((2, 2), 30.0), "ZDR": np.full((2, 2), 0.5)}
        np.testing.assert_array_equal(fhc_HCL(**kwargs), classify_hydrometeors(**kwargs))
        np.testing.assert_array_equal(fhc_hcl(**kwargs), classify_hydrometeors(**kwargs))


class WindAndGeometryBoundaryContracts(unittest.TestCase):
    def test_antenna_vector_edges_are_interpolated_once_for_arrays_and_xarray(self):
        import xarray as xr
        from pycwr.core.transforms import antenna_vectors_to_cartesian, antenna_vectors_to_cartesian_cwr

        ranges = np.array([1000.0, 2000.0, 3000.0])
        azimuth = np.array([0.0, 90.0, 180.0, 270.0])
        elevation = np.ones(azimuth.shape)
        expected = antenna_vectors_to_cartesian_cwr(ranges, azimuth, elevation, edges=True)
        for convert in (np.asarray, xr.DataArray):
            for transform in (antenna_vectors_to_cartesian, antenna_vectors_to_cartesian_cwr):
                with self.subTest(input_type=convert.__name__, transform=transform.__name__):
                    result = transform(convert(ranges), convert(azimuth), convert(elevation), edges=True)
                    for coordinate, reference in zip(result, expected):
                        self.assertEqual(coordinate.shape, (5, 4))
                        np.testing.assert_allclose(coordinate, reference)

    def test_lonlat_wind_product_persists_retrieval_and_quality_fields(self):
        radar = make_radar()
        volume = radar.add_product_WIND_VOLUME_lonlat(
            [120.03], [30.03], [400.0], sweeps=[0, 1, 2],
            az_num=9, bin_num=3, azimuth_step=6, range_step=3,
            horizontal_radius_m=8000.0, horizontal_min_neighbors=1, vertical_tolerance_m=1000.0,
        )
        self.assertTrue(np.isfinite(volume.u).all())
        np.testing.assert_allclose(volume.u, 8.0, atol=1e-5)
        np.testing.assert_allclose(volume.v, -4.0, atol=1e-5)
        for suffix, name in (("u", "u"), ("v", "v"), ("quality_score", "quality_score")):
            product = radar.product["WIND_VOLUME_geo_" + suffix]
            self.assertEqual(product.dims, ("z_wind_geo", "lon_wind", "lat_wind"))
            np.testing.assert_array_equal(product, volume[name])

    def test_generic_geographic_projection_roundtrip_and_vector_edges(self):
        from pycwr.core.transforms import geographic_to_cartesian, cartesian_to_geographic, cartesian_vectors_to_geographic

        for projection in ("pyart_aeqd", "aeqd"):
            params = {"proj": projection, "lon_0": 120.0, "lat_0": 30.0, "R": 6370997.0}
            lon = np.array([119.9, 120.0, 120.1])
            lat = np.array([29.9, 30.0, 30.1])
            x, y = geographic_to_cartesian(lon, lat, params)
            lon2, lat2 = cartesian_to_geographic(x, y, params)
            np.testing.assert_allclose(lon2, lon, atol=1e-10)
            np.testing.assert_allclose(lat2, lat, atol=1e-10)
            lon_edges, lat_edges = cartesian_vectors_to_geographic(np.array([-1000.0, 1000.0]), np.array([0.0]), params, edges=True)
            self.assertEqual(lon_edges.shape, (1, 3))
            np.testing.assert_allclose(lon_edges[0, 1], 120.0)
            np.testing.assert_allclose(lat_edges[0, 1], 30.0)

    def test_vad_all_missing_field_and_empty_range_do_not_invent_wind(self):
        radar = make_radar()
        for field in radar.fields:
            field.V.values[:] = np.nan
        vad = radar.retrieve_vad(sweeps=[0])
        self.assertTrue(np.isnan(vad.u).all())
        self.assertTrue((vad.valid_count == 0).all())
        empty = radar.retrieve_vad(sweeps=[0], max_range_km=0.0)
        self.assertEqual(empty.sizes["gate"], 0)

    def test_wind_retrieval_rejects_invalid_public_options(self):
        radar = make_radar()
        for kwargs in ({"gate_step": 0}, {"sweeps": []}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                radar.retrieve_vad(**kwargs)
        with self.assertRaises(KeyError):
            radar.retrieve_vad(field_name="missing")
        with self.assertRaises(ValueError):
            radar.retrieve_vvp(0, az_num=4)
        for levels in ([], [200.0, 100.0], [np.nan]):
            with self.subTest(levels=levels), self.assertRaises(ValueError):
                radar.retrieve_wind_volume_xy([0.0], [0.0], levels)

    def test_effective_radius_rejects_nonfinite_values(self):
        from pycwr.core import RadarGridC
        from pycwr.core.transforms import resolve_effective_earth_radius

        for radius in (0.0, -1.0, np.nan, np.inf):
            with self.subTest(radius=radius), self.assertRaises(ValueError):
                resolve_effective_earth_radius(radius)
            with self.subTest(backend="RadarGridC", radius=radius), self.assertRaises(ValueError):
                RadarGridC.antenna_to_cartesian(1000.0, 90.0, 1.0, 100.0, effective_earth_radius=radius)

    def test_temperature_profile_sorts_altitudes_and_handles_agl(self):
        from pycwr.retrieve import interpolate_temperature_profile

        result = interpolate_temperature_profile(
            np.array([[100.0, 600.0], [1100.0, np.nan]]), [1000.0, 0.0], [10.0, 20.0],
            radar_altitude=100.0, height_reference="agl",
        )
        np.testing.assert_allclose(result, [[20.0, 15.0], [10.0, np.nan]], equal_nan=True)


class NetworkCompositionContracts(unittest.TestCase):
    def test_nearest_and_max_ignore_missing_radar_values(self):
        from pycwr.interp.RadarInterp import compose_network_volume

        volumes = [np.array([[[10.0, -999.0]]]), np.array([[[30.0, 40.0]]])]
        x = np.array([[0.0, 1.0]])
        y = np.zeros_like(x)
        for method, expected in (("nearest", [[[10.0, 40.0]]]), ("max", [[[30.0, 40.0]]])):
            with self.subTest(method=method):
                result, count = compose_network_volume(volumes, [[0.0, 0.0], [1000.0, 0.0]], x, y, method=method)
                np.testing.assert_array_equal(result, expected)
                np.testing.assert_array_equal(count, [[[2, 1]]])

    def test_network_rejects_empty_or_inconsistent_volumes(self):
        from pycwr.interp.RadarInterp import compose_network_volume

        grid = np.zeros((1, 1))
        with self.assertRaises(ValueError):
            compose_network_volume([], [], grid, grid)
        with self.assertRaises(ValueError):
            compose_network_volume([np.ones((1, 1, 1)), np.ones((2, 1, 1))], [[0, 0], [0, 0]], grid, grid)


class RemainingSciencePublicContracts(unittest.TestCase):
    def test_python_ppi_interpolation_recovers_plane_and_handles_missing_rows(self):
        from pycwr.core import RadarGrid

        # A plane increasing by 9 across azimuth and 4 across range.
        value = RadarGrid.interp_ppi(30.0, 1250.0, 0.0, 90.0, 1000.0, 2000.0, 3.0, 7.0, 12.0, 16.0)
        self.assertAlmostEqual(value, 7.0)
        fallback = RadarGrid.interp_ppi(30.0, 1250.0, 0.0, 90.0, 1000.0, 2000.0, -999.0, -999.0, 12.0, 16.0)
        self.assertAlmostEqual(fallback, 13.0)
        self.assertEqual(RadarGrid.interp_ppi(30.0, 1250.0, 0.0, 90.0, 1000.0, 2000.0,
                                             -999.0, -999.0, -999.0, -999.0), -999.0)
        self.assertEqual(RadarGrid.interp_ppi(30.0, 1000.0, 0.0, 90.0, 1000.0, 1000.0,
                                             3.0, 7.0, 12.0, 16.0), -999.0)

    def test_python_grid_range_support_blind_fill_and_composite_maximum(self):
        from pycwr.core import RadarGrid
        from pycwr.core.transforms import antenna_to_cartesian_cwr

        azimuth = np.array([0.0, 90.0, 180.0, 270.0])
        ranges = np.array([1000.0, 2000.0])
        x, y, _ = antenna_to_cartesian_cwr(np.array([500.0, 1500.0, 3000.0]), 0.0, 1.0, 100.0)
        x, y = x[None, :], y[None, :]
        lower, upper = np.full((4, 2), 10.0), np.full((4, 2), 20.0)
        grid = RadarGrid.ppi_to_grid(azimuth, ranges, 1.0, lower, 100.0, x, y)
        np.testing.assert_allclose(grid, [[-999.0, 10.0, -999.0]])
        blind_filled = RadarGrid.ppi_to_grid(azimuth, ranges, 1.0, lower, 100.0, x, y, blind_method="nearest_gate")
        np.testing.assert_allclose(blind_filled, [[10.0, 10.0, -999.0]])
        composite = RadarGrid.get_CR_xy([azimuth, azimuth], [ranges, ranges], np.array([1.0, 2.0]),
                                       [lower, upper], 100.0, x, y)
        np.testing.assert_allclose(composite, [[-999.0, 20.0, -999.0]])
        missing = RadarGrid.get_CR_xy([azimuth], [ranges], np.array([1.0]),
                                     [np.full((4, 2), -999.0)], 100.0, x, y)
        np.testing.assert_array_equal(missing, np.full((1, 3), -999.0))

    def test_ordered_az_view_preserves_data_then_inplace_changes_dimension(self):
        from pycwr.core.NRadar import AzimuthSortedPRD

        radar = make_radar()
        original = radar.fields[0].copy(deep=True)
        expected = original.sortby("azimuth").swap_dims({"time": "azimuth"})
        view = radar.ordered_az()
        self.assertIsInstance(view, AzimuthSortedPRD)
        self.assertIs(view, radar.ordered_az())
        self.assertEqual(view.effective_earth_radius, radar.effective_earth_radius)
        np.testing.assert_array_equal(radar.fields[0].azimuth, original.azimuth)
        np.testing.assert_array_equal(view.fields[0].V, expected.V)
        np.testing.assert_array_equal(view.fields[0].time, expected.time)
        self.assertIsNone(radar.ordered_az(inplace=True))
        self.assertEqual(radar.fields[0].V.dims, ("azimuth", "range"))
        np.testing.assert_array_equal(radar.fields[0].V, expected.V)
        self.assertIsNot(view, radar.ordered_az())
        empty = AzimuthSortedPRD()
        self.assertEqual(empty.fields, [])
        self.assertEqual(len(empty.product.data_vars), 0)

    def test_hydrometeor_labels_match_one_based_class_ids(self):
        from pycwr.retrieve import available_hydrometeor_classes, hydrometeor_class_name

        classes = available_hydrometeor_classes()
        self.assertEqual(len(classes), 10)
        self.assertEqual(len(set(classes)), 10)
        self.assertEqual(hydrometeor_class_name(1), "Drizzle")
        self.assertEqual(hydrometeor_class_name(9), "Hail")
        self.assertEqual(hydrometeor_class_name(10), "Big Drops")
        for class_id, name in enumerate(classes, start=1):
            self.assertEqual(hydrometeor_class_name(class_id), name)
        for invalid in (0, 11):
            with self.subTest(class_id=invalid), self.assertRaises(ValueError):
                hydrometeor_class_name(invalid)

    def test_velocity_selection_prefers_corrected_but_honors_explicit_field(self):
        from pycwr.retrieve import select_velocity_field

        radar = make_radar()
        self.assertEqual(select_velocity_field(radar, 0), "V")
        radar.fields[0]["Vc"] = radar.fields[0].V.copy()
        self.assertEqual(select_velocity_field(radar, 0), "Vc")
        self.assertEqual(select_velocity_field(radar, 0, field_name="V"), "V")
        with self.assertRaises(KeyError):
            select_velocity_field(radar, 0, field_name="missing")
        radar.fields[0] = radar.fields[0].drop_vars(["V", "Vc"])
        with self.assertRaises(KeyError):
            select_velocity_field(radar, 0)


if __name__ == "__main__":
    unittest.main()
