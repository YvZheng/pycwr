"""Public numerical API regressions for coordinate shapes and missing data."""

import unittest

import numpy as np


class ScatteredInterpolationTests(unittest.TestCase):
    def _interpolators(self):
        from pycwr.interp.RadarInterp import radar_interp2d, radar_interp2d_var

        return ((radar_interp2d, {"around_r": 100.0}), (radar_interp2d_var, {}))

    def test_coordinate_matrix_and_tuple_grid_agree(self):
        points = np.array([[0.0, 0.0], [10.0, 0.0], [20.0, 0.0]])
        target = np.array([[2.0, 0.0], [5.0, 0.0], [8.0, 0.0]])
        for interp, kwargs in self._interpolators():
            with self.subTest(interp=interp.__name__):
                matrix = interp(points, [10.0, 20.0, 30.0], target, **kwargs)
                grid = interp(tuple(points.T), [10.0, 20.0, 30.0], tuple(target.T), **kwargs)
                self.assertEqual(matrix.shape, (3,))
                np.testing.assert_allclose(matrix, grid)

    def test_scalar_target_and_broadcast_grid(self):
        for interp, kwargs in self._interpolators():
            with self.subTest(interp=interp.__name__):
                scalar = interp([[0.0, 0.0]], [10], (0.0, 0.0), **kwargs)
                grid = interp([[0.0, 0.0]], [10], (np.zeros((2, 1)), np.zeros((1, 3))), **kwargs)
                self.assertEqual(scalar.shape, ())
                self.assertEqual(float(scalar), 10.0)
                np.testing.assert_array_equal(grid, np.full((2, 3), 10.0))

    def test_missing_neighbors_do_not_poison_valid_data(self):
        points = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [np.nan, 0.0]])
        values = np.ma.array([10.0, np.nan, 999.0, 30.0], mask=[False, False, True, False])
        for interp, kwargs in self._interpolators():
            with self.subTest(interp=interp.__name__):
                result = interp(points, values, np.array([[0.0, 0.0]]), **kwargs)
                np.testing.assert_array_equal(result, [10.0])

    def test_all_missing_sources_use_fill_value(self):
        for interp, kwargs in self._interpolators():
            with self.subTest(interp=interp.__name__):
                result = interp([[0.0, 0.0]], [np.nan], np.array([[0.0, 0.0]]), fill_value=-999.0, **kwargs)
                np.testing.assert_array_equal(result, [-999.0])

    def test_complex_data_are_preserved(self):
        for interp, kwargs in self._interpolators():
            with self.subTest(interp=interp.__name__):
                result = interp([[0.0, 0.0]], [1.0 + 2.0j], np.array([[0.0, 0.0]]), **kwargs)
                np.testing.assert_array_equal(result, [1.0 + 2.0j])

    def test_cressman_zero_weight_uses_requested_fill(self):
        from pycwr.interp.RadarInterp import radar_interp2d

        with np.errstate(all="raise"):
            result = radar_interp2d([[1.0, 0.0]], [10.0], np.array([[0.0, 0.0]]),
                                    around_r=2.0, influence_radius=1.0,
                                    method="cressman", fill_value=-999.0)
        np.testing.assert_array_equal(result, [-999.0])

    def test_cressman_excludes_points_outside_influence_radius(self):
        from pycwr.interp.RadarInterp import radar_interp2d

        result = radar_interp2d([[0.0, 0.0], [2.0, 0.0]], [10.0, 100.0],
                                np.array([[0.0, 0.0]]), around_r=3.0,
                                influence_radius=1.0, method="cressman")
        np.testing.assert_array_equal(result, [10.0])

    def test_invalid_method_is_rejected_without_neighbors(self):
        for interp, kwargs in self._interpolators():
            with self.subTest(interp=interp.__name__), self.assertRaises(ValueError):
                interp([[0.0, 0.0]], [np.nan], np.array([[0.0, 0.0]]), method="invalid", **kwargs)


class MissingProductDataTests(unittest.TestCase):
    def test_echo_top_does_not_interpolate_into_missing_gate(self):
        from pycwr.core.RadarProduct import derive_et

        volume = np.array([[30.0, 30.0, 30.0], [-999.0, np.nan, 10.0]])[:, None, :]
        et = derive_et(volume, [1000.0, 2000.0])
        np.testing.assert_allclose(et, [[1000.0, 1000.0, 1600.0]])

    def test_products_treat_masked_cells_as_missing(self):
        from pycwr.core.RadarProduct import derive_cr, derive_et, derive_vil

        volume = np.ma.array([[[30.0, 99.0]], [[99.0, 99.0]]],
                             mask=[[[False, True]], [[True, True]]])
        expected = volume.filled(np.nan)
        np.testing.assert_allclose(derive_cr(volume), derive_cr(expected), equal_nan=True)
        for derive in (derive_et, derive_vil):
            with self.subTest(derive=derive.__name__):
                np.testing.assert_allclose(derive(volume, [1000.0, 2000.0]),
                                           derive(expected, [1000.0, 2000.0]), equal_nan=True)


class MaskedRetrievalTests(unittest.TestCase):
    def test_hid_does_not_classify_masked_reflectivity(self):
        from pycwr.retrieve.HID import classify_hydrometeors

        result = classify_hydrometeors(np.ma.array([30.0, 30.0], mask=[False, True]), ZDR=[1.0, 1.0])
        self.assertTrue(np.isfinite(result[0]))
        self.assertTrue(np.isnan(result[1]))

    def test_vad_counts_only_unmasked_velocity_samples(self):
        from pycwr.retrieve.WindField import fit_vad_ring

        azimuth = np.arange(0.0, 360.0, 10.0)
        velocity = np.ma.array(10.0 * np.sin(np.deg2rad(azimuth)), mask=np.arange(36) < 30)
        result = fit_vad_ring(azimuth, 0.0, velocity)
        self.assertEqual(result["valid_count"], 6)
        self.assertTrue(np.isnan(result["u"]))

    def test_vad_sample_threshold_includes_geometry_and_weight_validity(self):
        from pycwr.retrieve.WindField import fit_vad_ring

        azimuth = np.arange(0.0, 360.0, 10.0)
        velocity = 10.0 * np.sin(np.deg2rad(azimuth))
        weights = np.zeros(azimuth.size)
        weights[::12] = 1.0
        elevation = np.full(azimuth.size, np.nan)
        elevation[::12] = 0.0
        for kwargs in ({"elevation": 0.0, "weights": weights}, {"elevation": elevation}):
            with self.subTest(kwargs=kwargs):
                result = fit_vad_ring(azimuth, radial_velocity=velocity, **kwargs)
                self.assertEqual(result["valid_count"], 3)
                self.assertTrue(np.isnan(result["u"]))

    def test_vad_supports_scalar_weights(self):
        from pycwr.retrieve.WindField import fit_vad_ring

        azimuth = np.arange(0.0, 360.0, 10.0)
        velocity = 10.0 * np.sin(np.deg2rad(azimuth))
        result = fit_vad_ring(azimuth, 0.0, velocity, weights=1.0)
        self.assertAlmostEqual(result["u"], 10.0)

    def test_attenuation_resets_at_masked_gate(self):
        from pycwr.qc import correct_attenuation_kdp, pia_from_kdp

        masked = np.ma.array([[1.0, 999.0, 1.0]], mask=[[False, True, False]])
        np.testing.assert_allclose(pia_from_kdp(masked, dr=0.1), [[0.016, np.nan, 0.016]], equal_nan=True)
        ref = np.ma.array([[20.0, 999.0, 20.0]], mask=[[False, True, False]])
        corrected = correct_attenuation_kdp(ref, np.ones((1, 3)), dr=0.1)[0]
        np.testing.assert_allclose(corrected, [[20.016, np.nan, 20.016]], equal_nan=True)

    def test_phase_smoothing_does_not_use_masked_underlying_values(self):
        from pycwr.qc import smooth_phidp

        phase = np.ma.array([[1.0, 999.0, 3.0]], mask=[[False, True, False]])
        result = smooth_phidp(phase, median_window=3, fit_window=3)
        expected = smooth_phidp(phase.filled(np.nan), median_window=3, fit_window=3)
        np.testing.assert_allclose(result, expected, equal_nan=True)

    def test_qc_pipeline_respects_masked_reflectivity_and_kdp(self):
        from pycwr.qc import run_dualpol_qc

        ref = np.ma.array(np.full((4, 12), 30.0), mask=False)
        kdp = np.ma.array(np.ones((4, 12)), mask=False)
        ref.mask[0, 1] = True
        kdp.mask[1, 1] = True
        masked = run_dualpol_qc(ref, kdp=kdp)
        expected = run_dualpol_qc(ref.filled(np.nan), kdp=kdp.filled(np.nan))
        for field in ("ref_corrected", "pia", "kdp_used", "qc_mask"):
            np.testing.assert_allclose(masked[field], expected[field], equal_nan=True)


if __name__ == "__main__":
    unittest.main()
