"""QC examples and primitive regression tests."""

import unittest

import numpy as np


class AttenuationTests(unittest.TestCase):
    def test_pia_from_kdp_is_monotonic_and_clips_negative_kdp(self):
        from pycwr.qc import pia_from_kdp

        kdp = np.array([[0.5, 1.0, -0.2, 2.0, np.nan]])
        pia = pia_from_kdp(kdp, dr=0.1, gamma=0.08)

        self.assertTrue(np.all(np.diff(pia[0, :4]) >= -1e-12))
        self.assertTrue(np.isnan(pia[0, 4]))

    def test_correct_attenuation_kdp_increases_reflectivity(self):
        from pycwr.qc import correct_attenuation_kdp

        ref = np.array([[20.0, 22.0, 24.0, 26.0]])
        zdr = np.array([[0.5, 0.6, 0.7, 0.8]])
        kdp = np.array([[0.3, 0.4, 0.5, 0.6]])

        zc, pia, zdrc, pia_zdr = correct_attenuation_kdp(ref, kdp, dr=0.1, zdr=zdr)

        self.assertEqual(zc.shape, ref.shape)
        self.assertEqual(pia.shape, ref.shape)
        self.assertTrue(np.all(zc >= ref))
        self.assertTrue(np.all(zdrc >= zdr))
        self.assertTrue(np.all(np.diff(pia[0]) >= -1e-12))
        self.assertTrue(np.all(np.diff(pia_zdr[0]) >= -1e-12))

    def test_correct_attenuation_hb_supports_leading_dimensions(self):
        from pycwr.qc import correct_attenuation_HB

        ref = np.array(
            [
                [[20.0, 21.0, 22.0], [30.0, 31.0, 32.0]],
                [[15.0, 16.0, np.nan], [25.0, 26.0, 27.0]],
            ]
        )
        zc, pia = correct_attenuation_HB(ref, gate_length=0.075)

        self.assertEqual(zc.shape, ref.shape)
        self.assertEqual(pia.shape, ref.shape)
        self.assertTrue(np.isnan(zc[1, 0, 2]))
        self.assertTrue(np.isnan(pia[1, 0, 2]))


class PolarimetricPrimitiveTests(unittest.TestCase):
    def test_smooth_phidp_removes_spike_and_enforces_monotonicity(self):
        from pycwr.qc import smooth_phidp

        phidp = np.array([[0.0, 1.0, 40.0, 3.0, 4.0, 5.0, 6.0]])
        smoothed = smooth_phidp(phidp, median_window=3, fit_window=3, enforce_monotonic=True)

        self.assertLess(smoothed[0, 2], 10.0)
        self.assertTrue(np.all(np.diff(smoothed[0][np.isfinite(smoothed[0])]) >= -1e-12))

    def test_kdp_from_phidp_recovers_linear_slope(self):
        from pycwr.qc import kdp_from_phidp

        dr = 0.1
        kdp_true = 1.0
        gates = np.arange(11)
        phidp = (2.0 * kdp_true * dr) * gates
        kdp = kdp_from_phidp(phidp[np.newaxis, :], dr=dr, fit_window=5)

        self.assertTrue(np.allclose(kdp[0, 2:-2], kdp_true, atol=1e-6))

    def test_phidp_texture_highlights_spikes(self):
        from pycwr.qc import phidp_texture

        phidp = np.array([[1.0, 1.0, 1.0, 10.0, 1.0, 1.0, 1.0]])
        texture = phidp_texture(phidp, window=3)

        self.assertGreater(texture[0, 3], texture[0, 0])

    def test_build_meteo_mask_combines_classic_thresholds(self):
        from pycwr.qc import build_meteo_mask

        ref = np.array([[10.0, 10.0, -5.0]])
        rhohv = np.array([[0.95, 0.60, 0.95]])
        texture = np.array([[5.0, 5.0, 25.0]])
        snr = np.array([[10.0, 10.0, 10.0]])

        mask = build_meteo_mask(ref, rhohv=rhohv, phidp_texture=texture, snr=snr)

        np.testing.assert_array_equal(mask, np.array([[True, False, False]]))

    def test_build_clear_air_mask_identifies_weak_nonprecip_echo(self):
        from pycwr.qc import build_clear_air_mask

        ref = np.array([[8.0, 8.0, 25.0]])
        rhohv = np.array([[0.92, 0.99, 0.99]])
        texture = np.array([[4.0, 4.0, 4.0]])
        snr = np.array([[10.0, 24.0, 24.0]])

        mask = build_clear_air_mask(ref, rhohv=rhohv, phidp_texture=texture, snr=snr)

        np.testing.assert_array_equal(mask, np.array([[True, False, False]]))

    def test_despeckle_mask_removes_small_isolated_echoes(self):
        from pycwr.qc import despeckle_mask

        mask = np.zeros((6, 6), dtype=bool)
        mask[2:5, 2:5] = True
        mask[0, 0] = True

        despeckled = despeckle_mask(mask, min_size=4)

        self.assertFalse(despeckled[0, 0])
        self.assertTrue(np.all(despeckled[2:5, 2:5]))


class PipelineTests(unittest.TestCase):
    def _build_synthetic_prd(self):
        from pycwr.core.NRadar import PRD

        nrays = 4
        nbins = 12
        fields = {
            "dBZ": np.tile(np.linspace(15.0, 35.0, nbins), (nrays, 1)),
            "ZDR": np.tile(np.linspace(0.2, 1.0, nbins), (nrays, 1)),
            "PhiDP": np.tile(np.linspace(0.0, 20.0, nbins), (nrays, 1)),
            "KDP": np.tile(np.linspace(0.3, 1.2, nbins), (nrays, 1)),
            "CC": np.ones((nrays, nbins)) * 0.98,
            "SNRH": np.ones((nrays, nbins)) * 20.0,
        }
        return PRD(
            fields=fields,
            scan_type="ppi",
            time=np.array(
                [np.datetime64("2026-03-20T00:00:00") + np.timedelta64(i, "s") for i in range(nrays)]
            ),
            range=np.arange(nbins, dtype=float) * 75.0,
            azimuth=np.linspace(0.0, 270.0, nrays),
            elevation=np.ones(nrays) * 1.0,
            latitude=30.0,
            longitude=120.0,
            altitude=100.0,
            sweep_start_ray_index=np.array([0]),
            sweep_end_ray_index=np.array([nrays - 1]),
            fixed_angle=np.array([1.0]),
            bins_per_sweep=np.array([nbins]),
            nyquist_velocity=np.array([15.0]),
            frequency=5.6,
            unambiguous_range=np.array([150000.0]),
            nrays=nrays,
            nsweeps=1,
            sitename="TEST",
        )

    def test_run_dualpol_qc_returns_expected_products(self):
        from pycwr.qc import run_dualpol_qc

        ref = np.tile(np.linspace(15.0, 35.0, 12), (4, 1))
        zdr = np.tile(np.linspace(0.2, 1.0, 12), (4, 1))
        phidp = np.tile(np.linspace(0.0, 20.0, 12), (4, 1))
        kdp = np.tile(np.linspace(0.3, 1.2, 12), (4, 1))
        rhohv = np.ones_like(ref) * 0.98
        snr = np.ones_like(ref) * 20.0

        result = run_dualpol_qc(ref, zdr=zdr, phidp=phidp, kdp=kdp, rhohv=rhohv, snr=snr, dr=0.075)

        self.assertEqual(set(result.keys()), {
            "clear_air_mask",
            "ref_corrected",
            "zdr_corrected",
            "pia",
            "pia_zdr",
            "phidp_smooth",
            "kdp_used",
            "phidp_texture",
            "meteo_mask",
            "qc_mask",
            "qc_reference_notes",
        })
        self.assertTrue(np.all(np.isnan(result["ref_corrected"]) | (result["ref_corrected"] >= ref)))

    def test_run_dualpol_qc_can_mask_clear_air_echoes(self):
        from pycwr.qc import run_dualpol_qc

        ref = np.tile(np.array([8.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 35.0, 36.0, 37.0]), (3, 1))
        zdr = np.tile(np.linspace(0.2, 1.0, 10), (3, 1))
        phidp = np.tile(np.linspace(0.0, 18.0, 10), (3, 1))
        rhohv = np.tile(np.array([0.92, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98]), (3, 1))
        snr = np.tile(np.array([10.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0]), (3, 1))

        labeled = run_dualpol_qc(
            ref,
            zdr=zdr,
            phidp=phidp,
            rhohv=rhohv,
            snr=snr,
            dr=0.075,
            use_existing_kdp=False,
            clear_air_mode="label",
        )
        masked = run_dualpol_qc(
            ref,
            zdr=zdr,
            phidp=phidp,
            rhohv=rhohv,
            snr=snr,
            dr=0.075,
            use_existing_kdp=False,
            clear_air_mode="mask",
        )

        self.assertTrue(labeled["clear_air_mask"][0, 0])
        self.assertTrue(labeled["qc_mask"][1, 1])
        self.assertFalse(masked["qc_mask"][0, 0])
        self.assertTrue(np.isnan(masked["ref_corrected"][0, 0]))

    def test_apply_dualpol_qc_adds_fields_without_mutating_original_when_copying(self):
        from pycwr.qc import apply_dualpol_qc

        prd = self._build_synthetic_prd()
        processed = apply_dualpol_qc(prd, inplace=False)

        self.assertNotIn("Zc", prd.fields[0].data_vars)
        for field_name in (
            "Zc",
            "ZDRc",
            "PhiDP_smooth",
            "KDPc",
            "PhiDP_texture",
            "METEO_MASK",
            "CLEAR_AIR_MASK",
            "QC_MASK",
            "PIA",
            "PIA_ZDR",
        ):
            self.assertIn(field_name, processed.fields[0].data_vars)

    def test_apply_dualpol_qc_inplace_mutates_target(self):
        from pycwr.qc import apply_dualpol_qc

        prd = self._build_synthetic_prd()
        returned = apply_dualpol_qc(prd, inplace=True)

        self.assertIs(returned, prd)
        self.assertIn("Zc", prd.fields[0].data_vars)

    def test_prd_apply_dualpol_qc_method_writes_results(self):
        prd = self._build_synthetic_prd()

        processed = prd.apply_dualpol_qc(inplace=False)

        self.assertIn("Zc", processed.fields[0].data_vars)
        self.assertIn("QC_MASK", processed.fields[0].data_vars)
        self.assertNotIn("Zc", prd.fields[0].data_vars)

    def test_qc_fields_export_to_pyart(self):
        prd = self._build_synthetic_prd()
        processed = prd.apply_dualpol_qc(inplace=False)

        radar = processed.to_pyart_radar(use_external=False, force_rebuild=True)

        self.assertIn("corrected_reflectivity", radar.fields)
        self.assertIn("corrected_differential_reflectivity", radar.fields)

    def test_apply_dualpol_qc_requires_dualpol_phase_information(self):
        from pycwr.qc import apply_dualpol_qc

        prd = self._build_synthetic_prd()
        prd.fields[0] = prd.fields[0].drop_vars(["KDP", "PhiDP"])

        with self.assertRaises(ValueError):
            apply_dualpol_qc(prd, inplace=False)


if __name__ == "__main__":
    unittest.main()
