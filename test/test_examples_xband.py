"""X-band usage examples for attenuation and dual-pol QC."""

import unittest

import numpy as np


class XBandQcExamples(unittest.TestCase):
    def test_resolve_x_band_coefficients_example(self):
        from pycwr.qc import resolve_kdp_coefficients

        gamma, beta = resolve_kdp_coefficients(band="X")

        self.assertAlmostEqual(gamma, 0.28)
        self.assertAlmostEqual(beta, 0.04)

    def test_run_dualpol_qc_x_band_example(self):
        from pycwr.qc import run_dualpol_qc

        ref = np.tile(np.linspace(20.0, 38.0, 10), (3, 1))
        zdr = np.tile(np.linspace(0.3, 1.1, 10), (3, 1))
        phidp = np.tile(np.linspace(0.0, 24.0, 10), (3, 1))
        rhohv = np.ones_like(ref) * 0.985
        snr = np.ones_like(ref) * 18.0

        result = run_dualpol_qc(
            ref,
            zdr=zdr,
            phidp=phidp,
            rhohv=rhohv,
            snr=snr,
            dr=0.05,
            band="X",
            use_existing_kdp=False,
        )

        self.assertIn("ref_corrected", result)
        self.assertIn("kdp_used", result)
        self.assertEqual(result["ref_corrected"].shape, ref.shape)
        self.assertTrue(np.all(np.isnan(result["ref_corrected"]) | (result["ref_corrected"] >= ref)))

    def test_apply_dualpol_qc_x_band_on_prd_example(self):
        from pycwr.core.NRadar import PRD

        nrays = 3
        nbins = 10
        prd = PRD(
            fields={
                "dBZ": np.tile(np.linspace(20.0, 38.0, nbins), (nrays, 1)),
                "ZDR": np.tile(np.linspace(0.3, 1.1, nbins), (nrays, 1)),
                "PhiDP": np.tile(np.linspace(0.0, 24.0, nbins), (nrays, 1)),
                "CC": np.ones((nrays, nbins)) * 0.985,
                "SNRH": np.ones((nrays, nbins)) * 18.0,
            },
            scan_type="ppi",
            time=np.arange(nrays),
            range=np.arange(nbins, dtype=float) * 50.0,
            azimuth=np.linspace(0.0, 240.0, nrays),
            elevation=np.ones(nrays) * 1.0,
            latitude=30.0,
            longitude=120.0,
            altitude=100.0,
            sweep_start_ray_index=np.array([0]),
            sweep_end_ray_index=np.array([nrays - 1]),
            fixed_angle=np.array([1.0]),
            bins_per_sweep=np.array([nbins]),
            nyquist_velocity=np.array([15.0]),
            frequency=9.4,
            unambiguous_range=np.array([150000.0]),
            nrays=nrays,
            nsweeps=1,
            sitename="X-BAND",
        )

        processed = prd.apply_dualpol_qc(inplace=False, band="X", use_existing_kdp=False)

        self.assertIn("Zc", processed.fields[0].data_vars)
        self.assertIn("KDPc", processed.fields[0].data_vars)
        self.assertIn("CLEAR_AIR_MASK", processed.fields[0].data_vars)
        self.assertIn("QC_MASK", processed.fields[0].data_vars)
        self.assertNotIn("Zc", prd.fields[0].data_vars)


if __name__ == "__main__":
    unittest.main()
