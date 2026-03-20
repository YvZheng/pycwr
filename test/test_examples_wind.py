"""Wind-retrieval examples and regression tests."""

import json
import unittest
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")


class WindRetrievalTests(unittest.TestCase):
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

    @staticmethod
    def _build_uniform_wind_prd(u=8.0, v=-4.0, include_vc=False):
        from pycwr.core.NRadar import PRD

        sweeps = np.array([1.0, 2.0], dtype=np.float64)
        nrays_per_sweep = 72
        nbins = 24
        gate_range = (np.arange(nbins, dtype=np.float64) + 1.0) * 250.0
        azimuth_one = np.linspace(0.0, 360.0, nrays_per_sweep, endpoint=False, dtype=np.float64)
        azimuth = np.concatenate([azimuth_one, azimuth_one])
        elevation = np.concatenate(
            [np.full(nrays_per_sweep, sweeps[0], dtype=np.float64), np.full(nrays_per_sweep, sweeps[1], dtype=np.float64)]
        )
        sin_az = np.sin(np.deg2rad(azimuth))[:, np.newaxis]
        cos_az = np.cos(np.deg2rad(azimuth))[:, np.newaxis]
        cos_el = np.cos(np.deg2rad(elevation))[:, np.newaxis]
        clean_velocity = (u * sin_az + v * cos_az) * cos_el + 0.5
        noisy_velocity = clean_velocity + 2.0 * np.sin(np.linspace(0.0, 6.0 * np.pi, azimuth.size))[:, np.newaxis]
        clean_velocity = np.repeat(clean_velocity, nbins, axis=1)
        noisy_velocity = np.repeat(noisy_velocity, nbins, axis=1)

        fields = {
            "V": noisy_velocity.astype(np.float32) if include_vc else clean_velocity.astype(np.float32),
            "dBZ": np.full((azimuth.size, nbins), 20.0, dtype=np.float32),
        }
        if include_vc:
            fields["Vc"] = clean_velocity.astype(np.float32)

        start = np.datetime64("2026-03-20T00:00:00")
        time = np.array([start + np.timedelta64(i, "s") for i in range(azimuth.size)])
        sweep_start = np.array([0, nrays_per_sweep], dtype=np.int32)
        sweep_end = np.array([nrays_per_sweep - 1, azimuth.size - 1], dtype=np.int32)
        bins_per_sweep = np.array([nbins, nbins], dtype=np.int32)
        nyquist = np.array([25.0, 25.0], dtype=np.float64)
        unambiguous_range = np.array([150000.0, 150000.0], dtype=np.float64)
        return PRD(
            fields=fields,
            scan_type="ppi",
            time=time,
            range=gate_range,
            azimuth=azimuth,
            elevation=elevation,
            latitude=30.0,
            longitude=120.0,
            altitude=100.0,
            sweep_start_ray_index=sweep_start,
            sweep_end_ray_index=sweep_end,
            fixed_angle=sweeps,
            bins_per_sweep=bins_per_sweep,
            nyquist_velocity=nyquist,
            frequency=5.6,
            unambiguous_range=unambiguous_range,
            nrays=azimuth.size,
            nsweeps=sweeps.size,
            sitename="SYNTH",
        )

    @staticmethod
    def _sparsify_velocity_field(prd):
        sparse = prd.fields[0]["V"].values.copy()
        sparse[:, ::3] = np.nan
        sparse[::4, :] = np.nan
        sparse[10:15, 5:10] = sparse[10:15, 5:10] + 15.0
        prd.fields[0]["V"].values[:] = sparse
        prd._invalidate_cached_views()
        return prd

    def test_fit_vad_ring_recovers_uniform_horizontal_wind(self):
        from pycwr.retrieve.WindField import fit_vad_ring

        azimuth = np.linspace(0.0, 360.0, 72, endpoint=False)
        elevation = np.full(azimuth.shape, 1.5)
        u_true = 8.0
        v_true = -4.0
        radial_velocity = (u_true * np.sin(np.deg2rad(azimuth)) + v_true * np.cos(np.deg2rad(azimuth))) * np.cos(
            np.deg2rad(elevation)
        ) + 0.5

        result = fit_vad_ring(azimuth, elevation, radial_velocity)

        self.assertAlmostEqual(result["u"], u_true, places=6)
        self.assertAlmostEqual(result["v"], v_true, places=6)
        self.assertAlmostEqual(result["radial_offset"], 0.5, places=6)
        self.assertGreaterEqual(result["azimuth_coverage_deg"], 350.0)

    def test_retrieve_vad_on_synthetic_prd_recovers_truth(self):
        prd = self._build_uniform_wind_prd()

        dataset = prd.retrieve_vad(sweeps=0, field_name="V")
        valid = np.isfinite(dataset["u"].values)

        self.assertGreater(valid.sum(), 0)
        self.assertTrue(np.allclose(dataset["u"].values[valid], 8.0, atol=1.0e-5))
        self.assertTrue(np.allclose(dataset["v"].values[valid], -4.0, atol=1.0e-5))
        self.assertIn("Browning and Wexler", json.loads(dataset.attrs["references"])[0])

    def test_retrieve_vvp_on_synthetic_prd_recovers_truth(self):
        prd = self._build_uniform_wind_prd()

        dataset = prd.retrieve_vvp(0, field_name="V", az_num=9, bin_num=5)
        valid = np.isfinite(dataset["u"].values)

        self.assertGreater(valid.sum(), 0)
        self.assertAlmostEqual(float(np.nanmean(dataset["u"].values)), 8.0, places=5)
        self.assertAlmostEqual(float(np.nanmean(dataset["v"].values)), -4.0, places=5)
        self.assertIn("Waldteufel and Corbin", json.loads(dataset.attrs["references"])[0])

    def test_default_velocity_field_prefers_vc(self):
        prd = self._build_uniform_wind_prd(include_vc=True)

        dataset = prd.retrieve_vad(sweeps=0)
        used = json.loads(dataset.attrs["velocity_field_used_by_sweep"])

        self.assertEqual(used["0"], "Vc")
        valid = np.isfinite(dataset["u"].values)
        self.assertTrue(np.allclose(dataset["u"].values[valid], 8.0, atol=1.0e-5))
        self.assertTrue(np.allclose(dataset["v"].values[valid], -4.0, atol=1.0e-5))

    def test_retrieve_vwp_remains_usable_with_sparse_velocity_data(self):
        prd = self._sparsify_velocity_field(self._build_uniform_wind_prd())

        profile = prd.retrieve_vwp(sweeps=0, height_step=300.0, max_gap_bins=1)

        valid = np.isfinite(profile["u"].values)
        self.assertGreater(valid.sum(), 0)
        self.assertAlmostEqual(float(np.nanmean(profile["u"].values)), 8.0, places=1)
        self.assertAlmostEqual(float(np.nanmean(profile["v"].values)), -4.0, places=1)

    def test_add_product_vwp_writes_profile_variables(self):
        prd = self._build_uniform_wind_prd()

        profile = prd.add_product_VWP(sweeps=0, height_step=300.0)

        self.assertEqual(profile.attrs["method"], "VWP")
        self.assertIn("VWP_u", prd.product)
        self.assertIn("VWP_speed", prd.product)
        self.assertIn("z_vwp", prd.product.coords)

    def test_backward_compatible_helpers_still_work(self):
        from pycwr.retrieve import VAD, VVP

        azimuth = np.linspace(0.0, 360.0, 72, endpoint=False)
        elevation = 1.0
        velocity = np.repeat(
            (
                (
                    5.0 * np.sin(np.deg2rad(azimuth)) + 3.0 * np.cos(np.deg2rad(azimuth))
                )
                * np.cos(np.deg2rad(elevation))
            )[:, np.newaxis],
            12,
            axis=1,
        )

        u, v = VAD(azimuth, elevation, velocity[:, 0])
        field_u, field_v = VVP(azimuth, elevation, velocity, az_num=9, bin_num=5)

        self.assertAlmostEqual(u, 5.0, places=6)
        self.assertAlmostEqual(v, 3.0, places=6)
        self.assertEqual(field_u.shape, velocity.shape)
        self.assertEqual(field_v.shape, velocity.shape)

    def test_easy_plot_wind_products_smoke(self):
        from pycwr.draw import plot_vvp, plot_wind_profile

        prd = self._build_uniform_wind_prd()
        profile_result = plot_wind_profile(prd)
        vvp_result = plot_vvp(prd, sweep=0, az_num=31, bin_num=5, background_field="dBZ")

        self.assertIsNotNone(profile_result.artist)
        self.assertIsNotNone(vvp_result.artist)

    def test_real_sample_wind_retrieval_smoke(self):
        from pycwr.io import read_auto

        sample = self._sample_z9046()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        vad = prd.retrieve_vad(sweeps=0, max_range_km=40.0, gate_step=4)
        vvp = prd.retrieve_vvp(0, max_range_km=20.0, az_num=91, bin_num=5, azimuth_step=12, range_step=4)
        profile = prd.retrieve_vwp(sweeps=[0, 1, 2], max_range_km=40.0, gate_step=4, height_step=500.0)

        self.assertEqual(vad.attrs["method"], "VAD")
        self.assertEqual(vvp.attrs["method"], "VVP")
        self.assertEqual(profile.attrs["method"], "VWP")
        self.assertIn("velocity_field_used_by_sweep", vad.attrs)
        self.assertIn("velocity_field_used", vvp.attrs)
        self.assertGreater(int(np.isfinite(vad["u"].values).sum()), 0)
        self.assertGreater(int(np.isfinite(vvp["u"].values).sum()), 0)
        self.assertGreater(int(np.isfinite(profile["u"].values).sum()), 0)


if __name__ == "__main__":
    unittest.main()
