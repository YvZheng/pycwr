import unittest
from pathlib import Path

import numpy as np


class PAReaderRegressionTests(unittest.TestCase):
    @staticmethod
    def _sample_path():
        return (
            Path(__file__).resolve().parents[2]
            / "Z_RADR_I_ZA460_20240808142000_O_DOR-XPD-CAP-FMT.BIN"
        )

    @staticmethod
    def _sample_zip_path():
        return (
            Path(__file__).resolve().parents[2]
            / "ZA460"
            / "Z_RADR_I_ZA460_20240808142000_O_DOR-XPD-CAP-FMT.BIN.zip"
        )

    def _require_sample(self):
        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("PA sample file is not available in this workspace")
        return sample

    def _require_zip_sample(self):
        sample = self._sample_zip_path()
        if not sample.exists():
            self.skipTest("PA zip sample file is not available in this workspace")
        return sample

    def test_radar_format_identifies_pa_sample(self):
        from pycwr.io.util import radar_format

        sample = self._require_sample()
        self.assertEqual(radar_format(str(sample)), "PA")

    def test_read_auto_and_read_pa_match_on_real_sample(self):
        from pycwr.io import read_PA, read_auto

        sample = self._require_sample()
        auto = read_auto(str(sample))
        direct = read_PA(str(sample))

        for prd in (auto, direct):
            self.assertEqual(str(prd.scan_info.scan_type.values), "ppi")
            self.assertEqual(int(prd.nsweeps), 39)
            self.assertEqual(int(prd.nrays), 10530)
            self.assertEqual(prd.fields[0]["dBZ"].shape, (270, 1472))
            np.testing.assert_array_equal(prd.scan_info["rays_per_sweep"].values, np.full((39,), 270))
            self.assertEqual(
                sorted(prd.available_fields()),
                ["CC", "KDP", "PhiDP", "SNRH", "V", "W", "ZDR", "dBT", "dBZ"],
            )
            summaries = prd.sweep_summary()
            self.assertEqual(len(summaries), 39)
            self.assertTrue(all(item["rays"] == 270 for item in summaries))
            self.assertAlmostEqual(float(summaries[0]["fixed_angle"]), 0.75, places=2)
            self.assertAlmostEqual(float(summaries[1]["fixed_angle"]), 2.25, places=2)

        np.testing.assert_array_equal(auto.scan_info.fixed_angle.values, direct.scan_info.fixed_angle.values)
        np.testing.assert_array_equal(auto.scan_info["rays_per_sweep"].values, direct.scan_info["rays_per_sweep"].values)
        np.testing.assert_allclose(auto.fields[0]["dBZ"].values, direct.fields[0]["dBZ"].values, equal_nan=True)

    def test_read_auto_accepts_single_file_pa_zip_sample(self):
        from pycwr.io import read_PA, read_auto

        sample = self._require_sample()
        zip_sample = self._require_zip_sample()
        zipped = read_auto(str(zip_sample))
        direct = read_PA(str(sample))

        self.assertEqual(int(zipped.nsweeps), int(direct.nsweeps))
        self.assertEqual(int(zipped.nrays), int(direct.nrays))
        np.testing.assert_array_equal(zipped.scan_info.fixed_angle.values, direct.scan_info.fixed_angle.values)
        np.testing.assert_array_equal(zipped.scan_info["rays_per_sweep"].values, direct.scan_info["rays_per_sweep"].values)
        np.testing.assert_allclose(zipped.fields[0]["dBT"].values, direct.fields[0]["dBT"].values, equal_nan=True)
        np.testing.assert_allclose(zipped.fields[1]["dBZ"].values, direct.fields[1]["dBZ"].values, equal_nan=True)

    def test_pa_sample_internal_pyart_export_smoke(self):
        from pycwr.core.PyartRadar import Radar as InternalRadar
        from pycwr.io import read_PA

        sample = self._require_sample()
        prd = read_PA(str(sample))
        radar = prd.to_pyart_radar(use_external=False, force_rebuild=True)

        self.assertIsInstance(radar, InternalRadar)
        self.assertEqual(radar.nrays, 10530)
        self.assertEqual(radar.ngates, 1472)
        self.assertIn("reflectivity", radar.fields)
        self.assertIn("velocity", radar.fields)

    def test_za460_sparse_pa_moments_are_preserved_through_prd(self):
        from pycwr.io import read_PA
        from pycwr.io.PAFile import PABaseData

        sample = self._require_sample()
        base = PABaseData(str(sample))
        prd = read_PA(str(sample))

        def finite_count(dataset, field_name):
            return int(np.isfinite(np.asarray(dataset[field_name].values, dtype=np.float64)).sum())

        # This sample stores total power widely, but most derived moments are
        # thresholded out by the source file itself. Ensure pycwr does not drop
        # additional data after raw decode.
        self.assertGreater(finite_count(prd.fields[0], "dBT"), 10_000)
        for field_name in ["dBZ", "V", "W", "ZDR", "CC", "PhiDP", "KDP", "SNRH"]:
            self.assertEqual(finite_count(prd.fields[0], field_name), 0)

        self.assertGreater(finite_count(prd.fields[1], "dBZ"), 0)
        self.assertGreater(finite_count(prd.fields[38], "dBZ"), 0)

        base_counts = []
        for sweep in (0, 1, 38):
            start = int(base.sweep_start_ray_index[sweep])
            end = int(base.sweep_end_ray_index[sweep])
            sweep_counts = {}
            for field_name in ["dBT", "dBZ", "V", "W", "ZDR", "CC", "PhiDP", "KDP", "SNRH"]:
                total = 0
                for iray in range(start, end + 1):
                    total += int(np.isfinite(np.asarray(base.radial[iray]["fields"][field_name], dtype=np.float64)).sum())
                sweep_counts[field_name] = total
            base_counts.append(sweep_counts)

        for sweep_index, expected_counts in zip((0, 1, 38), base_counts):
            for field_name, expected in expected_counts.items():
                self.assertEqual(finite_count(prd.fields[sweep_index], field_name), expected)
