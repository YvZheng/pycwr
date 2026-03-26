import importlib
import tempfile
import unittest
from pathlib import Path

import numpy as np


class WSR98DSplitSweepRegressionTests(unittest.TestCase):
    @staticmethod
    def _sample_path():
        return (
            Path(__file__).resolve().parents[2]
            / "Z_RADR_I_ZA202_20250819095753_O_DOR_YLD2-D_CAP_FMT.bin.bz2"
        )

    def _require_pyart(self):
        try:
            return importlib.import_module("pyart")
        except ImportError as exc:  # pragma: no cover
            self.skipTest(f"Py-ART is unavailable: {exc}")

    def test_reader_merges_split_low_level_sweeps(self):
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("ZA202 split-sweep sample is not available in this workspace")

        prd = read_auto(str(sample))
        self.assertEqual(int(prd.nsweeps), 9)
        self.assertEqual(int(prd.nrays), 3290)
        self.assertEqual(prd.resolve_field_name("V", sweep=0), "Vc")
        self.assertEqual(prd.resolve_field_name("W", sweep=0), "W")
        self.assertIn("dBZ", prd.fields[0])
        self.assertIn("Zc", prd.fields[0])
        self.assertIn("Vc", prd.fields[0])
        self.assertIn("W", prd.fields[0])
        self.assertNotIn("V", prd.fields[0])
        self.assertEqual(prd.fields[0]["dBZ"].shape, prd.fields[0]["Vc"].shape)
        self.assertEqual(prd.fields[1]["dBZ"].shape, prd.fields[1]["Vc"].shape)
        self.assertIsNot(prd.fields[0]["dBZ"].variable, prd.fields[0]["Zc"].variable)

    def test_pyart_and_xradar_exports_keep_standard_views_without_merging_storage(self):
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("ZA202 split-sweep sample is not available in this workspace")

        prd = read_auto(str(sample))
        radar = prd.to_pyart_radar(use_external=False, force_rebuild=True)
        sweeps = prd.to_xradar_sweeps(force_rebuild=True)

        self.assertIn("velocity", radar.fields)
        self.assertIn("reflectivity", radar.fields)
        self.assertIn("corrected_reflectivity", radar.fields)
        self.assertIn("sweep_0", sweeps)
        self.assertIn("DBZ", sweeps["sweep_0"].data_vars)
        self.assertIn("VR", sweeps["sweep_0"].data_vars)
        self.assertIn("WRADH", sweeps["sweep_0"].data_vars)

    def test_wsr98d_roundtrip_preserves_corrected_variables(self):
        from pycwr.io import read_WSR98D, read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("ZA202 split-sweep sample is not available in this workspace")

        prd = read_auto(str(sample))
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "za202_roundtrip.bin"
            prd.to_wsr98d(str(out), overwrite=True)
            back = read_WSR98D(str(out))

        self.assertEqual(int(back.nsweeps), int(prd.nsweeps))
        self.assertEqual(int(back.nrays), int(prd.nrays))
        self.assertIn("dBZ", back.available_fields())
        self.assertIn("Zc", back.available_fields())
        self.assertIn("Vc", back.available_fields())
        self.assertNotIn("V", back.available_fields())
        np.testing.assert_allclose(
            np.asarray(back.get_sweep_field(0, "Vc", range_mode=None).values, dtype=np.float64),
            np.asarray(prd.get_sweep_field(0, "Vc", range_mode=None).values, dtype=np.float64),
            atol=0.5 + 1e-6,
            rtol=0.0,
            equal_nan=True,
        )

    def test_nexrad_exports_can_be_read_by_pyart(self):
        pyart = self._require_pyart()
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("ZA202 split-sweep sample is not available in this workspace")

        prd = read_auto(str(sample))
        with tempfile.TemporaryDirectory() as tmpdir:
            msg31 = Path(tmpdir) / "za202_msg31.ar2v"
            msg1 = Path(tmpdir) / "za202_msg1.ar2v"
            prd.to_nexrad_level2_msg31(str(msg31), overwrite=True)
            prd.to_nexrad_level2_msg1(str(msg1), overwrite=True)
            radar31 = pyart.io.read_nexrad_archive(str(msg31))
            radar1 = pyart.io.read_nexrad_archive(str(msg1))

        self.assertEqual(radar31.nsweeps, int(prd.nsweeps))
        self.assertEqual(radar1.nsweeps, int(prd.nsweeps))
        self.assertIn("reflectivity", radar31.fields)
        self.assertIn("velocity", radar31.fields)
        self.assertIn("reflectivity", radar1.fields)
        self.assertIn("velocity", radar1.fields)

    def test_qc_and_network_smoke_on_split_sweep_volume(self):
        from pycwr.core.NRadar import grid_3d_network_xy
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("ZA202 split-sweep sample is not available in this workspace")

        prd = read_auto(str(sample))
        processed = prd.apply_dualpol_qc(inplace=False)
        self.assertIn("Zc", processed.fields[0])
        self.assertIn("ZDRc", processed.fields[0])
        self.assertIn("QC_MASK", processed.fields[0])

        x = np.arange(-20_000.0, 20_001.0, 10_000.0)
        y = np.arange(-20_000.0, 20_001.0, 10_000.0)
        network = grid_3d_network_xy([prd], x, y, np.array([3000.0]), range_mode="native")
        self.assertIn("network_3d", network)
        self.assertEqual(network["network_3d"].shape, (1, x.size, y.size))


if __name__ == "__main__":
    unittest.main()
