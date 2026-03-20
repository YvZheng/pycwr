import importlib
import tempfile
import unittest
from pathlib import Path

import numpy as np

from pycwr.configure.default_config import CINRAD_field_mapping
from pycwr.io.NEXRADLevel2File import NEXRAD_MSG1_FIELD_SPECS, NEXRAD_MSG31_FIELD_SPECS
from pycwr.io.WSR98DFile import WSR98D_WRITE_FIELD_SPECS


class RadialExportTests(unittest.TestCase):
    @staticmethod
    def _repo_root():
        return Path(__file__).resolve().parents[2]

    @classmethod
    def _sample_path(cls):
        return (
            cls._repo_root()
            / "Z9046"
            / "Z_RADR_I_Z9046_20260317065928_O_DOR_SAD_CAP_FMT.bin.bz2"
        )

    @staticmethod
    def _wsr98d_tolerance(field_name):
        spec = WSR98D_WRITE_FIELD_SPECS[field_name]
        return (1.0 / float(spec["scale"])) + 1e-6

    @staticmethod
    def _nexrad_tolerance(spec):
        return (1.0 / float(spec["scale"])) + 1e-6

    def _require_pyart(self):
        try:
            return importlib.import_module("pyart")
        except ImportError as exc:  # pragma: no cover
            self.skipTest(f"Py-ART is unavailable: {exc}")

    def test_wsr98d_roundtrip_with_project_reader(self):
        from pycwr.io import read_WSR98D, read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "export_wsr98d.bin"
            returned = prd.to_wsr98d(str(out), overwrite=True)
            self.assertEqual(Path(returned), out)

            roundtrip = read_WSR98D(str(out))
            self.assertEqual(int(roundtrip.nsweeps), int(prd.nsweeps))
            self.assertEqual(int(roundtrip.nrays), int(prd.nrays))
            np.testing.assert_allclose(
                np.asarray(roundtrip.scan_info["fixed_angle"].values, dtype=np.float64),
                np.asarray(prd.scan_info["fixed_angle"].values, dtype=np.float64),
                atol=1e-3,
                rtol=0.0,
            )
            for field_name in ["dBT", "dBZ", "V", "W", "ZDR", "CC", "PhiDP", "KDP", "SNRH"]:
                for sweep in range(int(prd.nsweeps)):
                    source = np.asarray(prd.get_sweep_field(sweep, field_name, range_mode=None).values, dtype=np.float64)
                    target = np.asarray(roundtrip.get_sweep_field(sweep, field_name, range_mode=None).values, dtype=np.float64)
                    self.assertEqual(target.shape, source.shape)
                    np.testing.assert_allclose(
                        target,
                        source,
                        atol=self._wsr98d_tolerance(field_name),
                        rtol=0.0,
                        equal_nan=True,
                    )

    def test_module_level_write_wsr98d_matches_prd_export(self):
        from pycwr.io import read_WSR98D, read_auto, write_wsr98d

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "module_export_wsr98d.bin"
            returned = write_wsr98d(prd, str(out), overwrite=True)

            self.assertEqual(Path(returned), out)
            roundtrip = read_WSR98D(str(out))
            self.assertEqual(int(roundtrip.nsweeps), int(prd.nsweeps))
            self.assertEqual(int(roundtrip.nrays), int(prd.nrays))
            np.testing.assert_allclose(
                np.asarray(roundtrip.get_sweep_field(0, "dBZ", range_mode=None).values, dtype=np.float64),
                np.asarray(prd.get_sweep_field(0, "dBZ", range_mode=None).values, dtype=np.float64),
                atol=self._wsr98d_tolerance("dBZ"),
                rtol=0.0,
                equal_nan=True,
            )

    def test_nexrad_msg31_can_be_read_by_pyart_and_matches_prd(self):
        pyart = self._require_pyart()
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "export_msg31.ar2v"
            prd.to_nexrad_level2_msg31(str(out), overwrite=True)
            radar = pyart.io.read_nexrad_archive(str(out))

            self.assertEqual(radar.nsweeps, int(prd.nsweeps))
            self.assertEqual(radar.nrays, int(prd.nrays))
            np.testing.assert_allclose(
                np.asarray(radar.fixed_angle["data"], dtype=np.float64),
                np.asarray(prd.scan_info["fixed_angle"].values, dtype=np.float64),
                atol=(360.0 / 65536.0) + 1e-6,
                rtol=0.0,
            )
            for source_name, spec in NEXRAD_MSG31_FIELD_SPECS.items():
                pyart_name = CINRAD_field_mapping[source_name]
                self.assertIn(pyart_name, radar.fields)
                for sweep in range(int(prd.nsweeps)):
                    source = np.asarray(prd.get_sweep_field(sweep, source_name, range_mode=None).values, dtype=np.float64)
                    target = np.ma.filled(radar.get_field(sweep, pyart_name), np.nan).astype(np.float64)
                    self.assertGreaterEqual(target.shape[1], source.shape[1])
                    np.testing.assert_allclose(
                        target[:, : source.shape[1]],
                        source,
                        atol=self._nexrad_tolerance(spec),
                        rtol=0.0,
                        equal_nan=True,
                    )

    def test_nexrad_msg1_can_be_read_by_pyart_and_matches_core_fields(self):
        pyart = self._require_pyart()
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "export_msg1.ar2v"
            prd.to_nexrad_level2_msg1(str(out), overwrite=True)
            radar = pyart.io.read_nexrad_archive(str(out))

            self.assertEqual(radar.nsweeps, int(prd.nsweeps))
            self.assertEqual(radar.nrays, int(prd.nrays))
            self.assertIn("reflectivity", radar.fields)
            self.assertIn("velocity", radar.fields)
            self.assertIn("spectrum_width", radar.fields)
            for source_name, spec in NEXRAD_MSG1_FIELD_SPECS.items():
                pyart_name = CINRAD_field_mapping[source_name]
                for sweep in range(int(prd.nsweeps)):
                    source = np.asarray(prd.get_sweep_field(sweep, source_name, range_mode=None).values, dtype=np.float64)
                    target = np.ma.filled(radar.get_field(sweep, pyart_name), np.nan).astype(np.float64)
                    compare_width = min(target.shape[1], source.shape[1])
                    np.testing.assert_allclose(
                        target[:, :compare_width],
                        source[:, :compare_width],
                        atol=self._nexrad_tolerance(spec),
                        rtol=0.0,
                        equal_nan=True,
                    )

    def test_module_level_nexrad_writers_match_public_examples(self):
        pyart = self._require_pyart()
        from pycwr.io import read_auto, write_nexrad_level2_msg1, write_nexrad_level2_msg31

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        with tempfile.TemporaryDirectory() as tmpdir:
            msg31 = Path(tmpdir) / "module_export_msg31.ar2v"
            msg1 = Path(tmpdir) / "module_export_msg1.ar2v"

            returned_msg31 = write_nexrad_level2_msg31(prd, str(msg31), overwrite=True)
            returned_msg1 = write_nexrad_level2_msg1(prd, str(msg1), overwrite=True)

            self.assertEqual(Path(returned_msg31), msg31)
            self.assertEqual(Path(returned_msg1), msg1)

            radar31 = pyart.io.read_nexrad_archive(str(msg31))
            radar1 = pyart.io.read_nexrad_archive(str(msg1))

            self.assertEqual(radar31.nsweeps, int(prd.nsweeps))
            self.assertEqual(radar1.nsweeps, int(prd.nsweeps))
            self.assertIn("reflectivity", radar31.fields)
            self.assertIn("reflectivity", radar1.fields)

    def test_msg1_rejects_unsupported_requested_fields(self):
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "bad_msg1.ar2v"
            with self.assertRaises(ValueError):
                prd.to_nexrad_level2_msg1(
                    str(out),
                    field_names=["dBZ", "V", "W", "ZDR"],
                    overwrite=True,
                )


if __name__ == "__main__":
    unittest.main()
