"""Compatibility exports preserve data and ownership across Radar classes."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from test_public_science_contracts import make_radar


class CompatibilityExportContracts(unittest.TestCase):
    def test_legacy_pyart_alias_preserves_fields_and_cache(self):
        prd = make_radar()
        exported = prd.ToPyartRadar(use_external=False, field_names=["dBZ"])
        self.assertEqual(set(exported.fields), {"reflectivity"})
        np.testing.assert_allclose(exported.fields["reflectivity"]["data"], 30.0)
        self.assertEqual(exported.nsweeps, 3)
        self.assertIs(exported, prd.to_pyart_radar(use_external=False, field_names=["dBZ"]))

    def test_clone_new_radar_class_owns_data_and_metadata(self):
        from pycwr.core.PyartRadar import Radar
        from pycwr.core.interop import clone_radar_to_class

        class RadarSubclass(Radar):
            pass

        source = make_radar().to_pyart_radar(use_external=False)
        clone = clone_radar_to_class(source, RadarSubclass)
        self.assertIsInstance(clone, RadarSubclass)
        np.testing.assert_array_equal(clone.fields["reflectivity"]["data"], source.fields["reflectivity"]["data"])
        clone.fields["reflectivity"]["data"][0, 0] = -99
        clone.latitude["data"][0] = -45
        self.assertEqual(source.fields["reflectivity"]["data"][0, 0], 30)
        self.assertEqual(source.latitude["data"][0], 30)
        self.assertIs(clone_radar_to_class(source, Radar), source)

    def test_msg31_method_writes_archive_and_protects_existing_file(self):
        radar = make_radar()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "volume.ar2v"
            radar.to_nexrad_level2_msg31(path, field_names=["dBZ"])
            original = path.read_bytes()
            self.assertTrue(original.startswith(b"AR2V"))
            self.assertGreater(len(original), 24)
            with self.assertRaises(FileExistsError):
                radar.to_nexrad_level2_msg31(path, field_names=["dBZ"])
            self.assertEqual(path.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
