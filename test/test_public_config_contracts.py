"""Configuration and lazy-mapping contracts independent of radar fixtures."""
import tempfile
import unittest
import warnings
from pathlib import Path
from unittest.mock import Mock

from pycwr.configure import pyart_config as config
from pycwr.configure.pyart_default_config import spectrum_width_limit, velocity_limit
from pycwr.configure.pyart_lazydict import LazyLoadDict


class ConfigurationContracts(unittest.TestCase):
    def tearDown(self):
        config.load_config()

    def test_failed_config_load_preserves_active_metadata(self):
        before = config.get_metadata("reflectivity")
        self.assertTrue(before)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.py"
            path.write_text("DEFAULT_METADATA = {}\n", encoding="utf-8")
            with self.assertRaises(AttributeError):
                config.load_config(str(path))
        self.assertEqual(config.get_metadata("reflectivity"), before)

    def test_metadata_and_mapping_are_independent_copies(self):
        data = config.get_metadata("reflectivity")
        data["units"] = "changed"
        self.assertNotEqual(config.get_metadata("reflectivity")["units"], "changed")
        mapping = config.get_field_mapping("nexrad_archive")
        self.assertTrue(mapping)
        mapping.clear()
        self.assertTrue(config.get_field_mapping("nexrad_archive"))
        self.assertEqual(config.get_metadata("not-a-field"), {})
        self.assertEqual(config.get_field_limits("not-a-field"), (None, None))
        self.assertIsInstance(config.get_fillvalue(), float)
        self.assertEqual(config.get_field_name("reflectivity"), "reflectivity")

    def test_fallback_colormap_respects_rcparams_without_warning(self):
        import matplotlib
        with matplotlib.rc_context({"image.cmap": "plasma"}), warnings.catch_warnings():
            warnings.simplefilter("error")
            self.assertEqual(config.get_field_colormap("not-a-field"), "plasma")
        self.assertTrue(config.get_field_colormap("reflectivity"))

    def test_file_metadata_priority_and_field_selection(self):
        metadata = config.FileMetadata(
            "nexrad_archive", field_names={"REF": "reflectivity", "VEL": "velocity"},
            additional_metadata={"reflectivity": {"units": "custom"}},
            include_fields=["reflectivity", "velocity"], exclude_fields=["velocity"],
        )
        self.assertEqual(metadata("reflectivity"), {"units": "custom"})
        self.assertEqual(metadata.get_field_name("REF"), "reflectivity")
        self.assertIsNone(metadata.get_field_name("VEL"))
        self.assertIsNone(metadata.get_field_name("unknown"))
        self.assertEqual(metadata("unknown"), {})
        metadata.get_metadata("reflectivity")["units"] = "changed"
        self.assertEqual(metadata("reflectivity")["units"], "custom")
        self.assertTrue(metadata("time"))
        raw = config.FileMetadata("unknown", file_field_names=True)
        self.assertEqual(raw.get_field_name("raw"), "raw")
        self.assertTrue(raw("reflectivity"))
        self.assertIsNone(config.FileMetadata("unknown", include_fields=[]).get_field_name("raw"))

    def test_limits_select_sweep_and_default_without_container(self):
        container = Mock(nsweeps=2)
        container.get_nyquist_vel.return_value = 15.0
        self.assertEqual(velocity_limit(container, 1), (-15.0, 15.0))
        container.get_nyquist_vel.assert_called_with(1, check_uniform=False)
        self.assertEqual(spectrum_width_limit(container, 10), (0.0, 15.0))
        container.get_nyquist_vel.assert_called_with(0, check_uniform=False)
        self.assertEqual(velocity_limit(), (-30.0, 30.0))
        self.assertEqual(spectrum_width_limit(), (0, 30.0))
        self.assertEqual(config.get_field_limits("velocity", container), (-15.0, 15.0))


class LazyMappingContracts(unittest.TestCase):
    def test_lazy_copy_cache_and_mutation(self):
        loader = Mock(return_value=42)
        mapping = LazyLoadDict({"existing": 1})
        mapping.set_lazy("lazy", loader)
        copy = mapping.copy()
        self.assertEqual(set(mapping), {"existing", "lazy"})
        self.assertEqual(len(mapping), 2)
        loader.assert_not_called()
        self.assertIn("LazyLoad", str(mapping))
        self.assertEqual(copy["lazy"], 42)
        self.assertEqual(copy["lazy"], 42)
        self.assertEqual(loader.call_count, 1)
        self.assertTrue(mapping.has_key("lazy"))
        self.assertEqual(loader.call_count, 2)
        copy["existing"] = 9
        self.assertEqual(mapping["existing"], 1)
        mapping.set_lazy("existing", lambda: 2)
        del mapping["existing"]
        self.assertNotIn("existing", mapping)
        mapping.set_lazy("override", lambda: 5)
        mapping["override"] = 8
        self.assertEqual(mapping["override"], 8)
        del mapping["override"]
        with self.assertRaises(KeyError):
            mapping["missing"]
