"""Boundary regressions for binary IO, timestamps and interop."""
import datetime
import importlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np


class ReaderBoundaryTests(unittest.TestCase):
    def test_format_detection_preserves_borrowed_stream(self):
        from pycwr.io.util import radar_format

        stream = io.BytesIO(b'RSTM' + b'\0' * 128)
        stream.seek(17)
        self.assertEqual(radar_format(stream), 'WSR98D')
        self.assertEqual(stream.tell(), 17)
        self.assertFalse(stream.closed)

    def test_pa_type_marker_without_magic_is_not_a_radar(self):
        from pycwr.io.util import radar_format

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'unrelated.bin'
            path.write_bytes(b'notradar' + b'\x10\0\0\0' + b'\0' * 128)
            self.assertIsNone(radar_format(path))

    def test_readers_close_file_after_invalid_header(self):
        for module_name, class_name in (
            ('WSR98DFile', 'WSR98DBaseData'), ('PAFile', 'PABaseData'),
            ('CCFile', 'CCBaseData'), ('SCFile', 'SCBaseData'),
        ):
            with self.subTest(reader=class_name):
                module = importlib.import_module('pycwr.io.' + module_name)
                stream = io.BytesIO(b'invalid')
                with mock.patch.object(module, '_prepare_for_read', return_value=stream):
                    with self.assertRaises(ValueError):
                        getattr(module, class_name)('invalid.bin')
                self.assertTrue(stream.closed)

    def test_sab_closes_file_when_read_fails(self):
        from pycwr.io.SABFile import SABBaseData

        stream = io.BytesIO(b'payload')
        with mock.patch('pycwr.io.SABFile._prepare_for_read', return_value=stream), \
                mock.patch('pycwr.io.SABFile._read_all', side_effect=ValueError('limit')):
            with self.assertRaises(ValueError):
                SABBaseData('sample.bin')
        self.assertTrue(stream.closed)


class TimeConversionBoundaryTests(unittest.TestCase):
    units = 'seconds since 2026-01-01T00:00:00Z'

    def test_python_datetime_scalar(self):
        from pycwr.io.util import date2num

        self.assertEqual(date2num(datetime.datetime(2026, 1, 1, 0, 0, 2), self.units), 2.0)

    def test_numpy_datetime_scalar_and_array(self):
        from pycwr.io.util import date2num

        self.assertEqual(date2num(np.datetime64('2026-01-01T00:00:02.5'), self.units), 2.5)
        dates = np.array([['2026-01-01T00:00:01', '2026-01-01T00:00:02']], dtype='datetime64[s]')
        np.testing.assert_array_equal(date2num(dates, self.units), [[1.0, 2.0]])

    def test_invalid_units_are_explicitly_rejected(self):
        from pycwr.io.util import date2num

        with self.assertRaises(ValueError):
            date2num([], 'hours since 2026-01-01T00:00:00Z')


class ExportBoundaryTests(unittest.TestCase):
    def test_wsr98d_preserves_fractional_gate_spacing(self):
        from test_examples_sections import SectionExtractionTests
        from pycwr.io import read_auto

        ranges = np.arange(1, 9) * 62.5
        prd = SectionExtractionTests()._build_ppi_prd(ranges=ranges)
        for sweep in range(prd.nsweeps):
            prd.fields[sweep] = prd.fields[sweep].assign_coords(
                time=np.datetime64('2026-01-01') + np.arange(8) * np.timedelta64(1, 's'))
        prd.scan_info['start_time'] = np.datetime64('2026-01-01')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'fractional.bin'
            prd.to_wsr98d(path, field_names=['dBZ'])
            result = read_auto(path)
        np.testing.assert_allclose(result.fields[0].range.values, ranges, atol=1e-6)

    def test_xradar_datetime_sweeps_can_be_written_to_netcdf(self):
        import xarray as xr
        from test_examples_sections import SectionExtractionTests

        prd = SectionExtractionTests()._build_ppi_prd()
        times = np.datetime64("2026-01-01T00:00:00", "ns") + np.arange(8) * np.timedelta64(123456789, "ns")
        for sweep in range(prd.nsweeps):
            prd.fields[sweep] = prd.fields[sweep].assign_coords(time=times)
        prd.scan_info["start_time"] = times[0]
        dataset = prd.to_xradar_sweeps()["sweep_0"]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sweep.nc"
            dataset.to_netcdf(path)
            with xr.open_dataset(path) as restored:
                np.testing.assert_array_equal(restored.time.values, times)
                np.testing.assert_allclose(restored.DBZ.values, dataset.DBZ.values, equal_nan=True)

    def test_msg31_preserves_extreme_valid_values_with_uniform_volume_encoding(self):
        try:
            import pyart
        except ImportError:
            self.skipTest("Py-ART is unavailable")
        from test_examples_sections import SectionExtractionTests
        from pycwr.io import write_nexrad_level2_msg31

        prd = SectionExtractionTests()._build_ppi_prd()
        for sweep in range(prd.nsweeps):
            prd.fields[sweep] = prd.fields[sweep].assign_coords(
                time=np.datetime64("2026-01-01") + np.arange(8) * np.timedelta64(1, "s"))
        prd.scan_info["start_time"] = np.datetime64("2026-01-01")
        prd.fields[0]["ZDR"].values[0, 0] = -9.0
        prd.fields[1]["ZDR"].values[1, 1] = 12.0
        prd.fields[1]["ZDR"].values[2, 1] = np.nan
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.ar2v"
            write_nexrad_level2_msg31(prd, path)
            radar = pyart.io.read_nexrad_archive(str(path))
            raw = pyart.io.nexrad_level2.NEXRADLevel2File(str(path))
            try:
                encodings = {(record["ZDR"]["scale"], record["ZDR"]["offset"],
                              record["ZDR"]["word_size"]) for record in raw.radial_records}
                self.assertEqual(len(encodings), 1)
                self.assertEqual(next(iter(encodings))[0], 16.0)
            finally:
                raw.close()
            for sweep in range(prd.nsweeps):
                for name, target, tolerance in (("dBZ", "reflectivity", 0.5),
                                                 ("ZDR", "differential_reflectivity", 1 / 16)):
                    np.testing.assert_allclose(radar.get_field(sweep, target).filled(np.nan),
                                               prd.fields[sweep][name].values,
                                               equal_nan=True, atol=tolerance, rtol=0)

    def test_wsr98d_encoding_rejects_unrepresentable_values(self):
        from pycwr.io.WSR98DFile import WSR98D_WRITE_FIELD_SPECS, _fit_wsr98d_encoding

        with self.assertRaisesRegex(ValueError, "encoding range"):
            _fit_wsr98d_encoding([-1e20, 1e20], WSR98D_WRITE_FIELD_SPECS["dBZ"])

    def test_wsr98d_quantization_keeps_valid_values_out_of_missing_codes(self):
        from pycwr.io.WSR98DFile import (
            WSR98D_WRITE_FIELD_SPECS, _encode_quantized, _fit_wsr98d_encoding,
        )

        for field, values in (
            ("SQI", [0.0, 0.005, 1.0, np.nan]),
            ("SQI", [0.0, 0.005, 1.0, np.inf]),
            ("SQI", [0.0, 0.005, 1.0, -np.inf]),
            ("KDP", [-12.5, 0.0, 25.0, np.nan]),
            ("dBZ", [-40.0, 0.0, 100.0, np.nan]),
        ):
            with self.subTest(field=field):
                source = np.array(values)
                spec = _fit_wsr98d_encoding(source, WSR98D_WRITE_FIELD_SPECS[field])
                encoded = _encode_quantized(source, spec["scale"], spec["offset"], spec["bin_length"])
                codes = np.frombuffer(encoded, dtype='u1' if spec["bin_length"] == 1 else '<u2')
                self.assertTrue(np.all(codes[np.isfinite(source)] >= 5))
                self.assertEqual(codes[-1], 3)
                decoded = (codes[:-1].astype(float) - spec["offset"]) / spec["scale"]
                np.testing.assert_allclose(decoded, source[:-1], atol=1 / spec["scale"], rtol=0)

    def test_wsr98d_half_steps_use_same_quantization_when_fitting_and_encoding(self):
        from pycwr.io.WSR98DFile import _encode_quantized, _fit_wsr98d_encoding

        for upper in (25.05, 6553.05):
            values = np.array([-0.05, 0.05, 0.15, upper, np.nan])
            spec = _fit_wsr98d_encoding(values, {"scale": 10, "offset": 0, "bin_length": 1})
            encoded = _encode_quantized(values, spec["scale"], spec["offset"], spec["bin_length"])
            codes = np.frombuffer(encoded, dtype='u1' if spec["bin_length"] == 1 else '<u2')
            expected = np.rint(values[:-1] * spec["scale"]) + spec["offset"]
            np.testing.assert_array_equal(codes[:-1], expected)
            self.assertTrue(np.all(codes[:-1] >= 5))
            self.assertEqual(codes[-1], 3)

    def test_nexrad_rejects_fractional_gate_geometry(self):
        from pycwr.io.NEXRADLevel2File import _range_geometry

        for ranges in ([62.5, 125.0, 187.5], [0.5, 250.5, 500.5], [62.5]):
            with self.subTest(ranges=ranges):
                with self.assertRaisesRegex(ValueError, "integer.*met"):
                    _range_geometry(ranges)
        self.assertEqual(_range_geometry([250.0, 500.0, 750.0]), (250, 250))

    def test_writers_reject_nonuniform_range_grids(self):
        from pycwr.io.WSR98DFile import _range_geometry as wsr_geometry
        from pycwr.io.NEXRADLevel2File import _range_geometry as nexrad_geometry

        for geometry in (wsr_geometry, nexrad_geometry):
            with self.subTest(writer=geometry.__module__):
                with self.assertRaisesRegex(ValueError, 'uniform'):
                    geometry([250.0, 500.0, 900.0])


if __name__ == '__main__':
    unittest.main()
