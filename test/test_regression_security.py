import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock


class ReaderSecurityTests(unittest.TestCase):
    def test_read_all_enforces_maximum_size(self):
        from pycwr.io.util import _read_all

        with self.assertRaises(ValueError):
            _read_all(io.BytesIO(b"0123456789"), "test payload", max_bytes=8)

    def test_wsr98d_reader_uses_larger_payload_cap_for_real_world_archives(self):
        from pycwr.io.WSR98DFile import WSR98DBaseData, WSR98D_MAX_DECODED_PAYLOAD_BYTES

        reader = WSR98DBaseData.__new__(WSR98DBaseData)
        reader.fid = io.BytesIO(b"")

        with mock.patch("pycwr.io.WSR98DFile._read_all", return_value=b"") as read_all:
            radial, status, azimuth, elevation, seconds, microseconds = reader._parse_radial()

        self.assertEqual(radial, [])
        self.assertEqual(status.size, 0)
        self.assertEqual(read_all.call_args.kwargs["max_bytes"], WSR98D_MAX_DECODED_PAYLOAD_BYTES)

    def test_pa_signature_uses_explicit_validation(self):
        from pycwr.io.PAFile import PABaseData

        reader = PABaseData.__new__(PABaseData)
        reader.fid = io.BytesIO(b"\x00" * 12)

        with self.assertRaises(ValueError):
            reader._check_standard_basedata()

    def test_sab_rejects_out_of_bounds_field_offsets(self):
        from pycwr.io.SABFile import SABBaseData

        reader = SABBaseData.__new__(SABBaseData)
        radial_buf = b"\x00" * 2432
        header = {
            "GatesNumberOfReflectivity": 10,
            "PtrOfReflectivity": 4096,
            "GatesNumberOfDoppler": 10,
            "PtrOfVelocity": 0,
            "PtrOfSpectrumWidth": 0,
        }

        with mock.patch("pycwr.io.SABFile._unpack_from_buf", return_value=(header, 128)):
            with self.assertRaises(ValueError):
                reader._parse_radial_single(radial_buf)

    def test_radar_format_returns_none_for_unrelated_archive(self):
        from pycwr.io.util import radar_format

        with tempfile.TemporaryDirectory() as tmpdir:
            archive = Path(tmpdir) / "package.tar.gz"
            archive.write_bytes(b"not-a-radar-archive")
            self.assertIsNone(radar_format(str(archive)))

    def test_radar_format_rejects_false_positive_sab_signature(self):
        from pycwr.io.util import radar_format

        with tempfile.TemporaryDirectory() as tmpdir:
            bogus = Path(tmpdir) / "module.pyc"
            payload = bytearray(b"\x00" * 64)
            payload[14:16] = b"\x01\x00"
            bogus.write_bytes(bytes(payload))
            self.assertIsNone(radar_format(str(bogus)))


class WebViewerSecurityTests(unittest.TestCase):
    def _create_app(self, allowed_root):
        try:
            from pycwr.GraphicalInterface.web_app import create_app
        except ImportError as exc:  # pragma: no cover
            self.skipTest(f"Flask viewer dependencies are unavailable: {exc}")
        return create_app(allowed_roots=[allowed_root], auth_token="test-token")

    def test_api_requires_token(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            app = self._create_app(tmpdir)
            client = app.test_client()
            response = client.get("/api/files", query_string={"dir": tmpdir})
            self.assertEqual(response.status_code, 403)

    def test_api_files_rejects_directory_outside_allowed_roots(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            app = self._create_app(tmpdir)
            client = app.test_client()
            outside_dir = Path(tmpdir).parent
            response = client.get(
                "/api/files",
                query_string={"dir": str(outside_dir), "token": "test-token"},
            )
            self.assertEqual(response.status_code, 403)

    def test_api_metadata_rejects_file_outside_allowed_roots(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            app = self._create_app(tmpdir)
            client = app.test_client()
            outside_file = Path(tmpdir).parent / "outside.bin"
            outside_file.write_bytes(b"not-a-radar-file")
            try:
                response = client.get(
                    "/api/metadata",
                    query_string={"path": str(outside_file), "token": "test-token"},
                )
                self.assertEqual(response.status_code, 403)
            finally:
                outside_file.unlink(missing_ok=True)

    def test_api_files_allows_scanning_within_allowed_roots(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            allowed = Path(tmpdir)
            (allowed / "sample.bin").write_bytes(b"RSTM" + b"\x00" * 128)
            app = self._create_app(tmpdir)
            client = app.test_client()
            response = client.get(
                "/api/files",
                query_string={"dir": str(allowed), "token": "test-token"},
            )
            self.assertEqual(response.status_code, 200)
            payload = response.get_json()
            self.assertTrue(payload["ok"])
            self.assertEqual(len(payload["files"]), 1)
            self.assertEqual(payload["files"][0]["name"], "sample.bin")

    def test_api_tree_recursively_includes_nested_supported_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            allowed = Path(tmpdir)
            nested = allowed / "nested" / "day1"
            nested.mkdir(parents=True)
            (nested / "sample.bin").write_bytes(b"RSTM" + b"\x00" * 128)
            app = self._create_app(tmpdir)
            client = app.test_client()
            response = client.get(
                "/api/tree",
                query_string={"dir": str(allowed), "token": "test-token"},
            )

            self.assertEqual(response.status_code, 200)
            payload = response.get_json()
            self.assertTrue(payload["ok"])
            self.assertEqual(payload["tree"]["file_count"], 1)
            nested_node = payload["tree"]["children"][0]["children"][0]["children"][0]
            self.assertEqual(nested_node["type"], "file")
            self.assertEqual(nested_node["name"], "sample.bin")

    def test_api_catalog_groups_files_by_station_and_time(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            allowed = Path(tmpdir)
            (allowed / "Z_RADR_I_Z9046_20260317065928_O_DOR_SAD_CAP_FMT.bin.bz2").write_bytes(b"RSTM" + b"\x00" * 128)
            (allowed / "Z_RADR_I_Z9046_20260317070448_O_DOR_SAD_CAP_FMT.bin.bz2").write_bytes(b"RSTM" + b"\x00" * 128)
            (allowed / "Z_RADR_I_Z9250_20260317065944_O_DOR_SAD_CAP_FMT.bin.bz2").write_bytes(b"RSTM" + b"\x00" * 128)
            app = self._create_app(tmpdir)
            client = app.test_client()
            response = client.get(
                "/api/catalog",
                query_string={"dir": str(allowed), "token": "test-token"},
            )

            self.assertEqual(response.status_code, 200)
            payload = response.get_json()
            self.assertTrue(payload["ok"])
            stations = payload["catalog"]["stations"]
            self.assertEqual([item["station_id"] for item in stations], ["Z9046", "Z9250"])
            self.assertEqual(stations[0]["file_count"], 2)
            self.assertEqual(stations[0]["files"][0]["scan_time"], "2026-03-17T06:59:28")

    def test_default_app_allows_local_directory_scan(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                from pycwr.GraphicalInterface.web_app import create_app
            except ImportError as exc:  # pragma: no cover
                self.skipTest(f"Flask viewer dependencies are unavailable: {exc}")

            sample = Path(tmpdir) / "sample.bin"
            sample.write_bytes(b"RSTM" + b"\x00" * 128)
            app = create_app(auth_token="test-token")
            client = app.test_client()
            response = client.get(
                "/api/files",
                query_string={"dir": tmpdir, "token": "test-token"},
            )

            self.assertEqual(response.status_code, 200)
            payload = response.get_json()
            self.assertTrue(payload["ok"])
            self.assertEqual(len(payload["files"]), 1)
            self.assertEqual(payload["files"][0]["name"], "sample.bin")

    def test_api_files_requires_supported_filename_suffix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            allowed = Path(tmpdir)
            sample = allowed / "sample.dat"
            sample.write_bytes(b"RSTM" + b"\x00" * 128)
            app = self._create_app(tmpdir)
            client = app.test_client()
            response = client.get(
                "/api/files",
                query_string={"dir": str(allowed), "token": "test-token"},
            )

            self.assertEqual(response.status_code, 200)
            payload = response.get_json()
            self.assertTrue(payload["ok"])
            self.assertEqual(len(payload["files"]), 0)

    def test_api_tree_ignores_non_radar_archives(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            allowed = Path(tmpdir)
            (allowed / "package.tar.gz").write_bytes(b"not-a-radar-archive")
            app = self._create_app(tmpdir)
            client = app.test_client()
            response = client.get(
                "/api/tree",
                query_string={"dir": str(allowed), "token": "test-token"},
            )

            self.assertEqual(response.status_code, 200)
            payload = response.get_json()
            self.assertTrue(payload["ok"])
            self.assertEqual(payload["tree"]["file_count"], 0)
            self.assertEqual(payload["tree"]["children"], [])

    def test_api_tree_ignores_cache_and_compiled_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            allowed = Path(tmpdir)
            cache_dir = allowed / "__pycache__"
            cache_dir.mkdir()
            bogus = cache_dir / "NRadar.cpython-312.pyc"
            payload = bytearray(b"\x00" * 64)
            payload[14:16] = b"\x01\x00"
            bogus.write_bytes(bytes(payload))
            app = self._create_app(tmpdir)
            client = app.test_client()
            response = client.get(
                "/api/tree",
                query_string={"dir": str(allowed), "token": "test-token"},
            )

            self.assertEqual(response.status_code, 200)
            payload = response.get_json()
            self.assertTrue(payload["ok"])
            self.assertEqual(payload["tree"]["file_count"], 0)
            self.assertEqual(payload["tree"]["children"], [])

    def test_plot_section_returns_png(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            allowed = Path(tmpdir)
            sample = allowed / "sample.bin"
            sample.write_bytes(b"placeholder")
            app = self._create_app(tmpdir)
            client = app.test_client()

            from test_examples_sections import SectionExtractionTests

            prd = SectionExtractionTests()._build_ppi_prd()
            with mock.patch("pycwr.GraphicalInterface.web_app.RadarFileCache.get", return_value=prd):
                response = client.get(
                    "/plot/section.png",
                    query_string={
                        "path": str(sample),
                        "token": "test-token",
                        "field": "dBZ",
                        "start_x_km": -10.0,
                        "start_y_km": 0.0,
                        "end_x_km": 10.0,
                        "end_y_km": 0.0,
                    },
                )

            self.assertEqual(response.status_code, 200)
            self.assertEqual(response.mimetype, "image/png")
            self.assertTrue(response.data.startswith(b"\x89PNG"))


if __name__ == "__main__":
    unittest.main()
