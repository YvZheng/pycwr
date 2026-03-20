import tempfile
import unittest
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")


class HydrometeorExampleTests(unittest.TestCase):
    @staticmethod
    def _sample_path():
        return (
            Path(__file__).resolve().parents[2]
            / "Z9046"
            / "Z_RADR_I_Z9046_20260317065928_O_DOR_SAD_CAP_FMT.bin.bz2"
        )

    def test_temperature_profile_interpolation_example(self):
        from pycwr.retrieve import interpolate_temperature_profile

        gate_altitude = np.array([[0.0, 500.0], [1000.0, 1500.0]], dtype=np.float64)
        profile_height = np.array([0.0, 1000.0, 2000.0], dtype=np.float64)
        profile_temperature = np.array([20.0, 10.0, 0.0], dtype=np.float64)

        interpolated = interpolate_temperature_profile(gate_altitude, profile_height, profile_temperature)
        np.testing.assert_allclose(
            interpolated,
            np.array([[20.0, 15.0], [10.0, 5.0]], dtype=np.float64),
            rtol=0.0,
            atol=1e-6,
        )

    def test_array_classifier_supports_profile_free_mode(self):
        from pycwr.retrieve import classify_hydrometeors

        dbz = np.array([[20.0, 35.0], [45.0, 55.0]], dtype=np.float64)
        zdr = np.array([[0.2, 0.6], [1.0, 2.5]], dtype=np.float64)
        kdp = np.array([[0.0, 0.2], [0.8, 2.0]], dtype=np.float64)
        cc = np.array([[0.99, 0.98], [0.95, 0.9]], dtype=np.float64)

        hcl, confidence = classify_hydrometeors(
            dBZ=dbz,
            ZDR=zdr,
            KDP=kdp,
            CC=cc,
            T=None,
            band="C",
            return_confidence=True,
        )
        self.assertEqual(hcl.shape, dbz.shape)
        self.assertEqual(confidence.shape, dbz.shape)
        self.assertTrue(np.isfinite(hcl).all())
        self.assertTrue(np.isfinite(confidence).all())
        self.assertTrue(np.all((hcl >= 1.0) & (hcl <= 10.0)))
        self.assertTrue(np.all((confidence >= 0.0) & (confidence <= 1.0)))

    def test_prd_classification_with_profile_and_plot(self):
        from pycwr.GraphicalInterface.web_colors import build_web_style
        from pycwr.draw import plot_ppi
        from pycwr.draw._plot_core import resolve_field_style
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        self.assertNotIn("HCL", prd.available_fields(sweep=0))

        classified = prd.classify_hydrometeors(
            inplace=False,
            band="C",
            profile_height=np.array([0.0, 2000.0, 4000.0, 8000.0, 12000.0], dtype=np.float64),
            profile_temperature=np.array([24.0, 12.0, 2.0, -16.0, -40.0], dtype=np.float64),
            confidence_field="HCL_CONF",
            temperature_field="HCL_T",
        )

        self.assertNotIn("HCL", prd.available_fields(sweep=0))
        self.assertIn("HCL", classified.available_fields(sweep=0))
        self.assertIn("HCL_CONF", classified.available_fields(sweep=0))
        self.assertIn("HCL_T", classified.available_fields(sweep=0))

        hcl = np.asarray(classified.get_sweep_field(0, "HCL").values, dtype=np.float64)
        conf = np.asarray(classified.get_sweep_field(0, "HCL_CONF").values, dtype=np.float64)
        temp = np.asarray(classified.get_sweep_field(0, "HCL_T").values, dtype=np.float64)
        reference = np.asarray(classified.get_sweep_field(0, "ZDR").values, dtype=np.float64)

        self.assertEqual(hcl.shape, reference.shape)
        self.assertEqual(conf.shape, reference.shape)
        self.assertEqual(temp.shape, reference.shape)
        self.assertGreater(np.isfinite(hcl).sum(), 0)
        self.assertTrue(np.nanmin(hcl) >= 1.0)
        self.assertTrue(np.nanmax(hcl) <= 10.0)
        self.assertTrue(np.nanmin(conf) >= 0.0)
        self.assertTrue(np.nanmax(conf) <= 1.0)

        hcl_field = classified.get_sweep_field(0, "HCL")
        style = resolve_field_style(classified, 0, "HCL", hcl_field)
        web_style = build_web_style(classified, 0, "HCL", hcl_field)
        self.assertEqual(style.colorbar.ticklabels[0], "1 毛毛雨")
        self.assertEqual(style.colorbar.ticklabels[-1], "10 大水滴")
        self.assertEqual(web_style.colorbar.ticklabels[0], "1 毛毛雨")
        self.assertEqual(web_style.colorbar.ticklabels[-1], "10 大水滴")

        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "hcl.png"
            result = plot_ppi(classified, field="HCL", sweep=0, save=out, cbar=False)
            self.assertIsNotNone(result.artist)
            self.assertTrue(out.exists())

    def test_prd_inplace_classification_without_profile(self):
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample))
        returned = prd.add_hydrometeor_classification(
            band="C",
            confidence_field="HCL_CONF",
        )
        self.assertIs(returned, prd)
        self.assertIn("HCL", prd.available_fields(sweep=0))
        self.assertIn("HCL_CONF", prd.available_fields(sweep=0))
        attrs = prd.fields[0]["HCL"].attrs
        self.assertEqual(attrs["uses_temperature"], "false")
        self.assertIn("hydrometeor_classes", attrs)

    def test_web_viewer_metadata_exposes_hcl_legend(self):
        try:
            from pycwr.GraphicalInterface import create_app
        except ImportError as exc:  # pragma: no cover
            self.skipTest(f"Flask viewer dependencies are unavailable: {exc}")
        from pycwr.io import read_auto

        sample = self._sample_path()
        if not sample.exists():
            self.skipTest("sample radar file is not available in this workspace")

        prd = read_auto(str(sample)).classify_hydrometeors(inplace=False, band="C")
        with tempfile.TemporaryDirectory() as tmpdir:
            exported = Path(tmpdir) / "hid_roundtrip_wsr98d.bin"
            prd.to_wsr98d(str(exported), overwrite=True)

            app = create_app(allowed_roots=[tmpdir], auth_token="sample-token")
            client = app.test_client()
            response = client.get(
                "/api/metadata",
                query_string={"path": str(exported), "token": "sample-token"},
            )

            self.assertEqual(response.status_code, 200)
            payload = response.get_json()
            self.assertTrue(payload["ok"])
            self.assertIn("HCL", payload["fields"])
            self.assertEqual(payload["field_legends"]["HCL"][0]["label_zh"], "毛毛雨")
            self.assertEqual(payload["field_legends"]["HCL"][-1]["label_zh"], "大水滴")


if __name__ == "__main__":
    unittest.main()
