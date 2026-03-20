import unittest
from pathlib import Path

import numpy as np


WSR98D_BASELINES = {
    "Z_RADR_I_Z9046_20260317065928_O_DOR_SAD_CAP_FMT.bin.bz2": {
        "station": "Z9046",
        "nrays": 3264,
        "nsweeps": 9,
        "shape": (361, 920),
        "fixed_angle_head": [0.5, 1.5, 2.4000000953674316, 3.4000000953674316, 4.300000190734863],
        "dBZ": {"nanmin": -22.5, "nanmax": 65.5, "nanmean": 15.491404386221905, "nancount": 140976},
        "V": {"nanmin": -27.5, "nanmax": 27.5, "nanmean": -0.252108468483785, "nancount": 167427},
    },
    "Z_RADR_I_Z9047_20260317065837_O_DOR_SAD_CAP_FMT.bin.bz2": {
        "station": "Z9047",
        "nrays": 3263,
        "nsweeps": 9,
        "shape": (361, 920),
        "fixed_angle_head": [0.5, 1.5, 2.4000000953674316, 3.4000000953674316, 4.300000190734863],
        "dBZ": {"nanmin": -30.5, "nanmax": 75.0, "nanmean": 11.793126100538071, "nancount": 148687},
        "V": {"nanmin": -28.0, "nanmax": 28.0, "nanmean": -1.646808578950748, "nancount": 207536},
    },
    "Z_RADR_I_Z9250_20260317065944_O_DOR_SAD_CAP_FMT.bin.bz2": {
        "station": "Z9250",
        "nrays": 3260,
        "nsweeps": 9,
        "shape": (359, 920),
        "fixed_angle_head": [0.5, 1.5, 2.4000000953674316, 3.4000000953674316, 4.300000190734863],
        "dBZ": {"nanmin": -30.5, "nanmax": 58.0, "nanmean": 18.826319487671036, "nancount": 68417},
        "V": {"nanmin": -26.5, "nanmax": 26.5, "nanmean": 0.10244883060555147, "nancount": 143987},
    },
}

SAB_BASELINES = {
    "Z_RADR_I_Z9517_20170811180500_O_DOR_SA_CAP.bin.bz2": {
        "station_dir": ("SAB", "20170812"),
        "nrays": 3313,
        "nsweeps": 9,
        "shape": (369, 920),
        "fixed_angle_head": [0.5, 1.45, 2.4, 3.35, 4.3],
        "dBZ": {"nanmin": -31.5, "nanmax": 59.5, "nanmean": 11.112748250030059, "nancount": 148191},
        "V": {"nanmin": -27.0, "nanmax": 27.0, "nanmean": -0.4214563990533406, "nancount": 202155},
    },
    "Z_RADR_I_Z9517_20170811185100_O_DOR_SA_CAP.bin.bz2": {
        "station_dir": ("SAB", "20170812"),
        "nrays": 3314,
        "nsweeps": 9,
        "shape": (369, 920),
        "fixed_angle_head": [0.5, 1.45, 2.4, 3.35, 4.3],
        "dBZ": {"nanmin": -31.0, "nanmax": 62.5, "nanmean": 11.399538805676729, "nancount": 151273},
        "V": {"nanmin": -27.0, "nanmax": 27.0, "nanmean": -0.08793493210237602, "nancount": 212966},
    },
    "Z_RADR_I_Z9517_20170811210500_O_DOR_SA_CAP.bin.bz2": {
        "station_dir": ("SAB", "20170812"),
        "nrays": 3313,
        "nsweeps": 9,
        "shape": (369, 920),
        "fixed_angle_head": [0.5, 1.45, 2.4, 3.35, 4.3],
        "dBZ": {"nanmin": -29.0, "nanmax": 64.0, "nanmean": 15.995291028050078, "nancount": 162796},
        "V": {"nanmin": -27.0, "nanmax": 27.0, "nanmean": 0.0913032333662895, "nancount": 202935},
    },
}


class WSR98DRegressionTests(unittest.TestCase):
    def test_wsr98d_sample_regressions(self):
        from pycwr.io import read_auto

        root = Path(__file__).resolve().parents[2]
        for filename, expected in WSR98D_BASELINES.items():
            sample = root / expected["station"] / filename
            if not sample.exists():
                self.skipTest(f"sample radar file is not available: {sample}")

            with self.subTest(sample=filename):
                prd = read_auto(str(sample))
                field0 = prd.fields[0]

                self.assertEqual(int(prd.nrays), expected["nrays"])
                self.assertEqual(int(prd.nsweeps), expected["nsweeps"])
                self.assertEqual(field0["dBZ"].shape, expected["shape"])
                np.testing.assert_allclose(
                    prd.scan_info.fixed_angle.values[:5],
                    expected["fixed_angle_head"],
                    rtol=0.0,
                    atol=1e-6,
                )

                for field_name in ("dBZ", "V"):
                    arr = field0[field_name].values
                    field_expected = expected[field_name]
                    self.assertAlmostEqual(float(np.nanmin(arr)), field_expected["nanmin"])
                    self.assertAlmostEqual(float(np.nanmax(arr)), field_expected["nanmax"])
                    self.assertAlmostEqual(float(np.nanmean(arr)), field_expected["nanmean"])
                    self.assertEqual(int(np.isnan(arr).sum()), field_expected["nancount"])


class SABRegressionTests(unittest.TestCase):
    def test_sab_sample_regressions(self):
        from pycwr.io import read_auto

        root = Path(__file__).resolve().parents[2]
        for filename, expected in SAB_BASELINES.items():
            sample = root.joinpath(*expected["station_dir"], filename)
            if not sample.exists():
                self.skipTest(f"sample radar file is not available: {sample}")

            with self.subTest(sample=filename):
                prd = read_auto(str(sample))
                field0 = prd.fields[0]

                self.assertEqual(int(prd.nrays), expected["nrays"])
                self.assertEqual(int(prd.nsweeps), expected["nsweeps"])
                self.assertEqual(field0["dBZ"].shape, expected["shape"])
                np.testing.assert_allclose(
                    prd.scan_info.fixed_angle.values[:5],
                    expected["fixed_angle_head"],
                    rtol=0.0,
                    atol=1e-6,
                )

                for field_name in ("dBZ", "V"):
                    arr = field0[field_name].values
                    field_expected = expected[field_name]
                    self.assertAlmostEqual(float(np.nanmin(arr)), field_expected["nanmin"])
                    self.assertAlmostEqual(float(np.nanmax(arr)), field_expected["nanmax"])
                    self.assertAlmostEqual(float(np.nanmean(arr)), field_expected["nanmean"])
                    self.assertEqual(int(np.isnan(arr).sum()), field_expected["nancount"])


if __name__ == "__main__":
    unittest.main()
