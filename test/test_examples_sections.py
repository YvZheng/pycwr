"""Section and RHI extraction examples plus regression coverage."""

import unittest

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


class SectionExtractionTests(unittest.TestCase):
    def _build_ppi_prd(self, azimuths=None, ranges=None, elevations=(1.0, 3.0)):
        from pycwr.core.NRadar import PRD

        azimuths = np.asarray(azimuths if azimuths is not None else np.arange(0.0, 360.0, 45.0), dtype=np.float64)
        ranges = np.asarray(ranges if ranges is not None else np.arange(1000.0, 9000.0, 1000.0), dtype=np.float64)
        nsweeps = len(elevations)
        nrays_per_sweep = azimuths.size
        nrays = nsweeps * nrays_per_sweep
        azimuth = np.tile(azimuths, nsweeps)
        elevation = np.repeat(np.asarray(elevations, dtype=np.float64), nrays_per_sweep)
        base_field = azimuth[:, None] + ranges[None, :] / 1000.0
        fields = {
            "dBZ": base_field.copy(),
            "ZDR": base_field.copy() * 0.1,
        }
        sweep_start = np.arange(nsweeps, dtype=np.int32) * nrays_per_sweep
        sweep_end = sweep_start + nrays_per_sweep - 1
        return PRD(
            fields=fields,
            scan_type="ppi",
            time=np.arange(nrays),
            range=ranges,
            azimuth=azimuth,
            elevation=elevation,
            latitude=30.0,
            longitude=120.0,
            altitude=100.0,
            sweep_start_ray_index=sweep_start,
            sweep_end_ray_index=sweep_end,
            fixed_angle=np.asarray(elevations, dtype=np.float64),
            bins_per_sweep=np.full(nsweeps, ranges.size, dtype=np.int32),
            nyquist_velocity=np.full(nsweeps, 15.0, dtype=np.float64),
            frequency=5.6,
            unambiguous_range=np.full(nsweeps, 150000.0, dtype=np.float64),
            nrays=nrays,
            nsweeps=nsweeps,
            sitename="TEST",
        )

    def _build_rhi_prd(self):
        from pycwr.core.NRadar import PRD

        ranges = np.arange(1000.0, 5000.0, 1000.0)
        elevation = np.array([0.5, 1.5, 2.5], dtype=np.float64)
        azimuth = np.full(elevation.shape, 45.0, dtype=np.float64)
        fields = {
            "dBZ": elevation[:, None] + ranges[None, :] / 1000.0,
        }
        return PRD(
            fields=fields,
            scan_type="rhi",
            time=np.arange(elevation.size),
            range=ranges,
            azimuth=azimuth,
            elevation=elevation,
            latitude=30.0,
            longitude=120.0,
            altitude=100.0,
            sweep_start_ray_index=np.array([0], dtype=np.int32),
            sweep_end_ray_index=np.array([elevation.size - 1], dtype=np.int32),
            fixed_angle=np.array([45.0], dtype=np.float64),
            bins_per_sweep=np.array([ranges.size], dtype=np.int32),
            nyquist_velocity=np.array([15.0], dtype=np.float64),
            frequency=5.6,
            unambiguous_range=np.array([150000.0], dtype=np.float64),
            nrays=elevation.size,
            nsweeps=1,
            sitename="RHI",
        )

    def test_extract_section_matches_linear_analytic_field(self):
        from pycwr.core.transforms import cartesian_to_antenna_cwr

        prd = self._build_ppi_prd()
        section = prd.extract_section(
            (1000.0, 2000.0),
            (5000.0, 2000.0),
            field_name="dBZ",
            point_units="m",
        )

        self.assertEqual(section["dBZ"].shape[0], 2)
        self.assertEqual(section.attrs["interpolation"], "linear")
        for strip_index, elevation in enumerate((1.0, 3.0)):
            azimuth, ranges, _ = cartesian_to_antenna_cwr(
                section["x"].isel(sweep=strip_index).values,
                section["y"].isel(sweep=strip_index).values,
                elevation,
                float(prd.scan_info.altitude.values),
                effective_earth_radius=prd.effective_earth_radius,
            )
            expected = azimuth + ranges / 1000.0
            np.testing.assert_allclose(
                section["dBZ"].isel(sweep=strip_index).values,
                expected,
                atol=1e-6,
                equal_nan=True,
            )

    def test_extract_rhi_handles_azimuth_wraparound(self):
        prd = self._build_ppi_prd(azimuths=np.array([5.0, 90.0, 180.0, 355.0]), elevations=(1.0,))
        sweep = prd.fields[0]
        custom = np.vstack(
            [
                np.full(sweep.range.size, 20.0),
                np.full(sweep.range.size, 40.0),
                np.full(sweep.range.size, 60.0),
                np.full(sweep.range.size, 10.0),
            ]
        )
        prd.fields[0]["dBZ"] = (("time", "range"), custom)
        section = prd.extract_rhi(azimuth=359.0, field_name="dBZ")

        np.testing.assert_allclose(section["dBZ"].isel(sweep=0).values, np.full(sweep.range.size, 14.0))

    def test_extract_rhi_supports_native_rhi_scan(self):
        prd = self._build_rhi_prd()
        section = prd.extract_rhi(field_name="dBZ")

        self.assertEqual(section.attrs["scan_type"], "rhi")
        self.assertEqual(section["dBZ"].shape, (3, 4))
        np.testing.assert_allclose(section["dBZ"].values, prd.fields[0]["dBZ"].values)
        np.testing.assert_array_equal(section["source_sweep"].values, np.zeros(3, dtype=np.int32))

    def test_legacy_section_helpers_remain_compatible(self):
        prd = self._build_ppi_prd()

        mesh_xy, mesh_z, field_data = prd.get_vcs_data((1000.0, 0.0), (6000.0, 0.0), "dBZ")
        rhi_xy, rhi_z, rhi_field = prd.get_RHI_data(45.0, "dBZ")

        self.assertEqual(len(mesh_xy), prd.nsweeps)
        self.assertEqual(len(mesh_z), prd.nsweeps)
        self.assertEqual(len(field_data), prd.nsweeps)
        self.assertEqual(len(rhi_xy), prd.nsweeps)
        self.assertEqual(rhi_xy[0].shape, rhi_z[0].shape)
        self.assertEqual(field_data[0].shape[0], 1)
        self.assertEqual(rhi_field[0].shape[0], 1)


class SectionPlotTests(unittest.TestCase):
    def _build_plot_prd(self):
        return SectionExtractionTests()._build_ppi_prd()

    def test_easy_plot_section_smoke(self):
        from pycwr.draw.easy import plot_section

        prd = self._build_plot_prd()
        result = plot_section(
            prd,
            start=(1.0, 2.0),
            end=(5.0, 2.0),
            field="dBZ",
            show=False,
        )
        try:
            self.assertIsNotNone(result.artist)
        finally:
            plt.close(result.fig)

    def test_easy_plot_rhi_smoke(self):
        from pycwr.draw.easy import plot_rhi

        prd = self._build_plot_prd()
        result = plot_rhi(
            prd,
            azimuth=45.0,
            field="dBZ",
            show=False,
        )
        try:
            self.assertIsNotNone(result.artist)
        finally:
            plt.close(result.fig)


if __name__ == "__main__":
    unittest.main()
