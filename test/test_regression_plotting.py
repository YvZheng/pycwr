"""Regression coverage for public and legacy plotting entry points."""

import unittest
import warnings

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import test_examples_sections as section_examples


class PlottingRegressionTests(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    @staticmethod
    def _radar():
        return section_examples.SectionExtractionTests()._build_ppi_prd()

    @staticmethod
    def _endpoints():
        from pycwr.core.transforms import cartesian_to_geographic_aeqd

        lon, lat = cartesian_to_geographic_aeqd(
            np.array([1000.0, 5000.0]), np.array([2000.0, 2000.0]), 120.0, 30.0
        )
        return (float(lon[0]), float(lat[0])), (float(lon[1]), float(lat[1]))

    @staticmethod
    def _map_dependencies():
        try:
            import cartopy.crs as ccrs
        except ImportError as exc:
            raise unittest.SkipTest("cartopy is unavailable") from exc
        from pycwr.draw._plot_core import MapOptions

        return ccrs, MapOptions(
            add_ocean=False, add_land=False, add_lakes=False,
            add_rivers=False, add_province_outline=False,
        )

    def test_easy_rhi_forwards_range_mode_and_resolves_field_alias(self):
        from pycwr.draw import plot_rhi

        radar = self._radar()
        for range_mode in ("native", "aligned"):
            with self.subTest(range_mode=range_mode):
                result = plot_rhi(radar, 45.0, field="reflectivity", range_mode=range_mode)
                expected = radar.extract_rhi(45.0, field_name="dBZ", range_mode=range_mode)
                np.testing.assert_allclose(result.artist.get_array(), expected["dBZ"].values)
                result.fig.canvas.draw()

    def test_easy_rhi_accepts_native_rhi_without_target_azimuth(self):
        from pycwr.draw import plot_rhi

        radar = section_examples.SectionExtractionTests()._build_rhi_prd()
        result = plot_rhi(radar)
        np.testing.assert_allclose(result.artist.get_array(), radar.fields[0]["dBZ"].values)
        result.fig.canvas.draw()

    def test_section_entry_points_resolve_metadata_field_alias(self):
        from pycwr.draw import plot_section, plot_section_lonlat
        from pycwr.draw.VerticalSectionPlot import VerticalSection

        radar = self._radar()
        start_lonlat, end_lonlat = self._endpoints()
        calls = (
            lambda: plot_section(radar, start=(1.0, 2.0), end=(5.0, 2.0), field="reflectivity").artist,
            lambda: plot_section_lonlat(radar, start_lonlat, end_lonlat, field="reflectivity").artist,
            lambda: VerticalSection(radar).section((1.0, 2.0), (5.0, 2.0), "reflectivity"),
            lambda: VerticalSection(radar).RHI(45.0, "reflectivity"),
        )
        for index, call in enumerate(calls):
            with self.subTest(entry_point=index):
                mesh = call()
                self.assertGreater(np.ma.count(mesh.get_array()), 0)
                mesh.figure.canvas.draw()

        fig, ax = plt.subplots()
        mesh = VerticalSection.GUI_section(
            fig, ax, None, radar, (1000.0, 2000.0), (5000.0, 2000.0), "reflectivity"
        )
        self.assertGreater(np.ma.count(mesh.get_array()), 0)
        fig.canvas.draw()

    def test_cartesian_sections_validate_units(self):
        from pycwr.draw import plot_section
        from pycwr.draw.VerticalSectionPlot import VerticalSection

        radar = self._radar()
        calls = (
            lambda: plot_section(radar, start=(1000.0, 2000.0), end=(5000.0, 2000.0), point_units="cm"),
            lambda: VerticalSection(radar).section((1000.0, 2000.0), (5000.0, 2000.0), "dBZ", point_units="cm"),
        )
        for call in calls:
            with self.assertRaisesRegex(ValueError, "point_units"):
                call()

    def test_legacy_lonlat_section_updates_labels_and_keeps_orientation(self):
        from pycwr.core.transforms import cartesian_to_geographic_aeqd
        from pycwr.draw.VerticalSectionPlot import VerticalSection

        start_lonlat, end_lonlat = self._endpoints()
        for orientation in ("vertical", "horizontal"):
            with self.subTest(orientation=orientation), warnings.catch_warnings():
                warnings.simplefilter("error")
                mesh = VerticalSection(self._radar()).section_map(
                    start_lonlat, end_lonlat, "reflectivity", orient=orientation
                )
                self.assertEqual(mesh.pycwr_colorbar.orientation, orientation)
                ax = mesh.axes
                ax.set_xlim(1.0, 3.0)
                ax.set_xticks([1.0, 2.0, 3.0])
                mesh.figure.canvas.draw()
                lon, lat = cartesian_to_geographic_aeqd(
                    np.array([2000.0, 3000.0, 4000.0]), np.full(3, 2000.0), 120.0, 30.0
                )
                self.assertEqual(
                    [tick.get_text() for tick in ax.get_xticklabels()],
                    ["(%.2f, %.2f)" % pair for pair in zip(lon, lat)],
                )
                self.assertIn("Longitude", ax.get_xlabel())

    def test_cartesian_products_support_native_and_aligned_ranges(self):
        from pycwr.draw.RadarPlot import Graph

        for range_mode in (None, "native", "aligned"):
            for product in ("CR", "CAPPI_200"):
                with self.subTest(range_mode=range_mode, product=product):
                    radar = self._radar()
                    fig, ax = plt.subplots()
                    graph = Graph(radar)
                    if product == "CR":
                        mesh = graph.plot_crf(ax, range_mode=range_mode, cbar=False)
                    else:
                        mesh = graph.plot_cappi(ax, level_height=200, range_mode=range_mode, cbar=False)
                    product_key = product if range_mode == "aligned" else product + "_native"
                    np.testing.assert_allclose(
                        mesh.get_array().filled(np.nan), radar.product[product_key].values,
                        equal_nan=True,
                    )
                    fig.canvas.draw()

    def test_rectangular_1d_cartesian_coordinates_follow_matplotlib_axes(self):
        from pycwr.draw.RadarPlot import plot_xy

        data = np.arange(6).reshape(2, 3)
        for use_edges in (False, True):
            with self.subTest(use_edges=use_edges):
                x = [-500.0, 500.0, 1500.0, 2500.0] if use_edges else [0.0, 1000.0, 2000.0]
                y = [-500.0, 500.0, 1500.0] if use_edges else [0.0, 1000.0]
                fig, ax = plt.subplots()
                mesh = plot_xy(ax, x, y, data, cbar=False)
                np.testing.assert_allclose(mesh.get_coordinates()[0, :, 0], [-0.5, 0.5, 1.5, 2.5])
                np.testing.assert_allclose(mesh.get_coordinates()[:, 0, 1], [-0.5, 0.5, 1.5])
                np.testing.assert_array_equal(mesh.get_array(), data)
                fig.canvas.draw()

    def test_legacy_color_limits_are_inferred_from_data(self):
        from pycwr.draw.SingleRadarPlot import RadarGraph
        from pycwr.draw.VerticalSectionPlot import VerticalSection

        x, y = np.meshgrid([0.0, 1000.0], [0.0, 1000.0])
        data = np.array([[20.0, 40.0], [60.0, 80.0]])
        mesh = RadarGraph.simple_plot_ppi_xy(x, y, data, continuously=True)
        self.assertEqual(mesh.get_clim(), (20.0, 80.0))
        fig, ax = plt.subplots()
        mesh = VerticalSection.SectionPlot_VCS(
            fig, ax, None, [x, x], [y, y + 1000.0], [data, data + 100.0], continuously=True
        )
        self.assertEqual(mesh.get_clim(), (20.0, 180.0))
        fig.canvas.draw()

    def test_auto_color_limits_handle_nonfinite_and_constant_data(self):
        from pycwr.draw.SingleRadarPlot import RadarGraph
        from pycwr.draw._plot_core import resolve_field_range

        x, y = np.meshgrid([0.0, 1000.0], [0.0, 1000.0])
        for data, limits in (
            (np.full((2, 2), np.nan), (0.0, 1.0)),
            (np.array([[np.inf, np.nan], [10.0, 20.0]]), (10.0, 20.0)),
            (np.full((2, 2), 20.0), (19.0, 21.0)),
        ):
            with warnings.catch_warnings():
                warnings.simplefilter("error")
                mesh = RadarGraph.simple_plot_ppi_xy(x, y, data, continuously=True)
                self.assertEqual(mesh.get_clim(), limits)
                self.assertEqual(resolve_field_range(None, 0, "custom_field", data), limits)
                mesh.figure.canvas.draw()

    def test_map_projection_support_and_requested_extent(self):
        from pycwr.draw import plot_ppi_map

        ccrs, options = self._map_dependencies()
        extent = (119.9, 120.1, 29.9, 30.1)
        projections = (
            ccrs.PlateCarree(), ccrs.Mercator(),
            ccrs.LambertConformal(central_longitude=120, central_latitude=30),
            ccrs.AzimuthalEquidistant(central_longitude=120, central_latitude=30),
        )
        for projection in projections:
            with self.subTest(projection=type(projection).__name__):
                result = plot_ppi_map(
                    self._radar(), projection=projection, map_options=options, extend=extent, cbar=False
                )
                result.fig.canvas.draw()
                if isinstance(projection, (ccrs.PlateCarree, ccrs.Mercator)):
                    np.testing.assert_allclose(result.ax.get_extent(ccrs.PlateCarree()), extent, atol=1e-6)
                self.assertIsNotNone(result.artist)

    def test_legacy_map_accepts_numpy_station_point_and_infers_color_range(self):
        from pycwr.draw.SingleRadarPlotMap import RadarGraphMap

        _, options = self._map_dependencies()
        data = np.array([[20.0, 40.0], [60.0, 80.0]])
        mesh = RadarGraphMap.simple_plot_ppi_map(
            data, _range=np.array([1000.0, 2000.0]), azimuth=np.array([0.0, 90.0]),
            elevation=np.array([1.0, 1.0]), main_point=np.array([120.0, 30.0]),
            map_options=options, continuously=True,
        )
        self.assertEqual(mesh.get_clim(), (20.0, 80.0))
        mesh.figure.canvas.draw()


if __name__ == "__main__":
    unittest.main()
