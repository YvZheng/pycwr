"""Sample-independent contracts for every public plotting entry point."""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest import mock
import warnings

import matplotlib
import numpy as np
import xarray as xr

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import test_examples_sections as section_examples
import test_examples_wind as wind_examples


class PublicPlotContracts(unittest.TestCase):
    def setUp(self):
        self.radar = section_examples.SectionExtractionTests()._build_ppi_prd()
        self.start = (1.0, 2.0)
        self.end = (5.0, 2.0)
        from pycwr.core.transforms import cartesian_to_geographic_aeqd
        lon, lat = cartesian_to_geographic_aeqd([1000.0, 5000.0], [2000.0, 2000.0], 120.0, 30.0)
        self.start_lonlat = (float(lon[0]), float(lat[0]))
        self.end_lonlat = (float(lon[1]), float(lat[1]))

    def tearDown(self):
        plt.close("all")

    @staticmethod
    def _map_options():
        try:
            import cartopy.crs as ccrs
        except ImportError as exc:
            raise unittest.SkipTest("Cartopy is unavailable") from exc
        from pycwr.draw._plot_core import MapOptions
        return ccrs, MapOptions(
            add_ocean=False, add_land=False, add_lakes=False,
            add_rivers=False, add_province_outline=False,
        )

    def _map_axes(self):
        ccrs, options = self._map_options()
        fig = plt.figure()
        return fig, fig.add_subplot(111, projection=ccrs.PlateCarree()), options

    def _assert_mesh(self, mesh, shape=None):
        self.assertIsNotNone(mesh)
        self.assertIsNotNone(mesh.axes)
        if shape is not None:
            self.assertEqual(mesh.get_array().size, int(np.prod(shape)))
        mesh.figure.canvas.draw()

    @staticmethod
    def _profile(counts=(10.0, 20.0, 30.0)):
        return xr.Dataset(
            {name: ("height", values) for name, values in {
                "u": [3.0, 6.0, 9.0], "v": [4.0, 8.0, 12.0],
                "wind_speed": [5.0, 10.0, 15.0], "sample_count": list(counts),
            }.items()},
            coords={"height": [500.0, 1000.0, 1500.0]},
        )

    @staticmethod
    def _vvp(values=(3.0, 6.0, 9.0)):
        values = np.asarray(values, dtype=float)
        return xr.Dataset({name: ("point", value) for name, value in {
            "u": values, "v": values, "wind_speed": np.hypot(values, values),
            "x": np.arange(values.size) * 1000.0, "y": np.arange(values.size) * 500.0,
        }.items()})

    def test_easy_cartesian_calls_reuse_axes_save_and_show(self):
        from pycwr.draw import plot_ppi, plot_rhi, plot_section, plot_section_lonlat
        cases = (
            (plot_ppi, {"field": "ZDR", "sweep": 1}),
            (plot_rhi, {"azimuth": 45.0, "range_mode": "aligned"}),
            (plot_section, {"start": self.start, "end": self.end}),
            (plot_section_lonlat, {"start_lonlat": self.start_lonlat, "end_lonlat": self.end_lonlat}),
        )
        with TemporaryDirectory() as directory:
            for function, kwargs in cases:
                with self.subTest(entry=function.__name__):
                    fig, ax = plt.subplots()
                    destination = Path(directory) / function.__name__ / "plot.png"
                    with mock.patch("matplotlib.pyplot.show") as show:
                        result = function(self.radar, ax=ax, save=destination, show=True, title="Contract", **kwargs)
                    self.assertIs(result.fig, fig)
                    self.assertIs(result.ax, ax)
                    self.assertEqual(ax.get_title(), "Contract")
                    self.assertGreater(destination.stat().st_size, 0)
                    show.assert_called_once_with()
                    self._assert_mesh(result.artist)
                    plt.close(fig)

    def test_easy_wind_calls_accept_precomputed_datasets_and_prd(self):
        from pycwr.draw import plot_vvp, plot_wind_profile
        for source in (self._profile(), wind_examples.WindRetrievalTests._build_uniform_wind_prd()):
            result = plot_wind_profile(source, show_barbs=False)
            self._assert_mesh(result.artist)
            self.assertEqual(result.ax.get_ylabel(), "Height (km)")
        for source, kwargs in (
            (self._vvp(), {}),
            (wind_examples.WindRetrievalTests._build_uniform_wind_prd(), {"az_num": 9, "bin_num": 5}),
        ):
            result = plot_vvp(source, **kwargs)
            self._assert_mesh(result.artist)

    def test_easy_dispatch_all_kinds(self):
        from pycwr.draw import plot
        cases = (
            (self.radar, "ppi", {}),
            (self.radar, "rhi", {"azimuth": 45.0}),
            (self.radar, "section", {"start": self.start, "end": self.end}),
            (self.radar, "section_lonlat", {"start_lonlat": self.start_lonlat, "end_lonlat": self.end_lonlat}),
            (self._profile(), "wind_profile", {}),
            (self._vvp(), "vvp", {}),
        )
        for source, kind, kwargs in cases:
            with self.subTest(kind=kind):
                result = plot(source, kind=kind, **kwargs)
                self._assert_mesh(result.artist)
                plt.close(result.fig)
        with self.assertRaisesRegex(ValueError, "unsupported kind"):
            plot(self.radar, kind="unknown")

    def test_easy_plotter_all_cartesian_methods(self):
        from pycwr.draw import EasyRadarPlotter
        plotter = EasyRadarPlotter(self.radar)
        for name, kwargs in (
            ("ppi", {}), ("rhi", {"azimuth": 45.0}),
            ("section", {"start": self.start, "end": self.end}),
            ("section_lonlat", {"start_lonlat": self.start_lonlat, "end_lonlat": self.end_lonlat}),
            ("quicklook", {"kind": "ppi"}),
        ):
            with self.subTest(method=name):
                result = getattr(plotter, name)(**kwargs)
                self._assert_mesh(result.artist)
                plt.close(result.fig)
        self._assert_mesh(EasyRadarPlotter(self._profile()).wind_profile().artist)
        self._assert_mesh(EasyRadarPlotter(self._vvp()).vvp().artist)

    def test_easy_all_mapped_wrappers_render_without_downloads(self):
        from pycwr.draw import EasyRadarPlotter, plot, plot_ppi_map
        _, options = self._map_options()
        calls = (
            lambda: plot_ppi_map(self.radar, map_options=options),
            lambda: plot(self.radar, kind="ppi_map", map_options=options),
            lambda: EasyRadarPlotter(self.radar).ppi_map(map_options=options),
        )
        with mock.patch("cartopy.io.Downloader._urlopen", side_effect=AssertionError("Unexpected map download")):
            for call in calls:
                result = call()
                self.assertIs(result.ax.pycwr_last_mappable, result.artist)
                self._assert_mesh(result.artist, self.radar.fields[0]["dBZ"].shape)
                plt.close(result.fig)

    def test_graph_scalar_field_and_product_methods(self):
        from pycwr.draw import Graph, RadarDisplay
        self.assertTrue(issubclass(RadarDisplay, Graph))
        graph = Graph(self.radar)
        calls = (
            lambda ax: graph.plot_ppi(ax, 1, "reflectivity", cbar_ticks=[0, 20, 40], cbar_ticklabels=["low", "medium", "high"]),
            lambda ax: graph.plot_rhi(ax, 0, "ZDR", continuously=True),
            lambda ax: graph.plot_vcs(ax, self.start, self.end, "dBZ"),
            lambda ax: graph.plot_crf(ax),
            lambda ax: graph.plot_cappi(ax, level_height=200.0),
        )
        for index, call in enumerate(calls):
            with self.subTest(method=index):
                fig, ax = plt.subplots()
                mesh = call(ax)
                self.assertIs(ax.pycwr_last_mappable, mesh)
                self._assert_mesh(mesh)
                plt.close(fig)

    def test_graph_native_rhi_sweep(self):
        from pycwr.draw import Graph
        radar = section_examples.SectionExtractionTests()._build_rhi_prd()
        fig, ax = plt.subplots()
        mesh = Graph(radar).plot_rhi(ax, 0, "dBZ", height_km=(0.0, 3.0))
        self._assert_mesh(mesh, radar.fields[0]["dBZ"].shape)
        np.testing.assert_allclose(ax.get_ylim(), [0.0, 3.0])

    def test_graph_wind_profile_missing_sample_counts_and_empty_profile(self):
        from pycwr.draw import Graph
        for profile in (self._profile(), self._profile((10.0, np.nan, np.inf)), self._profile().drop_vars("sample_count"), self._profile().isel(height=slice(0, 0))):
            with warnings.catch_warnings():
                warnings.simplefilter("error", RuntimeWarning)
                fig, ax = plt.subplots()
                artist = Graph(None).plot_wind_profile(ax, profile=profile)
                self._assert_mesh(artist)
                if "sample_count" in profile:
                    self.assertEqual(len(ax._pycwr_barb_axes.texts), int(np.isfinite(profile.sample_count).sum()))
                plt.close(fig)

    def test_graph_vvp_background_options_and_missing_calm_winds(self):
        from pycwr.draw import Graph
        for values in ((3.0, 6.0, 9.0), (0.0, 0.0, 0.0), (np.nan, np.nan), ()):
            for background in ("speed", None):
                with self.subTest(values=values, background=background), warnings.catch_warnings():
                    warnings.simplefilter("error", RuntimeWarning)
                    fig, ax = plt.subplots()
                    quiver = Graph(None).plot_vvp(ax, 0, retrieval=self._vvp(values), background=background)
                    self._assert_mesh(quiver)
                    self.assertEqual(quiver.U.size, np.isfinite(values).sum())
                    plt.close(fig)
        radar = wind_examples.WindRetrievalTests._build_uniform_wind_prd()
        fig, ax = plt.subplots()
        quiver = Graph(radar).plot_vvp(ax, 0, az_num=9, bin_num=5, background_field="dBZ", max_range_km=3.0, extent_km=(-3, 3, -3, 3))
        self._assert_mesh(quiver)
        np.testing.assert_allclose(ax.get_xlim(), [-3.0, 3.0])

    def test_cartesian_decorations_and_low_level_plot_functions(self):
        from pycwr.draw.RadarPlot import Graph, add_rings, plot_az_ranges, plot_xy
        fig, ax = plt.subplots()
        graph = Graph(self.radar)
        for rings in ([], np.array([])):
            self.assertIs(add_rings(ax, rings), ax)
            self.assertIs(graph.add_rings(ax, rings), ax)
            self.assertEqual(len(ax.lines), 0)
        graph.add_rings(ax, [1.0, 2.0])
        self.assertEqual(len(ax.lines), 8)
        line = graph.add_lines(ax, self.start, self.end)[0]
        np.testing.assert_array_equal(line.get_xdata(), [1.0, 5.0])
        fig.canvas.draw()
        field = self.radar.fields[0]["dBZ"]
        for call in (
            lambda ax: plot_xy(ax, field.x, field.y, field, cbar=False),
            lambda ax: plot_az_ranges(ax, field.range, field.azimuth, field.elevation, field, cbar=False),
        ):
            fig, ax = plt.subplots()
            self._assert_mesh(call(ax), field.shape)

    def test_vvp_zdr_background_resolves_colormap_without_pyart(self):
        import test_issue60_colormap as colormap_examples

        colormap_examples.Issue60ColormapTests()._run_fresh(colormap_examples._RADAR_SETUP + """
import xarray as xr
from pycwr.draw import Graph
retrieval = xr.Dataset({name: ('point', values) for name, values in {
    'u': [3., 6.], 'v': [4., 8.], 'wind_speed': [5., 10.],
    'x': [1000., 2000.], 'y': [1000., 2000.],
}.items()})
fig, ax = plt.subplots()
Graph(radar).plot_vvp(ax, 0, retrieval=retrieval, background_field='ZDR')
fig.canvas.draw()
assert ax._pycwr_last_mappable.cmap.name == 'copy_pyart_RefDiff'
assert 'pyart' not in __import__('sys').modules
plt.close(fig)
""")

    def test_graphmap_all_methods_without_external_samples(self):
        from pycwr.draw import GraphMap, RadarMapDisplay
        ccrs, options = self._map_options()
        self.assertTrue(issubclass(RadarMapDisplay, GraphMap))
        graph = GraphMap(self.radar, ccrs.PlateCarree())
        with mock.patch("cartopy.io.Downloader._urlopen", side_effect=AssertionError("Unexpected map download")):
            for name in ("plot_ppi_map", "plot_crf_map", "plot_cappi_map"):
                for range_mode in ("native", "aligned"):
                    with self.subTest(method=name, range_mode=range_mode):
                        fig, ax, _ = self._map_axes()
                        args = (ax, 0, "dBZ") if name == "plot_ppi_map" else ((ax, 200.0) if name == "plot_cappi_map" else (ax,))
                        returned = getattr(graph, name)(*args, range_mode=range_mode, map_options=options, extend=(119.95, 120.05, 29.95, 30.05))
                        self.assertIs(returned, ax)
                        self._assert_mesh(ax.pycwr_last_mappable)
                        line = graph.add_lines_map(ax, self.start_lonlat, self.end_lonlat)[0]
                        np.testing.assert_allclose(line.get_xdata(), [self.start_lonlat[0], self.end_lonlat[0]])
                        plt.close(fig)
        fig, ax = plt.subplots()
        self._assert_mesh(graph.plot_vcs_map(ax, self.start_lonlat, self.end_lonlat, "dBZ"))

    def test_legacy_cartesian_all_methods_and_custom_colormap(self):
        from pycwr.draw.SingleRadarPlot import RadarGraph
        radar = self.radar
        field = radar.fields[0]["dBZ"]
        self._assert_mesh(RadarGraph(NuistRadar=radar).plot(0, "dBZ", cmap="plasma", cmap_bins=7, dark=True))
        fig, ax = plt.subplots()
        mesh = RadarGraph.GUI_plot(radar, fig, ax, None, 0, "dBZ", cmap="plasma", cmap_bins=7)
        self.assertEqual(mesh.cmap.name, "plasma")
        self._assert_mesh(mesh)
        self._assert_mesh(RadarGraph.simple_plot_ppi(field))
        self._assert_mesh(RadarGraph.simple_plot_ppi(field.values, _range=field.range, azimuth=field.azimuth, elevation=field.elevation))
        self._assert_mesh(RadarGraph.simple_plot_ppi_xy(field.x, field.y, field))
        fig, ax = plt.subplots()
        self._assert_mesh(RadarGraph.plot_ppi(fig, ax, None, field.x, field.y, field))

    def test_legacy_map_all_methods_and_custom_colormap(self):
        from pycwr.draw.SingleRadarPlotMap import RadarGraphMap, ax_projection
        ccrs, options = self._map_options()
        self.assertIsInstance(ax_projection(None), ccrs.PlateCarree)
        self.assertIsInstance(ax_projection(options), ccrs.PlateCarree)
        field = self.radar.fields[0]["dBZ"]
        with mock.patch("cartopy.io.Downloader._urlopen", side_effect=AssertionError("Unexpected map download")):
            self._assert_mesh(RadarGraphMap(NuistRadar=self.radar).plot(0, "dBZ", map_options=options, cmap="plasma", cmap_bins=7))
            fig, ax, _ = self._map_axes()
            mesh = RadarGraphMap.GUI_plot(self.radar, fig, ax, None, 0, "dBZ", map_options=options, cmap="plasma", cmap_bins=7)
            self.assertEqual(mesh.cmap.name, "plasma")
            self._assert_mesh(mesh)
            self._assert_mesh(RadarGraphMap.simple_plot_ppi_map(field, map_options=options))
            self._assert_mesh(RadarGraphMap.simple_plot_ppi_xy_map(field.x, field.y, field.values, (120.0, 30.0), map_options=options))
            fig, ax, _ = self._map_axes()
            self._assert_mesh(RadarGraphMap.plot_ppi_map(fig, ax, None, field.lon, field.lat, field, map_options=options))

    def test_legacy_vertical_section_all_methods(self):
        from pycwr.draw.VerticalSectionPlot import VerticalSection
        section = VerticalSection(NuistRadar=self.radar)
        self._assert_mesh(section.RHI(45.0, "dBZ"))
        self._assert_mesh(section.section(self.start, self.end, "dBZ"))
        self._assert_mesh(section.section_map(self.start_lonlat, self.end_lonlat, "dBZ"))
        for name, args in (
            ("GUI_section", (self.radar, (1000.0, 2000.0), (5000.0, 2000.0), "dBZ")),
            ("GUI_section_map", (self.radar, self.start_lonlat, self.end_lonlat, "dBZ")),
            ("SectionPlot_VCS_map", (self.start_lonlat, self.end_lonlat, "dBZ", self.radar)),
        ):
            with self.subTest(method=name):
                fig, ax = plt.subplots()
                self._assert_mesh(getattr(section, name)(fig, ax, None, *args))
        fig, ax = plt.subplots()
        meshes = self.radar.get_vcs_data((1000.0, 2000.0), (5000.0, 2000.0), "dBZ")
        self._assert_mesh(section.SectionPlot_VCS(fig, ax, None, *meshes))
        x, y = section.get_points_from_ranges((1.0, 2.0), (5.0, 2.0), np.array([0.0, 2.0, 4.0]))
        np.testing.assert_allclose(x, [1.0, 3.0, 5.0])
        np.testing.assert_allclose(y, [2.0, 2.0, 2.0])
        with self.assertRaisesRegex(ValueError, "distinct"):
            section.get_points_from_ranges((1.0, 2.0), (1.0, 2.0), 0.0)

    def test_invalid_fields_units_axes_and_legacy_arguments(self):
        from pycwr.draw import Graph, GraphMap, plot_section
        from pycwr.draw.SingleRadarPlot import RadarGraph
        from pycwr.draw.SingleRadarPlotMap import RadarGraphMap
        from pycwr.draw.VerticalSectionPlot import VerticalSection
        fig, ax = plt.subplots()
        with self.assertRaises(KeyError):
            Graph(self.radar).plot_ppi(ax, 0, "missing_field")
        with self.assertRaises(ValueError):
            Graph(self.radar).plot_ppi(ax, 0, "dBZ", range_mode="unknown")
        with self.assertRaises(ValueError):
            plot_section(self.radar, start=self.start, end=self.end, point_units="miles")
        with self.assertRaises(ValueError):
            plot_section(self.radar, start=self.start, end=self.start)
        with self.assertRaises(TypeError):
            Graph(self.radar).plot_ppi(None, 0, "dBZ")
        for cls in (RadarGraph, RadarGraphMap, VerticalSection):
            with self.assertRaisesRegex(TypeError, "unexpected keyword"):
                cls(self.radar, unknown=True)
        ccrs, options = self._map_options()
        with self.assertRaises(TypeError):
            GraphMap(self.radar, ccrs.PlateCarree()).plot_ppi_map(ax, 0, "dBZ", map_options=options)
        fig, map_ax, _ = self._map_axes()
        with self.assertRaises(TypeError):
            Graph(self.radar).plot_ppi(map_ax, 0, "dBZ")

    def test_web_palette_continuous_colormap_interpolates_and_renders(self):
        from matplotlib.colors import LinearSegmentedColormap, to_rgba
        from pycwr.GraphicalInterface.web_colors import FieldPalette, PALETTES

        palette = FieldPalette(levels=(0.0, 1.0), labels=("low", "high"), colors=("black", "white"))
        cmap = palette.continuous_cmap("gradient")
        values = np.linspace(0.0, 1.0, 256)
        colors = cmap(values)
        self.assertIsInstance(cmap, LinearSegmentedColormap)
        np.testing.assert_allclose(colors[:, :3], np.repeat(values[:, None], 3, axis=1), atol=1e-12)
        np.testing.assert_allclose(colors[:, 3], 1.0)
        fig, ax = plt.subplots()
        artist = ax.imshow(values[None, :], cmap=cmap, aspect="auto")
        self._assert_mesh(artist)

        for field_name, field_palette in PALETTES.items():
            with self.subTest(field=field_name):
                field_cmap = field_palette.continuous_cmap(field_name)
                rgba = field_cmap(values)
                self.assertEqual(rgba.shape, (256, 4))
                self.assertTrue(np.isfinite(rgba).all())
                self.assertTrue(((rgba >= 0.0) & (rgba <= 1.0)).all())
                np.testing.assert_allclose(rgba[0], to_rgba(field_palette.colors[0]))
                np.testing.assert_allclose(rgba[-1], to_rgba(field_palette.colors[-1]))


if __name__ == "__main__":
    unittest.main()
