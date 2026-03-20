import unittest

import numpy as np


def _angle_diff_deg(a, b):
    return ((a - b + 180.0) % 360.0) - 180.0


class GeometryRoundTripTests(unittest.TestCase):
    def test_cartesian_xyz_to_antenna_random_round_trip(self):
        from pycwr.core.transforms import antenna_to_cartesian_cwr, cartesian_xyz_to_antenna

        rng = np.random.default_rng(20260319)
        ranges = rng.uniform(500.0, 460000.0, size=2048)
        azimuths = rng.uniform(0.0, 360.0, size=2048)
        elevations = rng.uniform(0.1, 25.0, size=2048)
        height = 183.0

        x, y, z = antenna_to_cartesian_cwr(ranges, azimuths, elevations, height)
        az2, r2, el2 = cartesian_xyz_to_antenna(x, y, z, height)

        self.assertTrue(np.allclose(r2, ranges, rtol=1e-7, atol=5e-5))
        self.assertTrue(np.allclose(el2, elevations, rtol=1e-7, atol=1e-7))
        self.assertTrue(np.all(np.abs(_angle_diff_deg(az2, azimuths)) < 1e-7))

    def test_default_effective_radius_keeps_legacy_transform_values(self):
        from pycwr.core.transforms import (
            DEFAULT_EFFECTIVE_EARTH_RADIUS,
            antenna_to_cartesian_cwr,
            cartesian_xy_elevation_to_range_z,
        )

        ranges = np.array([1000.0, 5000.0, 30000.0, 80000.0])
        azimuths = np.array([0.0, 45.0, 180.0, 315.0])
        elevations = np.array([0.5, 1.5, 3.0, 10.0])
        height = 123.0

        x0, y0, z0 = antenna_to_cartesian_cwr(ranges, azimuths, elevations, height)
        x1, y1, z1 = antenna_to_cartesian_cwr(
            ranges,
            azimuths,
            elevations,
            height,
            effective_earth_radius=DEFAULT_EFFECTIVE_EARTH_RADIUS,
        )
        az0, r0, z_back0 = cartesian_xy_elevation_to_range_z(x0, y0, elevations, height)
        az1, r1, z_back1 = cartesian_xy_elevation_to_range_z(
            x0,
            y0,
            elevations,
            height,
            effective_earth_radius=DEFAULT_EFFECTIVE_EARTH_RADIUS,
        )

        np.testing.assert_array_equal(x0, x1)
        np.testing.assert_array_equal(y0, y1)
        np.testing.assert_array_equal(z0, z1)
        np.testing.assert_array_equal(az0, az1)
        np.testing.assert_array_equal(r0, r1)
        np.testing.assert_array_equal(z_back0, z_back1)

    def test_generic_antenna_transform_delegates_to_cwr_zero_height(self):
        from pycwr.core.transforms import antenna_to_cartesian, antenna_to_cartesian_cwr

        ranges = np.array([1000.0, 5000.0, 30000.0, 80000.0])
        azimuths = np.array([0.0, 45.0, 180.0, 315.0])
        elevations = np.array([0.5, 1.5, 3.0, 10.0])

        x0, y0, z0 = antenna_to_cartesian(ranges, azimuths, elevations)
        x1, y1, z1 = antenna_to_cartesian_cwr(ranges, azimuths, elevations, 0.0)

        np.testing.assert_allclose(x0, x1, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(y0, y1, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(z0, z1, rtol=0.0, atol=0.0)

    def test_generic_vector_transform_delegates_to_cwr_zero_height(self):
        from pycwr.core.transforms import antenna_vectors_to_cartesian, antenna_vectors_to_cartesian_cwr

        ranges = np.array([1000.0, 2000.0, 3000.0])
        azimuths = np.array([0.0, 90.0])
        elevations = np.array([1.0, 2.0])

        x0, y0, z0 = antenna_vectors_to_cartesian(ranges, azimuths, elevations)
        x1, y1, z1 = antenna_vectors_to_cartesian_cwr(ranges, azimuths, elevations, h=0.0)

        np.testing.assert_allclose(x0, x1, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(y0, y1, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(z0, z1, rtol=0.0, atol=0.0)

    def test_cartesian_xy_elevation_to_range_z_round_trip(self):
        from pycwr.core.transforms import (
            antenna_to_cartesian_cwr,
            cartesian_xy_elevation_to_range_z,
        )

        ranges = np.array([1000.0, 5000.0, 30000.0, 80000.0])
        azimuths = np.array([0.0, 45.0, 180.0, 315.0])
        elevations = np.array([0.5, 1.5, 3.0, 10.0])
        height = 123.0

        x, y, z = antenna_to_cartesian_cwr(ranges, azimuths, elevations, height)
        az2, r2, z2 = cartesian_xy_elevation_to_range_z(x, y, elevations, height)

        self.assertTrue(np.allclose(r2, ranges, rtol=1e-7, atol=1e-5))
        self.assertTrue(np.allclose(z2, z, rtol=1e-7, atol=1e-5))
        self.assertTrue(np.all(np.abs(_angle_diff_deg(az2, azimuths)) < 1e-7))

    def test_custom_effective_radius_round_trip(self):
        from pycwr.core.transforms import antenna_to_cartesian_cwr, cartesian_xyz_to_antenna

        ranges = np.array([1000.0, 7000.0, 25000.0, 100000.0])
        azimuths = np.array([10.0, 90.0, 225.0, 359.0])
        elevations = np.array([0.5, 1.0, 2.5, 8.0])
        height = 321.0
        effective_earth_radius = 8100000.0

        x, y, z = antenna_to_cartesian_cwr(
            ranges,
            azimuths,
            elevations,
            height,
            effective_earth_radius=effective_earth_radius,
        )
        az2, r2, el2 = cartesian_xyz_to_antenna(
            x,
            y,
            z,
            height,
            effective_earth_radius=effective_earth_radius,
        )

        self.assertTrue(np.allclose(r2, ranges, rtol=1e-7, atol=1e-5))
        self.assertTrue(np.allclose(el2, elevations, rtol=1e-7, atol=1e-7))
        self.assertTrue(np.all(np.abs(_angle_diff_deg(az2, azimuths)) < 1e-7))

    def test_cartesian_xyz_to_antenna_round_trip(self):
        from pycwr.core.transforms import antenna_to_cartesian_cwr, cartesian_xyz_to_antenna

        ranges = np.array([1000.0, 7000.0, 25000.0, 100000.0])
        azimuths = np.array([10.0, 90.0, 225.0, 359.0])
        elevations = np.array([0.5, 1.0, 2.5, 8.0])
        height = 321.0

        x, y, z = antenna_to_cartesian_cwr(ranges, azimuths, elevations, height)
        az2, r2, el2 = cartesian_xyz_to_antenna(x, y, z, height)

        self.assertTrue(np.allclose(r2, ranges, rtol=1e-7, atol=1e-5))
        self.assertTrue(np.allclose(el2, elevations, rtol=1e-7, atol=1e-7))
        self.assertTrue(np.all(np.abs(_angle_diff_deg(az2, azimuths)) < 1e-7))

    def test_cython_inverse_matches_python(self):
        from pycwr.core.transforms import cartesian_xy_elevation_to_range_z

        try:
            from pycwr.core.RadarGridC import xye_to_antenna
        except ImportError:
            self.skipTest("RadarGridC extension is not available")

        x = 12345.0
        y = -5432.0
        elevation = 1.2
        height = 150.0

        az_py, r_py, z_py = cartesian_xy_elevation_to_range_z(x, y, elevation, height)
        az_cy, r_cy, z_cy = xye_to_antenna(x, y, elevation, height)

        self.assertAlmostEqual(az_cy, az_py, places=10)
        self.assertTrue(np.isclose(r_cy, r_py, rtol=1e-10, atol=1e-6))
        self.assertTrue(np.isclose(z_cy, z_py, rtol=1e-10, atol=1e-6))

    def test_cython_inverse_matches_python_for_custom_effective_radius(self):
        from pycwr.core.transforms import cartesian_xy_elevation_to_range_z

        try:
            from pycwr.core.RadarGridC import xye_to_antenna
        except ImportError:
            self.skipTest("RadarGridC extension is not available")

        x = 12345.0
        y = -5432.0
        elevation = 1.2
        height = 150.0
        effective_earth_radius = 8100000.0

        az_py, r_py, z_py = cartesian_xy_elevation_to_range_z(
            x,
            y,
            elevation,
            height,
            effective_earth_radius=effective_earth_radius,
        )
        az_cy, r_cy, z_cy = xye_to_antenna(
            x,
            y,
            elevation,
            height,
            effective_earth_radius=effective_earth_radius,
        )

        self.assertAlmostEqual(az_cy, az_py, places=10)
        self.assertTrue(np.isclose(r_cy, r_py, rtol=1e-10, atol=1e-6))
        self.assertTrue(np.isclose(z_cy, z_py, rtol=1e-10, atol=1e-6))

    def test_cython_xy_to_azimuth_quadrants(self):
        try:
            from pycwr.core.RadarGridC import xy_to_azimuth
        except ImportError:
            self.skipTest("RadarGridC extension is not available")

        points = [
            ((0.0, 1.0), 0.0),
            ((1.0, 0.0), 90.0),
            ((0.0, -1.0), 180.0),
            ((-1.0, 0.0), 270.0),
            ((-1.0, -1.0), 225.0),
        ]

        for (x, y), expected in points:
            self.assertAlmostEqual(xy_to_azimuth(x, y), expected, places=12)

    def test_cython_cartesian_inverse_matches_python_across_quadrants(self):
        from pycwr.core.transforms import cartesian_xyz_to_antenna

        try:
            from pycwr.core.RadarGridC import cartesian_to_antenna
        except ImportError:
            self.skipTest("RadarGridC extension is not available")

        coords = [
            (1200.0, 3400.0, 180.0),
            (1200.0, -3400.0, 180.0),
            (-1200.0, 3400.0, 180.0),
            (-1200.0, -3400.0, 180.0),
        ]
        height = 150.0

        for x, y, z in coords:
            az_py, r_py, el_py = cartesian_xyz_to_antenna(x, y, z, height)
            az_cy, r_cy, el_cy = cartesian_to_antenna(x, y, z, height)
            self.assertTrue(abs(_angle_diff_deg(az_cy, az_py)) < 1e-10)
            self.assertTrue(np.isclose(r_cy, r_py, rtol=1e-10, atol=1e-5))
            self.assertTrue(np.isclose(el_cy, el_py, rtol=1e-10, atol=1e-9))

    def test_geographic_aeqd_round_trip(self):
        from pycwr.core.transforms import (
            geographic_to_cartesian_aeqd,
            cartesian_to_geographic_aeqd,
        )

        lon = np.array([118.0, 118.5, 119.0])
        lat = np.array([31.0, 31.5, 32.0])
        x, y = geographic_to_cartesian_aeqd(lon, lat, 118.3, 31.2)
        lon2, lat2 = cartesian_to_geographic_aeqd(x, y, 118.3, 31.2)

        self.assertTrue(np.allclose(lon2, lon, atol=1e-10, rtol=0.0))
        self.assertTrue(np.allclose(lat2, lat, atol=1e-10, rtol=0.0))

    def test_geographic_aeqd_matches_pyproj(self):
        from pycwr.core.transforms import (
            cartesian_to_geographic_aeqd,
            geographic_to_cartesian_aeqd,
        )

        try:
            import pyproj
        except ImportError:
            self.skipTest("pyproj is not available")

        lon_0 = 118.3
        lat_0 = 31.2
        lon = np.array([117.8, 118.0, 118.4, 119.1])
        lat = np.array([30.8, 31.0, 31.6, 32.0])
        proj = pyproj.Proj({"proj": "aeqd", "lon_0": lon_0, "lat_0": lat_0, "R": 6370997.0})

        x_ref, y_ref = proj(lon, lat, inverse=False)
        x, y = geographic_to_cartesian_aeqd(lon, lat, lon_0, lat_0)
        np.testing.assert_allclose(x, x_ref, rtol=0.0, atol=1e-6)
        np.testing.assert_allclose(y, y_ref, rtol=0.0, atol=1e-6)

        lon2, lat2 = cartesian_to_geographic_aeqd(x_ref, y_ref, lon_0, lat_0)
        lon_ref, lat_ref = proj(x_ref, y_ref, inverse=True)
        np.testing.assert_allclose(lon2, lon_ref, rtol=0.0, atol=1e-10)
        np.testing.assert_allclose(lat2, lat_ref, rtol=0.0, atol=1e-10)

    def test_geographic_aeqd_handles_projection_center(self):
        from pycwr.core.transforms import (
            cartesian_to_geographic_aeqd,
            geographic_to_cartesian_aeqd,
        )

        lon_0 = 118.3
        lat_0 = 31.2
        x, y = geographic_to_cartesian_aeqd(lon_0, lat_0, lon_0, lat_0)
        lon, lat = cartesian_to_geographic_aeqd(0.0, 0.0, lon_0, lat_0)

        self.assertTrue(np.allclose(np.asarray(x), 0.0, atol=1e-12))
        self.assertTrue(np.allclose(np.asarray(y), 0.0, atol=1e-12))
        self.assertTrue(np.allclose(np.asarray(lon), lon_0, atol=1e-12))
        self.assertTrue(np.allclose(np.asarray(lat), lat_0, atol=1e-12))

    def test_rotation_transforms_preserve_range_norm(self):
        from pycwr.core.transforms import (
            antenna_to_cartesian_track_relative,
            antenna_to_cartesian_earth_relative,
            antenna_to_cartesian_aircraft_relative,
        )

        ranges_km = np.array([1.0, 5.0, 10.0])
        rot = np.array([20.0, 45.0, 90.0])
        roll = np.array([2.0, -3.0, 5.0])
        drift = np.array([10.0, 15.0, -20.0])
        heading = np.array([30.0, 60.0, 120.0])
        tilt = np.array([1.0, 4.0, 8.0])
        pitch = np.array([3.0, -2.0, 6.0])

        x, y, z = antenna_to_cartesian_track_relative(ranges_km, rot, roll, drift, tilt, pitch)
        self.assertTrue(np.allclose(np.sqrt(x ** 2 + y ** 2 + z ** 2), ranges_km * 1000.0, rtol=1e-12, atol=1e-9))

        x, y, z = antenna_to_cartesian_earth_relative(ranges_km, rot, roll, heading, tilt, pitch)
        self.assertTrue(np.allclose(np.sqrt(x ** 2 + y ** 2 + z ** 2), ranges_km * 1000.0, rtol=1e-12, atol=1e-9))

        x, y, z = antenna_to_cartesian_aircraft_relative(ranges_km, rot, tilt)
        self.assertTrue(np.allclose(np.sqrt(x ** 2 + y ** 2 + z ** 2), ranges_km * 1000.0, rtol=1e-12, atol=1e-9))

    def test_interpolate_edges_extrapolates_last_bin_correctly(self):
        from pycwr.core.transforms import _interpolate_range_edges, _interpolate_elevation_edges

        ranges = np.array([1000.0, 2000.0, 3000.0])
        elevations = np.array([0.5, 1.5, 2.5])

        self.assertTrue(np.allclose(_interpolate_range_edges(ranges), np.array([500.0, 1500.0, 2500.0, 3500.0])))
        self.assertTrue(np.allclose(_interpolate_elevation_edges(elevations), np.array([0.0, 1.0, 2.0, 3.0])))

    def test_prd_stores_and_uses_custom_effective_radius(self):
        from pycwr.core.NRadar import PRD
        from pycwr.core.transforms import antenna_vectors_to_cartesian_cwr

        effective_earth_radius = 8100000.0
        fields = {"dBZ": np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)}
        time = np.array(["2026-03-19T00:00:00", "2026-03-19T00:00:01"], dtype="datetime64[ns]")
        ranges = np.array([1000.0, 2000.0], dtype=np.float64)
        azimuth = np.array([0.0, 90.0], dtype=np.float64)
        elevation = np.array([1.0, 1.0], dtype=np.float64)

        prd = PRD(
            fields=fields,
            scan_type="ppi",
            time=time,
            range=ranges,
            azimuth=azimuth,
            elevation=elevation,
            latitude=31.0,
            longitude=118.0,
            altitude=100.0,
            sweep_start_ray_index=np.array([0]),
            sweep_end_ray_index=np.array([1]),
            fixed_angle=np.array([1.0]),
            bins_per_sweep=np.array([2]),
            nyquist_velocity=np.array([10.0]),
            frequency=9.4,
            unambiguous_range=np.array([100000.0]),
            nrays=2,
            nsweeps=1,
            sitename="TEST",
            effective_earth_radius=effective_earth_radius,
        )

        x, y, z = antenna_vectors_to_cartesian_cwr(
            ranges,
            azimuth,
            elevation,
            100.0,
            effective_earth_radius=effective_earth_radius,
        )

        self.assertEqual(prd.effective_earth_radius, effective_earth_radius)
        self.assertEqual(float(prd.scan_info["effective_earth_radius"].values), effective_earth_radius)
        self.assertTrue(np.allclose(prd.fields[0]["x"].values, x))
        self.assertTrue(np.allclose(prd.fields[0]["y"].values, y))
        self.assertTrue(np.allclose(prd.fields[0]["z"].values, z))


if __name__ == "__main__":
    unittest.main()
