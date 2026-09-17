"""Public IO contracts using small synthetic protocol fixtures, not field samples."""
import bz2
import gzip
import io
from pathlib import Path
import struct
import tempfile
import unittest
from unittest import mock
import zipfile

import numpy as np


def _packed(structure, **values):
    return struct.pack('<' + ''.join(kind for _, kind in structure), *[
        values.get(name, b'' if kind.endswith('s') else 0) for name, kind in structure
    ])


def _place(buffer, offset, payload):
    buffer[offset:offset + len(payload)] = payload


def cc_fixture(gate_counts=(4,), ray_counts=None):
    from pycwr.io.BaseDataProtocol.CCProtocol import dtype_cc as p
    ray_counts = ray_counts or (2,) * len(gate_counts)
    header = bytearray(1024)
    _place(header, 0, _packed(p.BaseDataHeader['RadarHeader1'], cRadarType=b'CINRAD/CC',
        lLongitudeValue=120 * 3600000, lLatitudeValue=30 * 3600000, lHeight=100000,
        ucSYear1=20, ucSYear2=26, ucSMonth=1, ucSDay=1, ucEYear1=20, ucEYear2=26,
        ucEMonth=1, ucEDay=1, ucESecond=10, ucScanMode=100 + len(gate_counts)))
    cuts = np.zeros(30, dtype=p.BaseDataHeader['CutConfigX30'])
    for sweep, (gates, rays) in enumerate(zip(gate_counts, ray_counts)):
        for key, value in dict(usMaxV=1500, usMaxL=15000, usBindWidth=125,
            usBinNumber=gates, usRecordNumber=rays, usAngle=100 + 100 * sweep).items():
            cuts[key][sweep] = value
    _place(header, p.CutSize_pos, cuts.tobytes())
    _place(header, p.HeaderSize2_pos, _packed(p.BaseDataHeader['RadarHeader2'], lWavelength=100000))
    radials = []
    for gates, rays in zip(gate_counts, ray_counts):
        for _ in range(rays):
            data = np.zeros(1, dtype=p.RadialData(gates))
            for field in ('dBZ', 'V', 'W'):
                data[field][0] = np.arange(gates) * 10
                if gates > 1:
                    data[field][0, 1] = -32768
            radials.append(data.tobytes().ljust(p.PerRadialSize, b'\0'))
    return bytes(header) + b''.join(radials)


def sc_fixture(rays=2):
    from pycwr.io.BaseDataProtocol.SCProtocol import dtype_sc as p
    header = bytearray(1024)
    _place(header, 0, _packed(p.BaseDataHeader['RadarSite'], radartype=b'CINRAD/SC',
        longitudevalue=12000, latitudevalue=3000, height=100000))
    _place(header, p.RadarObserationParamPos_1, _packed(p.BaseDataHeader['RadarObserationParam_1'],
        stype=101, syear=2026, smonth=1, sday=1, shour=8))
    _place(header, p.RadarObserationParamPos_2, _packed(p.BaseDataHeader['RadarObserationParam_2'],
        Eyear=2026, Emonth=1, Eday=1, Ehour=8, Esecond=10))
    cuts = np.zeros(30, dtype=p.BaseDataHeader['LayerParamX30'])
    for key, value in dict(MaxV=1600, MaxL=15000, binWidth=5000, binnumber=500,
                           recordnumber=rays, Swangles=100).items():
        cuts[key][0] = value
    _place(header, p.LayerParamPos, cuts.tobytes())
    radials = []
    for ray in range(rays):
        data = np.zeros(500, dtype=p.RadialData())
        data['dBZ'][:4] = [64, 66, 0, 68]
        data['dBT'][:4] = [64, 66, 0, 68]
        data['V'][:4] = [128, 136, 0, 144]
        data['W'][:4] = [16, 32, 0, 48]
        radial = _packed(p.RadialHeader(), sStrAz=int(ray * 65536 / rays),
                         sEndAz=int(ray * 65536 / rays), sStrEl=91, sEndEl=91) + data.tobytes()
        radials.append(radial.ljust(p.PerRadialSize, b'\0'))
    return bytes(header) + b''.join(radials)


def sab_fixture(split=False):
    from pycwr.io.BaseDataProtocol.SABProtocol import dtype_sab as p
    radials = []
    for ray, state in enumerate((3, 2, 0, 4) if split else (3, 4)):
        ref_count = 0 if split and ray >= 2 else 4
        dop_count = 0 if split and ray < 2 else 4
        azimuth = (ray % 2) * 32768 if ray < 2 else (1 - ray % 2) * 32768
        header = _packed(p.RadialHeader(), flag=1, JulianDate=20455, mSends=ray * 1000,
            URange=1500, AZ=azimuth, RadialNumber=ray % 2 + 1, RadialStatus=state,
            El=182, ElNumber=ray // 2 + 1, GateSizeOfReflectivity=250, GateSizeOfDoppler=250,
            GatesNumberOfReflectivity=ref_count, GatesNumberOfDoppler=dop_count,
            PtrOfReflectivity=100, PtrOfVelocity=100 + ref_count,
            PtrOfSpectrumWidth=100 + ref_count + dop_count,
            ResolutionOfVelocity=2, Nyquist=1500)
        data = np.zeros(1, dtype=p.RadialData(ref_count, dop_count))
        if ref_count:
            data['dBZ'][0] = [86, 88, 0, 90] if split and ray == 1 else [66, 68, 0, 70]
        if dop_count:
            data['V'][0] = data['W'][0] = [129, 131, 0, 133]
        radials.append((header + data.tobytes()).ljust(2432, b'\0'))
    return b''.join(radials)


def pa_fixture():
    from pycwr.io.BaseDataProtocol.PAProtocol import dtype_PA as p
    header = bytearray(416)
    _place(header, 0, _packed(p.BaseDataHeader['GenericHeaderBlock'], MagicWord=1297371986, GenericType=16))
    _place(header, 32, _packed(p.BaseDataHeader['SiteConfigurationBlock'], SiteCode=b'SYNTH', SiteName=b'Synthetic PA',
        Latitude=30, Longitude=120, Height=100, Frequency=3000))
    _place(header, 160, _packed(p.BaseDataHeader['TaskConfigurationBlock'], CutNumber=1, ScanType=0))
    cuts = np.zeros(1, dtype=p.BaseDataHeader['CutConfigurationBlock'])
    for key, value in dict(Elevation=1, LogResolution=250, DopplerResolution=250,
                          NyquistSpeed=15, MaximumRange=150000).items():
        cuts[key][0] = value
    radials = []
    for ray, state in enumerate((3, 4)):
        moments = b''.join(_packed(p.RadialData(), DataType=kind, Scale=2, Offset=66,
            BinLength=1, Length=4) + bytes([66, 68, 0, 70]) for kind in (2, 3, 4))
        radials.append(_packed(p.RadialHeader(), RadialState=state, RadialNumber=ray + 1,
            ElevationNumber=1, Azimuth=ray * 180, Elevation=1, Seconds=1767225600 + ray,
            MomentNumber=3, LengthOfData=len(moments)) + moments)
    return bytes(header) + cuts.tobytes() + b''.join(radials)


def wsr98d_fixture(split=False):
    from pycwr.io.BaseDataProtocol.WSR98DProtocol import dtype_98D as p
    header = bytearray(416)
    _place(header, 0, _packed(p.BaseDataHeader['GenericHeaderBlock'], MagicWord=1297371986, GenericType=1))
    _place(header, 32, _packed(p.BaseDataHeader['SiteConfigurationBlock'], SiteCode=b'SYNTH',
        SiteName=b'Synthetic WSR', Latitude=30, Longitude=120, Height=100, Frequency=3000))
    _place(header, 160, _packed(p.BaseDataHeader['TaskConfigurationBlock'], CutNumber=2 if split else 1,
        ScanType=0, VolumeStartTime=1767225600))
    cuts = np.zeros(2 if split else 1, dtype=p.BaseDataHeader['CutConfigurationBlock'])
    for key, value in dict(Elevation=1, LogResolution=250, DopplerResolution=250,
                          NyquistSpeed=15, MaximumRange=150000).items():
        cuts[key] = value
    radials = []
    for ray, state in enumerate((3, 2, 0, 4) if split else (3, 4)):
        kinds = ((2,) if ray < 2 else (3, 4)) if split else (2, 3, 4)
        codes = bytes([86, 88, 0, 90]) if split and ray == 1 else bytes([66, 68, 0, 70])
        moments = b''.join(_packed(p.RadialData(), DataType=kind, Scale=2, Offset=66,
            BinLength=1, Length=4) + codes for kind in kinds)
        azimuth = (ray % 2) * 180 if ray < 2 else (1 - ray % 2) * 180
        radials.append(_packed(p.RadialHeader(), RadialState=state, SequenceNumber=ray + 1,
            RadialNumber=ray % 2 + 1, ElevationNumber=ray // 2 + 1,
            Azimuth=azimuth, Elevation=1, Seconds=1767225600 + ray,
            MomentNumber=len(kinds), LengthOfData=len(moments)) + moments)
    return bytes(header) + cuts.tobytes() + b''.join(radials)


class SyntheticReaderContracts(unittest.TestCase):
    def _check_read(self, family, payload):
        import pycwr.io as radar_io
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / ('synthetic_' + family + '.bin')
            path.write_bytes(payload)
            self.assertEqual(radar_io.radar_format(path), family)
            automatic = radar_io.read_auto(path)
            direct = getattr(radar_io, 'read_' + family)(path)
            self.assertEqual(automatic.nrays, 2)
            self.assertEqual(automatic.nsweeps, 1)
            for field in ('dBZ', 'V', 'W'):
                np.testing.assert_allclose(automatic.fields[0][field], direct.fields[0][field], equal_nan=True)
            return automatic

    def test_cc_minimal_public_reader(self):
        radar = self._check_read('CC', cc_fixture())
        np.testing.assert_allclose(radar.fields[0].dBZ[0], [0, np.nan, 2, 3], equal_nan=True)
        np.testing.assert_allclose(radar.fields[0].azimuth, [0, 180])

    def test_sc_minimal_public_reader(self):
        radar = self._check_read('SC', sc_fixture())
        np.testing.assert_allclose(radar.fields[0].dBZ[0, :4], [0, 1, np.nan, 2], equal_nan=True)
        self.assertEqual(radar.fields[0].dBZ.shape, (2, 500))
        np.testing.assert_allclose(radar.fields[0].azimuth, [0, 180])

    def test_sab_minimal_public_reader(self):
        radar = self._check_read('SAB', sab_fixture())
        np.testing.assert_allclose(radar.fields[0].dBZ[0], [0, 1, np.nan, 2], equal_nan=True)

    def test_pa_minimal_public_reader(self):
        radar = self._check_read('PA', pa_fixture())
        np.testing.assert_allclose(radar.fields[0].dBZ[0], [0, 1, np.nan, 2], equal_nan=True)

    def test_wsr98d_minimal_public_reader(self):
        radar = self._check_read('WSR98D', wsr98d_fixture())
        np.testing.assert_allclose(radar.fields[0].dBZ[0], [0, 1, np.nan, 2], equal_nan=True)

    def test_cc_different_gate_counts_per_sweep(self):
        from pycwr.io import read_CC
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'mixed.bin'
            path.write_bytes(cc_fixture(gate_counts=(4, 2)))
            radar = read_CC(path)
        self.assertEqual([f.dBZ.shape for f in radar.fields], [(2, 4), (2, 2)])

    def test_reader_bridge_pyart_aliases(self):
        from importlib import import_module
        for family, payload in [('CC', cc_fixture()), ('SC', sc_fixture()),
                                ('SAB', sab_fixture()), ('PA', pa_fixture()),
                                ('WSR98D', wsr98d_fixture())]:
            with self.subTest(family=family):
                module = import_module('pycwr.io.' + family + 'File')
                raw = getattr(module, family + 'BaseData')(io.BytesIO(payload))
                bridge = getattr(module, family + '2NRadar')(raw)
                for exporter in (bridge.ToPyartRadar, bridge.to_pyart_radar):
                    radar = exporter(use_external=False)
                    self.assertEqual(radar.nrays, 2)
                    self.assertIn('reflectivity', radar.fields)

    def test_raw_and_bridge_metadata_getters_preserve_values_and_units(self):
        from importlib import import_module
        import datetime
        # SAB's legacy raw/bridge range getters use km; the other protocols use m.
        cases = [('CC', cc_fixture(), 15, 150000, 10),
                 ('SC', sc_fixture(), 16, 150000, 10),
                 ('SAB', sab_fixture(), 15, 150, 1),
                 ('PA', pa_fixture(), 15, 150000, 1),
                 ('WSR98D', wsr98d_fixture(), 15, 150000, 1)]
        for family, payload, speed, max_range, duration in cases:
            with self.subTest(family=family):
                module = import_module('pycwr.io.' + family + 'File')
                raw = getattr(module, family + 'BaseData')(io.BytesIO(payload))
                bridge = getattr(module, family + '2NRadar')(raw)
                for reader in (raw, bridge):
                    np.testing.assert_allclose(reader.get_nyquist_velocity(), [speed, speed])
                    np.testing.assert_allclose(reader.get_unambiguous_range(), [max_range, max_range])
                    np.testing.assert_allclose(reader.get_azimuth(), [0, 180])
                    times = reader.get_scan_time()
                    self.assertEqual(times.shape, (2,))
                    self.assertEqual(times[0], datetime.datetime(2026, 1, 1))
                    self.assertEqual(times[1] - times[0], datetime.timedelta(seconds=duration))
                    # PA's bridge exposes these indices as attributes, not getter methods.
                    if reader is raw or family != 'PA':
                        np.testing.assert_array_equal(reader.get_sweep_start_ray_index(), [0])
                        np.testing.assert_array_equal(reader.get_sweep_end_ray_index(), [1])
                        np.testing.assert_array_equal(reader.get_rays_per_sweep(), [2])
                if family == 'PA':
                    np.testing.assert_allclose(raw.get_elevation(), [1, 1])
                    np.testing.assert_allclose(bridge.get_dbz_range_per_radial(4), [250, 500, 750, 1000])
                else:
                    np.testing.assert_allclose(bridge.get_NRadar_nyquist_speed(), [speed])
                    np.testing.assert_allclose(bridge.get_NRadar_unambiguous_range(), [max_range])

    def test_cc_metadata_getters_follow_unequal_sweep_ray_counts(self):
        from pycwr.io.CCFile import CCBaseData, CC2NRadar
        raw = CCBaseData(io.BytesIO(cc_fixture(gate_counts=(4, 2), ray_counts=(2, 3))))
        for reader in (raw, CC2NRadar(raw)):
            np.testing.assert_array_equal(reader.get_sweep_start_ray_index(), [0, 2])
            np.testing.assert_array_equal(reader.get_sweep_end_ray_index(), [1, 4])
            np.testing.assert_array_equal(reader.get_rays_per_sweep(), [2, 3])
            np.testing.assert_allclose(reader.get_nyquist_velocity(), [15] * 5)
            np.testing.assert_allclose(reader.get_unambiguous_range(), [150000] * 5)
            np.testing.assert_allclose(reader.get_elevation(), [1, 1, 2, 2, 2])

    def test_split_sweep_interpolation_and_legacy_aliases(self):
        from importlib import import_module
        for family, payload in [('SAB', sab_fixture(split=True)), ('WSR98D', wsr98d_fixture(split=True))]:
            with self.subTest(family=family):
                module = import_module('pycwr.io.' + family + 'File')
                raw = getattr(module, family + 'BaseData')(io.BytesIO(payload))
                bridge = getattr(module, family + '2NRadar')(raw)
                self.assertEqual(bridge.get_reomve_radial_num(), [0, 1])
                if family == 'WSR98D':
                    np.testing.assert_array_equal(bridge.get_dbz_idx(), [0])
                    np.testing.assert_array_equal(bridge.get_v_idx(), [1])
                np.testing.assert_allclose(bridge.azimuth, [180, 0])
                np.testing.assert_allclose(bridge.fields['dBZ'],
                    [[10, 11, np.nan, 12], [0, 1, np.nan, 2]], equal_nan=True)
                np.testing.assert_array_equal(bridge.get_sweep_start_ray_index(), [0])
                np.testing.assert_array_equal(bridge.get_sweep_end_ray_index(), [1])
                # Explicit compatibility interpolation must remap the updated source by azimuth.
                raw.radial[0]['fields']['dBZ'] = np.array([30, 31, np.nan, 32], dtype=np.float32)
                if family == 'WSR98D':
                    bridge.interp_VCP26([0], [1])
                else:
                    bridge.interp_dBZ(0, 1)
                np.testing.assert_allclose(raw.radial[3]['fields']['dBZ'],
                    [30, 31, np.nan, 32], equal_nan=True)
                np.testing.assert_allclose(raw.radial[2]['fields']['dBZ'],
                    [10, 11, np.nan, 12], equal_nan=True)

    def test_cc_raw_decoder_uses_only_declared_gates(self):
        from pycwr.io.CCFile import CCBaseData
        raw = CCBaseData(io.BytesIO(cc_fixture()))
        for radial in raw.radial:
            self.assertEqual(radial["fields"]["dBZ"].size, 4)

    def test_readers_reject_truncated_payloads(self):
        import pycwr.io as radar_io
        for family, payload in [('CC', cc_fixture()), ('SC', sc_fixture()),
                                ('SAB', sab_fixture()), ('PA', pa_fixture())]:
            with self.subTest(family=family), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'truncated.bin'
                path.write_bytes(payload[:-1])
                with self.assertRaises(ValueError):
                    getattr(radar_io, 'read_' + family)(path)

    def test_compressed_pa_public_reader(self):
        from pycwr.io import read_auto
        payload = pa_fixture()
        with tempfile.TemporaryDirectory() as directory:
            for suffix, encoded in [('bz2', bz2.compress(payload)), ('gz', gzip.compress(payload))]:
                path = Path(directory) / ('sample.' + suffix)
                path.write_bytes(encoded)
                self.assertEqual(read_auto(path).nrays, 2)
            path = Path(directory) / 'sample.zip'
            with zipfile.ZipFile(path, 'w') as archive:
                archive.writestr('sample.bin', payload)
            self.assertEqual(read_auto(path).nrays, 2)

    def test_format_detection_keeps_stream_open_and_restores_cursor(self):
        from pycwr.io.util import radar_format
        for family, payload in [('CC', cc_fixture()), ('SC', sc_fixture()),
                                ('SAB', sab_fixture()), ('PA', pa_fixture())]:
            with self.subTest(family=family):
                stream = io.BytesIO(payload)
                stream.seek(19)
                self.assertEqual(radar_format(stream), family)
                self.assertFalse(stream.closed)
                self.assertEqual(stream.tell(), 19)

    def test_public_readers_accept_file_like_inputs(self):
        from pycwr.io import read_auto
        for family, payload in [('CC', cc_fixture()), ('SC', sc_fixture()),
                                ('SAB', sab_fixture()), ('PA', pa_fixture())]:
            with self.subTest(family=family):
                radar = read_auto(io.BytesIO(payload), station_lon=120, station_lat=30, station_alt=100)
                self.assertEqual(radar.nrays, 2)


class InteropPublicContracts(unittest.TestCase):
    def _radar(self):
        from pycwr.io import read_PA
        return read_PA(io.BytesIO(pa_fixture()))

    def test_internal_pyart_does_not_import_optional_pyart(self):
        from pycwr.core.interop import resolve_pyart_radar_class
        from pycwr.core.PyartRadar import Radar
        with mock.patch('pycwr.core.interop._import_external_pyart_radar_class', side_effect=RuntimeError('broken optional dependency')):
            self.assertEqual(resolve_pyart_radar_class(use_external=False), (Radar, False))

    def test_missing_external_pyart_strict_and_fallback(self):
        radar = self._radar()
        with mock.patch('pycwr.core.interop._import_external_pyart_radar_class', return_value=None):
            with self.assertRaises(ImportError):
                radar.to_pyart_radar(use_external=True, strict=True)
            self.assertEqual(radar.to_pyart_radar(use_external=True, strict=False).nrays, 2)

    def test_pyart_exports_location_fields_and_independent_arrays(self):
        prd = self._radar()
        radar = prd.to_pyart_radar(use_external=False)
        np.testing.assert_allclose(radar.longitude['data'], [120])
        self.assertIn('reflectivity', radar.fields)
        radar.fields['reflectivity']['data'][0, 0] = 99
        self.assertEqual(float(prd.fields[0].dBZ[0, 0]), 0.0)

    def test_xradar_optional_datatree_contract(self):
        prd = self._radar()
        with mock.patch('pycwr.core.interop._import_datatree_class', return_value=None):
            with self.assertRaises(ImportError):
                prd.to_xradar(strict=True, force_rebuild=True)
            sweeps = prd.to_xradar(strict=False, force_rebuild=True)
        self.assertEqual(list(sweeps), ['sweep_0'])
        self.assertEqual(float(sweeps['sweep_0'].longitude), 120)

    def test_internal_radar_sweep_iteration_and_geometry(self):
        radar = self._radar().to_pyart_radar(use_external=False)
        self.assertEqual(list(radar.iter_start()), [0])
        self.assertEqual(list(radar.iter_end()), [1])
        self.assertEqual(list(radar.iter_start_end()), [(0, 1)])
        self.assertEqual(list(radar.iter_slice()), [slice(0, 2)])
        np.testing.assert_array_equal(next(radar.iter_azimuth()), radar.get_azimuth(0))
        np.testing.assert_array_equal(next(radar.iter_elevation()), radar.get_elevation(0))
        np.testing.assert_array_equal(next(radar.iter_field('reflectivity')), radar.get_field(0, 'reflectivity'))
        self.assertEqual(radar.get_nyquist_vel(0), 15.0)
        for coordinates in (radar.get_gate_x_y_z(0), radar.get_gate_lat_lon_alt(0, reset_gate_coords=True)):
            for coordinate in coordinates:
                self.assertEqual(coordinate.shape, (2, 4))
                self.assertTrue(np.isfinite(coordinate).all())
        with self.assertRaises(IndexError):
            radar.get_start(1)
        with self.assertRaises(KeyError):
            radar.get_field(0, 'absent')
        output = io.StringIO()
        radar.info(out=output)
        self.assertIn('reflectivity', output.getvalue())

    def test_internal_radar_add_fields_and_extract_are_independent(self):
        radar = self._radar().to_pyart_radar(use_external=False)
        values = radar.get_field(0, 'reflectivity', copy=True)
        radar.add_field_like('reflectivity', 'copy', values)
        self.assertEqual(radar.fields['copy']['units'], radar.fields['reflectivity']['units'])
        with self.assertRaises(ValueError):
            radar.add_field('copy', {'data': values})
        with self.assertRaises(ValueError):
            radar.add_field('bad', {'data': np.zeros((1, 1))})
        subset = radar.extract_sweeps([0])
        np.testing.assert_array_equal(subset.fields['copy']['data'], values)
        subset.range['data'][0] = 999
        subset.latitude['data'][0] = 99
        self.assertEqual(radar.range['data'][0], 250)
        self.assertEqual(radar.latitude['data'][0], 30)

    def test_extract_sweeps_preserves_optional_ray_metadata(self):
        radar = self._radar().to_pyart_radar(use_external=False)
        radar.rotation = {'data': np.array([10.0, 20.0])}
        radar.heading = {'data': np.array([30.0, 40.0])}
        radar.ray_angle_res = {'data': np.array([180.0])}
        subset = radar.extract_sweeps([0])
        np.testing.assert_array_equal(subset.rotation['data'], [10, 20])
        np.testing.assert_array_equal(subset.heading['data'], [30, 40])
        np.testing.assert_array_equal(subset.ray_angle_res['data'], [180])

    def test_extract_sweeps_rejects_empty_fractional_and_invalid_indices(self):
        radar = self._radar().to_pyart_radar(use_external=False)
        for sweeps in ([], [0.5], [np.nan], [-1], [1]):
            with self.subTest(sweeps=sweeps), self.assertRaises(ValueError):
                radar.extract_sweeps(sweeps)

    def test_export_unknown_field_and_overwrite_contracts(self):
        prd = self._radar()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'volume.bin'
            prd.to_wsr98d(path)
            original = path.read_bytes()
            with self.assertRaises(FileExistsError):
                prd.to_wsr98d(path)
            self.assertEqual(path.read_bytes(), original)
            with self.assertRaises(ValueError):
                prd.to_wsr98d(Path(directory) / 'invalid.bin', field_names=['unknown'])
            with self.assertRaises(ValueError):
                prd.to_nexrad_level2_msg1(Path(directory) / 'invalid.ar2v', field_names=['ZDR'])


if __name__ == '__main__':
    unittest.main()
