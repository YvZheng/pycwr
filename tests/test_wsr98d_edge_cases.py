import numpy as np
import pytest
import sys, os, pathlib
sys.path.insert(0, os.fspath(pathlib.Path(__file__).resolve().parents[1]))

import io
import types
import importlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def import_wsr98d():
    # Inject lightweight package shells to avoid running pycwr/__init__
    pkg = types.ModuleType('pycwr')
    pkg.__path__ = [str(ROOT / 'pycwr')]
    sys.modules.setdefault('pycwr', pkg)
    # Now import via package path so relative imports resolve
    wsrmod = importlib.import_module('pycwr.io.WSR98DFile')
    try:
        errmod = importlib.import_module('pycwr.io.errors')
        err_cls = getattr(errmod, 'CINRADFormatError')
    except Exception:
        err_cls = AssertionError
    return wsrmod.WSR98D2NRadar, wsrmod.WSR98DBaseData, err_cls


WSR98D2NRadar, WSR98DBaseData, CINRADFormatError = import_wsr98d()


class FakeWSR98D:
    def __init__(self, radial, sweep_starts, sweep_ends, *,
                 log_res=250.0, dop_res=250.0, nyq=15.0, max_range=150000.0,
                 elevs=None, azis=None, sitename=b"TEST", task_name=b"TASKXX"):
        self.radial = radial
        self.sweep_start_ray_index = np.array(sweep_starts, dtype=int)
        self.sweep_end_ray_index = np.array(sweep_ends, dtype=int)
        self.header = {
            "CutConfig": {
                "LogResolution": np.array([log_res], dtype=float),
                "DopplerResolution": np.array([dop_res], dtype=float),
                "NyquistSpeed": np.array([nyq], dtype=float),
                "MaximumRange": np.array([max_range], dtype=float),
                "Elevation": np.array([0.5], dtype=float),
                "Azimuth": np.array([0.0], dtype=float),
            },
            "SiteConfig": {
                "Latitude": 30.0,
                "Longitude": 120.0,
                "Height": 50.0,
                "Frequency": 9400.0,
                "BeamWidthHori": 999.0,
                "BeamWidthVert": -1.0,
                "SiteName": sitename,
            },
            "TaskConfig": {
                "ScanType": 0,  # ppi
                "TaskName": task_name,
                "PulseWidth": 1000,
            },
        }

    def get_scan_type(self):
        return "ppi"

    def get_latitude_longitude_altitude_frequency(self):
        sc = self.header["SiteConfig"]
        return sc["Latitude"], sc["Longitude"], sc["Height"], sc["Frequency"]/1000.0

    def get_azimuth(self):
        return np.array([r.get("Azimuth", 0.0) for r in self.radial], dtype=float)

    def get_elevation(self):
        return np.array([r.get("Elevation", 0.5) for r in self.radial], dtype=float)

    def get_sitename(self):
        return self.header["SiteConfig"]["SiteName"].decode("utf-8", "ignore")


def _mk_ray(state, ngates=4, fields=("V", "dBZ")):
    # Minimal ray with given state and fields containing ngates floats
    arr = np.arange(ngates, dtype=np.float32)
    f = {k: arr.copy() for k in fields}
    return {"RadialState": state, "fields": f, "Azimuth": 10.0, "Elevation": 0.5}


def test_happy_path_minimal_volume():
    rays = [_mk_ray(0), _mk_ray(1), _mk_ray(2)]  # start, middle, end
    fake = FakeWSR98D(rays, sweep_starts=[0], sweep_ends=[2])
    obj = WSR98D2NRadar(fake)
    # golden samples
    assert obj.nsweeps == 1
    assert obj.nrays == 3
    assert obj.scan_type == "ppi"
    assert set(obj.fields.keys()) >= {"V", "dBZ"}
    # range length equals max bins (doppler res based)
    assert obj.range.size == obj.max_bins
    # azimuth/elevation shape matches rays
    assert obj.azimuth.shape[0] == obj.nrays
    assert obj.elevation.shape[0] == obj.nrays


def test_missing_sweep_end_raises_on_rays_per_sweep():
    # One start without a matching end
    rays = [_mk_ray(0), _mk_ray(1), _mk_ray(1)]
    fake = FakeWSR98D(rays, sweep_starts=[0], sweep_ends=[])
    obj = WSR98D2NRadar(fake)
    with pytest.raises(ValueError):
        # different-length arrays should error during subtraction
        _ = obj.get_rays_per_sweep()


def test_missing_rays_negative_count():
    # End index precedes start: negative rays count
    rays = [_mk_ray(2), _mk_ray(0)]
    fake = FakeWSR98D(rays, sweep_starts=[1], sweep_ends=[0])
    obj = WSR98D2NRadar(fake)
    counts = obj.get_rays_per_sweep()
    assert (counts <= 0).any()


def test_odd_beam_widths_propagate_to_pyart():
    pytest.importorskip("netCDF4")
    rays = [_mk_ray(0), _mk_ray(1), _mk_ray(2)]
    fake = FakeWSR98D(rays, sweep_starts=[0], sweep_ends=[2])
    obj = WSR98D2NRadar(fake)
    radar = obj.ToPyartRadar()
    bw_h = radar.instrument_parameters["radar_beam_width_h"]["data"][0]
    bw_v = radar.instrument_parameters["radar_beam_width_v"]["data"][0]
    assert bw_h == pytest.approx(999.0)
    assert bw_v == pytest.approx(-1.0)


def test_inconsistent_v_dbz_mapping_asserts():
    # Create two sweeps; first has only V, third has only dBZ
    # So indices are not adjacent -> assertion in constructor
    rays = []
    # Sweep 0 start
    rays.append(_mk_ray(0, fields=("V",)))
    rays.append(_mk_ray(1, fields=("V",)))
    rays.append(_mk_ray(2, fields=("V",)))
    # Sweep 1 start (gap) — only dBZ at this sweep start
    rays.append(_mk_ray(0, fields=("dBZ",)))
    rays.append(_mk_ray(2, fields=("dBZ",)))
    fake = FakeWSR98D(rays, sweep_starts=[0, 3], sweep_ends=[2, 4], task_name=b"TASKXX")
    with pytest.raises((CINRADFormatError, AssertionError)):
        WSR98D2NRadar(fake)


def test_repair_missing_sweep_indices():
    # Missing end repaired when flag set
    rays = [_mk_ray(0), _mk_ray(1), _mk_ray(1)]
    fake = FakeWSR98D(rays, sweep_starts=[0], sweep_ends=[])
    obj = WSR98D2NRadar(fake, repair_missing=True)
    counts = obj.get_rays_per_sweep()
    # repaired to length 1, end clamped to start
    assert counts.shape == (1,)
    assert counts[0] >= 1


def test_bad_magic_raises_cinrad_format_error():
    bio = io.BytesIO(b"NOTRSTM" + b"\x00" * 100)
    with pytest.raises(CINRADFormatError):
        WSR98DBaseData(bio)
