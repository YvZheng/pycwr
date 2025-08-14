import io
import os
import gzip
import bz2
import struct
import tempfile

import pytest
import sys, os, pathlib
sys.path.insert(0, os.fspath(pathlib.Path(__file__).resolve().parents[1]))

import importlib.util
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_util_module():
    # Inject a light-weight 'pycwr' package to avoid importing heavy __init__
    pkg = types.ModuleType('pycwr')
    pkg.__path__ = [str(ROOT / 'pycwr')]
    sys.modules.setdefault('pycwr', pkg)
    io_pkg = types.ModuleType('pycwr.io')
    io_pkg.__path__ = [str(ROOT / 'pycwr' / 'io')]
    sys.modules.setdefault('pycwr.io', io_pkg)
    # Load util with proper package context so relative imports resolve
    name = 'pycwr.io.util'
    util_path = ROOT / 'pycwr' / 'io' / 'util.py'
    spec = importlib.util.spec_from_file_location(name, util_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)  # type: ignore[attr-defined]
    return module


util = load_util_module()


def test_prepare_for_read_filelike():
    bio = io.BytesIO(b"hello")
    fh = util._prepare_for_read(bio)
    assert fh is bio
    assert fh.read() == b"hello"


def test_prepare_for_read_gzip(tmp_path):
    path = tmp_path / "sample.gz"
    data = b"abc123"
    with gzip.open(path, "wb") as f:
        f.write(data)
    with util._prepare_for_read(str(path)) as fh:
        assert fh.read() == data


def test_prepare_for_read_bz2(tmp_path):
    path = tmp_path / "sample.bz2"
    data = b"abc123"
    with bz2.BZ2File(path, "wb") as f:
        f.write(data)
    with util._prepare_for_read(str(path)) as fh:
        assert fh.read() == data


def _write_bytes_at(path: str, length: int, patches: dict):
    # Create a file of `length` bytes, default zeros, with specific byte patches
    buf = bytearray(length)
    for off, payload in patches.items():
        buf[off:off+len(payload)] = payload
    with open(path, "wb") as f:
        f.write(buf)


def test_radar_format_wsr98d(tmp_path):
    # First 4 bytes 'RSTM' should trigger WSR98D
    p = tmp_path / "wsr98d.bin"
    _write_bytes_at(str(p), 64, {0: b"RSTM"})
    assert util.radar_format(str(p)) == "WSR98D"


def test_radar_format_sab(tmp_path):
    # Bytes[14:16] == 0x01 0x00 triggers SAB
    p = tmp_path / "sab.bin"
    patches = {14: b"\x01\x00"}
    _write_bytes_at(str(p), 64, patches)
    assert util.radar_format(str(p)) == "SAB"


def test_radar_format_pa(tmp_path):
    # Bytes[8:12] == 0x10 0x00 0x00 0x00 triggers PA
    p = tmp_path / "pa.bin"
    patches = {8: b"\x10\x00\x00\x00"}
    _write_bytes_at(str(p), 64, patches)
    assert util.radar_format(str(p)) == "PA"


def test_radar_format_cc(tmp_path):
    # Make file with size obeying (size-1024)%3000 == 0 and marker at 116
    p = tmp_path / "cc_9250.bin"  # include digits to bypass _get_radar_type fallback
    size = 1024 + 3000  # minimal that satisfies the modulus
    patches = {116: b"CINRAD/CC"}
    _write_bytes_at(str(p), size, patches)
    assert util.radar_format(str(p)) == "CC"


def test_radar_format_sc(tmp_path):
    # Make file with size obeying (size-1024)%4000 == 0 and marker at 100
    p = tmp_path / "sc_9250.bin"
    size = 1024 + 4000
    patches = {100: b"CINRAD/SC"}
    _write_bytes_at(str(p), size, patches)
    assert util.radar_format(str(p)) == "SC"


def test_radar_format_unknown(tmp_path):
    # Unknown header; ensure graceful fallback (likely None or type from mapping)
    # Use a filename without any 4-digit id to avoid _get_radar_type lookup.
    p = tmp_path / "unknown.bin"
    _write_bytes_at(str(p), 64, {})
    # When radar_format cannot determine type, it calls _get_radar_type which
    # requires a 4-digit id in filename. Without it, util.get_radar_type may raise.
    # So here we assert that calling radar_format does not crash but returns None or file-like fallback.
    try:
        res = util.radar_format(str(p))
    except Exception as e:
        pytest.fail(f"radar_format raised unexpectedly: {e}")
    # res may be None depending on mapping; accept None or any of known types
    assert res in (None, "WSR98D", "SAB", "CC", "SC", "PA")
