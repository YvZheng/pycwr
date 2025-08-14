import io
import pytest
import sys, os, pathlib
sys.path.insert(0, os.fspath(pathlib.Path(__file__).resolve().parents[1]))

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IO_INIT_PATH = ROOT / 'pycwr' / 'io' / '__init__.py'
spec = importlib.util.spec_from_file_location('pycwr_io_init', IO_INIT_PATH)
ioinit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ioinit)  # type: ignore
read_auto = ioinit.read_auto


def test_read_auto_unsupported_type_raises():
    # Passing a file-like leads radar_format to return the object itself,
    # which is not a supported type string in read_auto -> TypeError.
    bio = io.BytesIO(b"not a known header")
    with pytest.raises(TypeError):
        read_auto(bio)
