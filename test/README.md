# `test/` Guide

This directory serves two roles:

- runnable regression checks for behavior that must stay stable
- copyable example scripts that show new users how to use each feature

## Naming Convention

- `test_examples_<topic>.py`: tutorial-style examples built as runnable tests
- `test_regression_<topic>.py`: numeric or behavior regression coverage
- `plot_*.py`: manual plotting demos that are useful for visual inspection

## How To Run

Run the whole automated suite:

```bash
python3 -m unittest discover -s test -p 'test_*.py'
```

Run one feature group:

```bash
python3 -m unittest discover -s test -p 'test_examples_public_api.py'
python3 -m unittest discover -s test -p 'test_examples_qc.py'
python3 -m unittest discover -s test -p 'test_examples_hid.py'
python3 -m unittest discover -s test -p 'test_regression_geometry.py'
```

If the required sample radar files are not available in the current workspace,
sample-driven tests skip automatically. Slow real-sample 3D network tests are
opt-in via `PYCWR_RUN_SLOW_SAMPLE_TESTS=1`.

Recommended pre-release order:

```bash
python3 -m unittest discover -s test -p 'test_regression_*.py'
python3 -m unittest discover -s test -p 'test_examples_*.py'
PYCWR_RUN_SLOW_SAMPLE_TESTS=1 python3 -m unittest discover -s test -p 'test_examples_interp.py'
```

## README Feature Matrix

This matrix tracks the public features documented in `README.md` and the test
entry points that validate them.

| README feature group | Primary tests | Default profile |
| --- | --- | --- |
| Reading, inspection, `PRD` summary, aligned/native reflectivity | `test_examples_public_api.py`, `test_regression_io.py` | fast |
| Single-radar products: `CR`, `CAPPI`, `CAPPI_3D`, `VIL`, `ET` on `xy` and `lonlat` grids | `test_examples_public_api.py`, `test_regression_compatibility.py` | fast |
| Plotting quickstart: PPI, section, RHI | `test_examples_public_api.py`, `test_examples_sections.py` | fast |
| Dual-pol QC and `clear_air_mode` | `test_examples_qc.py`, `test_examples_xband.py` | fast |
| Hydrometeor classification | `test_examples_hid.py` | fast |
| Wind retrieval: `VAD`, `VVP`, `VWP`, stored `VWP` product, wind plotting | `test_examples_wind.py`, `test_examples_public_api.py` | fast |
| Export and interop: `to_wsr98d`, `to_nexrad_level2_msg31`, `to_nexrad_level2_msg1`, module-level writer APIs, Py-ART/xradar | `test_radial_export.py`, `test_examples_public_api.py` | fast |
| Web viewer: `create_app`, `launch`, path-security boundaries | `test_examples_public_api.py`, `test_examples_hid.py`, `test_regression_security.py` | fast |
| Multi-radar compositing: external config, products, QC options, plotting, netCDF output | `test_examples_interp.py` | slow for real-sample workflows |

## Feature Map

- `test_examples_public_api.py`: first-stop examples for reading, inspecting, plotting, web-viewer smoke tests, and generating `xy`/`lonlat` products from a `PRD`
- `test_examples_qc.py`: dual-pol QC primitives and end-to-end QC workflow examples
- `test_examples_hid.py`: hydrometeor classification examples with and without temperature profiles
- `test_examples_xband.py`: X-band attenuation correction and QC examples
- `test_examples_sections.py`: section and RHI extraction examples
- `test_examples_interp.py`: multi-radar interpolation, compositing, and network-product examples
- `test_regression_io.py`: parser and sample I/O regression checks
- `test_regression_geometry.py`: radar/geographic geometry round trips and projection validation
- `test_regression_grid.py`: gridding and low-level Cartesian helper regression checks
- `test_regression_compatibility.py`: import-path and optional-dependency compatibility checks
- `test_regression_security.py`: path-boundary and malformed-input regression checks
- `test_radial_export.py`: round-trip export checks for `PRD.to_*` and module-level writer APIs
- `plot_draw_compare.py`: manual plotting comparison script for visual review

## Recommended Learning Order

1. Start with `test_examples_public_api.py`
2. Continue with `test_examples_sections.py`, `test_examples_qc.py`, and `test_examples_hid.py`
3. Read `test_examples_interp.py` for multi-radar workflows
4. Use `test_regression_geometry.py` when you need to understand the coordinate model
