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

## Feature Map

- `test_examples_public_api.py`: first-stop examples for reading, inspecting, exporting, and generating products from a `PRD`
- `test_examples_qc.py`: dual-pol QC primitives and end-to-end QC workflow examples
- `test_examples_hid.py`: hydrometeor classification examples with and without temperature profiles
- `test_examples_xband.py`: X-band attenuation correction and QC examples
- `test_examples_sections.py`: section and RHI extraction examples
- `test_examples_interp.py`: multi-radar interpolation and compositing examples
- `test_regression_io.py`: parser and sample I/O regression checks
- `test_regression_geometry.py`: radar/geographic geometry round trips and projection validation
- `test_regression_grid.py`: gridding and low-level Cartesian helper regression checks
- `test_regression_compatibility.py`: import-path and optional-dependency compatibility checks
- `test_regression_security.py`: path-boundary and malformed-input regression checks
- `plot_draw_compare.py`: manual plotting comparison script for visual review

## Recommended Learning Order

1. Start with `test_examples_public_api.py`
2. Continue with `test_examples_sections.py`, `test_examples_qc.py`, and `test_examples_hid.py`
3. Read `test_examples_interp.py` for multi-radar workflows
4. Use `test_regression_geometry.py` when you need to understand the coordinate model
