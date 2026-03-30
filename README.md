# pycwr

`pycwr` is a Python toolkit for operational Chinese weather radar workflows.
It covers radar base-data reading, geometry, plotting, quality control,
hydrometeor classification, single-radar wind retrieval, multi-radar
compositing, and export.

- Current release: `1.0.5`
- [中文说明](README_CN.md)
- [API reference](docs/api_reference.md)
- [Test guide](test/README.md)
- [Draw quickstart](docs/draw_quickstart.md)
- [Web viewer quickstart](docs/web_viewer_quickstart.md)
- [Contributors](CONTRIBUTORS.txt)

## Why 1.0.5 matters

`1.0.5` continues the first stable release line intended to be usable for production-style
usage and GitHub distribution.
This refresh also includes the `Wc` variable update in the public release line.

Highlights:

- compatibility-focused readers for common Chinese weather radar formats
- clearer `PRD` data model for sweep inspection and downstream processing
- aligned and native reflectivity access kept side by side when low-level
  reflectivity range is longer than Doppler range
- lighter default dependencies with plotting, QC, web viewer, and interop
  isolated into the optional full stack
- stronger geometry and regression coverage
- built-in single-radar wind retrieval workflows: `VAD`, `VVP`, and `VWP`
- improved multi-radar compositing and reference-style CR plotting

## Installation

Base install:

```bash
python -m pip install -r requirements-core.txt
python -m pip install .
```

Full install:

```bash
python -m pip install -r requirements-full.txt
python -m pip install ".[full]"
```

Notes:

- `pycwr 1.0.5` requires Python `>=3.9`
- base install is enough for readers, `PRD`, geometry, interpolation, and
  NetCDF-style export
- full install is recommended for plotting, map plotting, QC, Py-ART/xradar
  interop, and the web viewer
- upstream `arm_pyart` and `xradar` currently require Python `>=3.10`, so on
  Python `3.9` the full install still covers plotting, QC, and the web
  viewer, but not those two optional interop stacks
- `pandas` is pinned to `<3` in `1.0.5` for release stability

Rebuild the extension after editing `pycwr/core/RadarGridC.pyx`:

```bash
python setup.py build_ext --inplace
```

Build release artifacts:

```bash
python -m build
```

## 5-minute quickstart

Read one radar file and inspect the returned volume:

```python
from pycwr.io import read_auto

radar = read_auto("Z_RADR_I_Z9046_20260317065928_O_DOR_SAD_CAP_FMT.bin.bz2")
print(radar.summary())
print(radar.available_fields())
print(radar.sweep_summary()[0])
```

Get one field from one sweep:

```python
dBZ0 = radar.get_sweep_field(0, "dBZ")
velocity0 = radar.get_sweep_field(0, "V")
```

Plot a PPI:

```python
from pycwr.draw import plot_ppi

plot_ppi(radar, field="dBZ", sweep=0, show=True)
```

Extract a vertical section:

```python
from pycwr.draw import plot_section

plot_section(radar, start=(-50, 0), end=(50, 0), field="dBZ", show=True)
```

Generate a simple product:

```python
import numpy as np

x = np.arange(-150_000.0, 150_001.0, 1_000.0)
y = np.arange(-150_000.0, 150_001.0, 1_000.0)
radar.add_product_CR_xy(x, y)
print(radar.product)
```

## API map

| Module | What it is for | Recommended starting points |
| --- | --- | --- |
| `pycwr.io` | Read and write radar base data | `read_auto`, `read_WSR98D`, `read_SAB`, `read_CC`, `read_SC`, `read_PA` |
| `pycwr.core` | Central volume object, geometry, export helpers | `PRD`, `radar.summary()`, `radar.get_sweep_field()` |
| `pycwr.draw` | Plotting and quick-look figures | `plot_ppi`, `plot_ppi_map`, `plot_rhi`, `plot_section`, `plot_vvp`, `plot_wind_profile` |
| `pycwr.qc` | Dual-pol QC | `apply_dualpol_qc`, `run_dualpol_qc` |
| `pycwr.retrieve` | Hydrometeor and wind retrieval | `classify_hydrometeors`, `retrieve_vad`, `retrieve_vvp`, `retrieve_vwp` |
| `pycwr.interp` | Multi-radar compositing | `run_radar_network_3d` |
| `pycwr.GraphicalInterface` | Local web viewer | `create_app`, `launch` |

## Core object model

All readers return `pycwr.core.NRadar.PRD`.

The most important parts of `PRD` are:

- `fields`: one `xarray.Dataset` per sweep
- `scan_info`: site and scan metadata
- `extended_fields`: native sidecar fields when aligned and native ranges differ
- `product`: computed product dataset

Useful inspection helpers:

- `summary()`: compact full-volume summary
- `available_fields(sweep=None, range_mode="aligned")`
- `sweep_summary()`
- `get_sweep_field(sweep, field_name, range_mode="aligned", sort_by_azimuth=False)`
- `get_native_sweep_field(sweep, field_name)`
- `ordered_az(inplace=False)`

### Aligned vs native reflectivity

For some low sweeps, reflectivity may exist in two forms:

- aligned: historical shared grid used by old processing chains
- native: original reflectivity range before Doppler-driven truncation

Use:

- `range_mode="aligned"` for historical compatibility
- `range_mode="native"` when full low-level reflectivity coverage matters

Example:

```python
aligned = radar.get_sweep_field(0, "dBZ", range_mode="aligned")
native = radar.get_sweep_field(0, "dBZ", range_mode="native")
```

## Common workflows

### Plotting

Recommended public plotting APIs:

```python
from pycwr.draw import (
    plot,
    plot_ppi,
    plot_ppi_map,
    plot_rhi,
    plot_section,
    plot_section_lonlat,
    plot_vvp,
    plot_wind_profile,
)
```

These functions return an `EasyPlotResult` with `fig`, `ax`, and `artist`.

### Products

Common `PRD` product methods:

- `add_product_CR_xy`
- `add_product_CAPPI_xy`
- `add_product_CAPPI_3d_xy`
- `add_product_VIL_xy`
- `add_product_ET_xy`
- `add_product_CR_lonlat`
- `add_product_CAPPI_lonlat`
- `add_product_VIL_lonlat`
- `add_product_ET_lonlat`
- `add_product_VWP`

Example:

```python
radar.add_product_CAPPI_xy(x, y, 3000.0)
radar.add_product_VIL_xy(x, y, [1000.0, 2000.0, 3000.0])
```

### Quality control

```python
from pycwr.qc import apply_dualpol_qc

qc_radar = apply_dualpol_qc(radar, inplace=False, clear_air_mode="mask")
```

Common corrected fields include `Zc`, `ZDRc`, `PhiDPc`, `KDPc`, and mask fields
such as `QC_MASK` and `CLEAR_AIR_MASK` when enabled.

### Hydrometeor classification

```python
hcl_radar = radar.classify_hydrometeors(
    inplace=False,
    band="C",
    profile_height=[0.0, 2000.0, 4000.0, 8000.0],
    profile_temperature=[24.0, 12.0, 2.0, -16.0],
    confidence_field="HCL_CONF",
)
```

You can also classify directly from arrays with
`pycwr.retrieve.classify_hydrometeors(...)`.

### Wind retrieval

`pycwr` ships three single-radar wind workflows:

- `retrieve_vad`: ring-wise harmonic fit on one or more sweeps
- `retrieve_vvp`: local least-squares horizontal wind retrieval on one sweep
- `retrieve_vwp`: vertical wind profile built from multiple VAD layers

Example:

```python
vad = radar.retrieve_vad(sweeps=[0, 1, 2], max_range_km=40.0, gate_step=4)
vvp = radar.retrieve_vvp(0, max_range_km=20.0, az_num=91, bin_num=5)
vwp = radar.retrieve_vwp(sweeps=[0, 1, 2], max_range_km=40.0, height_step=500.0)
```

The stored profile product path is:

```python
radar.add_product_VWP(sweeps=[0, 1, 2], max_range_km=40.0, height_step=500.0)
```

### Export and interop

Common export helpers:

- `radar.to_wsr98d(...)`
- `radar.to_nexrad_level2_msg31(...)`
- `radar.to_nexrad_level2_msg1(...)`
- `radar.to_pyart_radar(...)`
- `radar.to_xradar(...)`
- `radar.to_cfgridded_netcdf(...)`

Use these when you need to move `pycwr` data into Py-ART, xradar, or other
NetCDF-style workflows.

### Multi-radar compositing

```python
from pycwr.interp import run_radar_network_3d
```

This entry point builds 3D network products on a regular lon/lat grid and can
write NetCDF directly. It is the recommended high-level interface for network
CR and CAPPI workflows.

### Web viewer

```python
from pycwr.GraphicalInterface import create_app, launch
```

Or use the script:

```bash
python scripts/LaunchGUI.py
```

The viewer is local-only by design and requires a token for API access.

## What to read next

- [API reference](docs/api_reference.md): detailed public interface notes
- [Test guide](test/README.md): runnable examples by feature
- [Draw quickstart](docs/draw_quickstart.md): plotting-focused entry points
- [Web viewer quickstart](docs/web_viewer_quickstart.md): local viewer setup

## Release notes for users

For `1.0.5`, the most important user-visible behavior rules are:

- radar reading returns one stable `PRD` object across supported formats
- low-level reflectivity can now be queried explicitly in aligned or native mode
- QC and hydrometeor classification can write corrected fields back into `PRD`
- single-radar wind retrieval is now part of the public API
- multi-radar compositing has a documented high-level workflow
- packaging is split into base and full dependency sets for lower installation
  friction
