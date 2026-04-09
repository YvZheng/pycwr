# pycwr

`pycwr` is a Python toolkit for operational Chinese weather radar workflows.
It covers radar base-data reading, geometry, plotting, quality control,
hydrometeor classification, single-radar wind retrieval, multi-radar
compositing, and export.

- Current release: `1.0.8`
- [中文说明](README_CN.md)
- [API reference](docs/api_reference.md)
- [Radar network quickstart](docs/radar_network_quickstart_en.md)
- [Read the Docs](https://pycwr.readthedocs.io/en/latest/)
- [Test guide](test/README.md)
- [Draw quickstart](docs/draw_quickstart.md)

## Why 1.0.8 matters

`1.0.8` continues the first stable release line intended to be usable for production-style
usage and GitHub distribution.
This refresh keeps the PA reader behavior aligned while improving the
maintainability of the PA decode path and preserving the plotting colormap
compatibility update in the public release line.

Highlights:

- compatibility-focused readers for common Chinese weather radar formats
- clearer `PRD` data model for sweep inspection and downstream processing
- aligned and native reflectivity access kept side by side when low-level
  reflectivity range is longer than Doppler range
- lighter default dependencies with plotting, QC, web viewer, and interop
  isolated into the optional full stack
- stronger geometry and regression coverage
- built-in single-radar wind retrieval workflows: `VAD`, `VVP`, `VWP`, and gridded horizontal wind volumes
- improved multi-radar compositing and reference-style CR plotting
- `read_auto` and `read_PA` now correctly recognize and build supported PA files
- `plot_ppi`, `plot_ppi_map`, and the legacy `Graph` / `GraphMap` plotting paths
  now auto-register `pycwr` colormap names such as `CN_ref` and `CN_vel`
- the PA reader implementation is cleaner internally, with explicit helpers for
  moment decoding, field padding, and range construction

## Installation

Install from PyPI:

```bash
python -m pip install pycwr
```

Install the optional full stack from PyPI:

```bash
python -m pip install "pycwr[full]"
```

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

- `pycwr 1.0.8` requires Python `>=3.9`
- `python -m pip install pycwr` is the recommended path for normal users
- base install is enough for readers, `PRD`, geometry, interpolation, and
  NetCDF-style export
- full install is recommended for plotting, map plotting, QC, Py-ART/xradar
  interop, and the web viewer
- upstream `arm_pyart` and `xradar` currently require Python `>=3.10`, so on
  Python `3.9` the full install still covers plotting, QC, and the web
  viewer, but not those two optional interop stacks
- `pandas` is pinned to `<3` in `1.0.8` for release stability
- source install is still useful when you are developing locally or rebuilding
  the Cython extension

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

Plot a mapped PPI:

```python
from pycwr.draw import plot_ppi_map

plot_ppi_map(radar, field="dBZ", sweep=0, show=True)
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
| `pycwr.retrieve` | Hydrometeor and wind retrieval | `classify_hydrometeors`, `retrieve_vad`, `retrieve_vvp`, `retrieve_vwp`, `retrieve_wind_volume_xy`, `retrieve_wind_volume_lonlat` |
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

#### Default colormap strategy

`pycwr.draw` chooses a default colormap from the plotted field when you leave
`cmap=None`.

- `dBZ` uses the reflectivity palette such as `CN_ref`
- `V` uses the velocity palette such as `CN_vel`
- `HCL` uses a discrete hydrometeor classification palette
- unknown fields fall back to `viridis`

This default strategy is shared by `plot_ppi`, `plot_ppi_map`, `plot_rhi`,
`plot_section`, and the legacy `Graph` / `GraphMap` classes.

Recommended usage:

```python
from pycwr.draw import plot_ppi, plot_ppi_map

plot_ppi(radar, field="dBZ", sweep=0, show=True)
plot_ppi_map(radar, field="dBZ", sweep=0, show=True)
```

Explicit colormap override:

```python
plot_ppi(radar, field="dBZ", sweep=0, cmap="CN_ref", cmap_bins=16, show=True)
plot_ppi_map(radar, field="V", sweep=0, cmap="CN_vel", min_max=(-27.0, 27.0), show=True)
```

Legacy class-style usage is also supported:

```python
from pycwr.draw.RadarPlot import GraphMap
from cartopy import crs as ccrs

display = GraphMap(radar, ccrs.PlateCarree())
display.plot_ppi_map(ax, 0, "dBZ", cmap="CN_ref")
```

Notes:

- `CN_ref`, `CN_vel`, and the other `pycwr` colormap names are registered
  automatically when plotting starts. You no longer need to "warm up" the
  plotting module first.
- If a field is missing in the radar file, plotting raises a field-not-found
  error. For example, some files may contain `dBZ` but not `V`.
- If you need the raw colormap registry for custom matplotlib work, you can
  explicitly import `pycwr.draw.colormap`.

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

`pycwr` ships four single-radar wind workflows:

- `retrieve_vad`: ring-wise harmonic fit on one or more sweeps
- `retrieve_vvp`: local least-squares horizontal wind retrieval on one sweep
- `retrieve_vwp`: vertical wind profile built from multiple VAD layers
- `retrieve_wind_volume_xy` / `retrieve_wind_volume_lonlat`: fixed-height gridded horizontal wind volume built from multiple VVP sweeps

Quick example:

```python
vad = radar.retrieve_vad(sweeps=[0, 1, 2], max_range_km=40.0, gate_step=4)
vvp = radar.retrieve_vvp(0, max_range_km=20.0, az_num=91, bin_num=5)
vwp = radar.retrieve_vwp(sweeps=[0, 1, 2], max_range_km=40.0, height_step=500.0)
wind = radar.retrieve_wind_volume_xy(
    XRange=np.arange(-20_000.0, 20_001.0, 10_000.0),
    YRange=np.arange(-20_000.0, 20_001.0, 10_000.0),
    level_heights=np.array([500.0, 1000.0, 1500.0]),
    sweeps=[0, 1, 2],
    max_range_km=30.0,
)
```

The stored profile product path is:

```python
radar.add_product_VWP(sweeps=[0, 1, 2], max_range_km=40.0, height_step=500.0)
```

Store the gridded horizontal wind volume in `radar.product`:

```python
radar.add_product_WIND_VOLUME_xy(
    XRange=np.arange(-20_000.0, 20_001.0, 10_000.0),
    YRange=np.arange(-20_000.0, 20_001.0, 10_000.0),
    level_heights=np.array([500.0, 1000.0, 1500.0]),
    sweeps=[0, 1, 2],
)
```

#### 3-D wind volume: what it is

`retrieve_wind_volume_xy` and `retrieve_wind_volume_lonlat` return a single-radar
fixed-height horizontal wind volume.

Treat it as:

- gridded `u/v` wind on multiple height levels
- one wind profile for each horizontal grid point
- a robust diagnostic product built from Doppler radial velocity

Do not treat it as:

- a full dynamic `u/v/w` retrieval
- a dual-Doppler or multi-radar variational wind analysis
- a substitute for vertical motion retrieval in convective dynamics studies

#### 3-D wind volume: retrieval chain

The operational logic is:

1. For each requested sweep, `pycwr` selects the velocity field, preferring
   corrected velocity when available.
2. A local VVP is solved on each sweep with a moving azimuth-range window.
3. Each sweep is summarized and, when `sweeps=None` or `sweeps="auto"`, the
   strongest sweeps are selected automatically.
4. The valid sweep-level VVP samples are interpolated horizontally onto the
   target grid with inverse-distance weighting.
5. Fixed-height wind values are reconstructed from the selected sweeps with
   local vertical aggregation or short-gap linear interpolation.
6. Diagnostics such as fit residual, support counts, and heuristic quality
   score are stored with the result.

This makes the product more stable on real radar volumes with missing velocity,
weak echo regions, or sweep-to-sweep range differences.

#### 3-D wind volume: recommended API

Use:

- `radar.retrieve_wind_volume_xy(...)` for Cartesian `x/y/z`
- `radar.retrieve_wind_volume_lonlat(...)` for regular lon/lat grids
- `radar.add_product_WIND_VOLUME_xy(...)` or
  `radar.add_product_WIND_VOLUME_lonlat(...)` to store the result in
  `radar.product`

Recommended defaults for production-style usage:

- use `sweeps=None` to enable auto sweep selection
- use `sweeps="all"` only when you explicitly want every available sweep
- set a realistic `max_range_km` instead of using the whole radar volume
- increase `workers` on multi-core machines when the target grid is not tiny

Example with explicit operational settings:

```python
wind = radar.retrieve_wind_volume_xy(
    XRange=np.arange(-60_000.0, 60_001.0, 5_000.0),
    YRange=np.arange(-60_000.0, 60_001.0, 5_000.0),
    level_heights=np.arange(500.0, 3_500.0, 500.0),
    sweeps=None,
    max_range_km=60.0,
    az_num=31,
    bin_num=5,
    azimuth_step=12,
    range_step=6,
    horizontal_radius_m=8_000.0,
    max_horizontal_radius_m=14_000.0,
    horizontal_min_neighbors=3,
    vertical_tolerance_m=400.0,
    max_vertical_gap_m=1_200.0,
    workers=4,
)
```

#### 3-D wind volume: key arguments

- `XRange`, `YRange`, `level_heights`: target output grid in meters
- `XLon`, `YLat`: target output grid in degrees for lon/lat mode
- `sweeps`: `None` or `"auto"` enables auto selection; `"all"` keeps all
  sweeps; a list such as `[0, 1, 2]` forces explicit sweeps
- `max_range_km`: strongest first-order control for retrieval stability
- `az_num`, `bin_num`: local VVP window size
- `azimuth_step`, `range_step`: output thinning of the local sweep VVP
- `horizontal_radius_m`: horizontal radius used when mapping sweep VVP samples
  to the target grid
- `max_horizontal_radius_m`: optional expanded search radius when local support
  is sparse
- `vertical_tolerance_m`: local height matching tolerance
- `max_vertical_gap_m`: maximum allowed vertical interpolation gap
- `workers`: process-based parallelism for sweep-level and grid-level work

Practical tuning rules:

- use a wider `az_num` when velocity is sparse or noisy
- keep `bin_num` modest; too small is noisy and too large oversmooths
- do not push `max_range_km` farther than the range where velocity remains
  reasonably continuous
- on sparse cases, it is better to slightly increase search radius than to
  force too many poor sweeps

#### 3-D wind volume: output fields and metadata

The returned dataset includes:

- `u`, `v`: eastward and northward wind components
- `wind_speed`, `wind_direction`
- `fit_rmse`: local retrieval residual
- `valid_count`: effective sample count
- `source_sweep_count`: number of sweeps contributing to the voxel
- `neighbor_count`: horizontal support count
- `effective_vertical_support`: number of vertically supporting sweep samples
- `sweep_spread_m`: vertical spread of supporting sweeps
- `quality_score`: heuristic score from `0` to `100`
- `quality_flag`: `0/1/2/3` for invalid, low-confidence, usable, and
  high-confidence voxels

Important attrs include:

- `selection_mode`
- `requested_sweeps`
- `selected_sweeps`
- `rejected_sweeps`
- `rejected_sweep_reasons`
- `workers_used`
- `parallel_mode`

These attrs are often more informative than the plot itself when diagnosing
why a real case contains gaps.

#### 3-D wind volume: quality interpretation

If a real case has large missing areas, that does not automatically mean the
algorithm failed.

Common causes are:

- no precipitation or very weak echo, so there is no valid Doppler velocity
- poor azimuth coverage in the local window
- sweep geometry that does not support the requested height well
- too aggressive a range limit or too small a horizontal search radius

When checking a case, inspect:

- `quality_score` and `quality_flag`
- `fit_rmse`
- `valid_count`
- `source_sweep_count`
- `selected_sweeps` and `rejected_sweep_reasons`

#### 3-D wind volume: performance and parallelism

The `workers` argument uses process-based parallelism, not Python threads.

In practice:

- small grids may see little benefit
- moderate and large grids usually benefit from `workers=2` or `workers=4`
- very large worker counts may stop scaling once process overhead dominates

For reproducible validation, compare a serial run and a parallel run on the
same case and confirm that `u`, `v`, and `quality_flag` match.

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
- [Radar network quickstart](docs/radar_network_quickstart_en.md): beginner-friendly multi-radar workflow
- [Test guide](test/README.md): runnable examples by feature
- [Draw quickstart](docs/draw_quickstart.md): plotting-focused entry points
- [Read the Docs](https://pycwr.readthedocs.io/en/latest/): hosted documentation site
- Web viewer: launch it locally with `python scripts/LaunchGUI.py`

## Release notes for users

For `1.0.8`, the most important user-visible behavior rules are:

- radar reading returns one stable `PRD` object across supported formats
- low-level reflectivity can now be queried explicitly in aligned or native mode
- QC and hydrometeor classification can write corrected fields back into `PRD`
- single-radar wind retrieval is now part of the public API, including gridded horizontal wind volumes
- multi-radar compositing has a documented high-level workflow
- packaging is split into base and full dependency sets for lower installation
  friction
- PA samples matching the supported format markers are routed through `read_PA`
  instead of being misidentified as `WSR98D`
- default plotting now resolves `pycwr` colormap names consistently across both
  the easy plotting API and legacy class-style plotting
