# pycwr API Reference

This document is a practical reference for the public APIs most users need in
`pycwr 1.0.1`. It does not try to mirror every internal helper. The focus is:

- where to enter the package
- which object each API returns
- which arguments control behavior
- how aligned and native reflectivity workflows differ

For runnable examples grouped by feature, see [../test/README.md](../test/README.md).

## Package layout

| Module | Main purpose | Recommended entry points |
| --- | --- | --- |
| `pycwr.io` | Read and write radar base data | `read_auto`, `read_WSR98D`, `read_SAB`, `read_CC`, `read_SC`, `read_PA` |
| `pycwr.core` | Central volume object, geometry, export helpers | `PRD`, `radar.summary()`, `radar.get_sweep_field()` |
| `pycwr.draw` | Plotting and quick-look figures | `plot_ppi`, `plot_ppi_map`, `plot_rhi`, `plot_section`, `plot_vvp`, `plot_wind_profile` |
| `pycwr.qc` | Dual-pol quality control | `apply_dualpol_qc`, `run_dualpol_qc` |
| `pycwr.retrieve` | Hydrometeor and wind retrieval | `classify_hydrometeors`, `retrieve_vad`, `retrieve_vvp`, `retrieve_vwp` |
| `pycwr.interp` | Multi-radar compositing | `run_radar_network_3d`, `radar_network_3d_to_netcdf` |
| `pycwr.GraphicalInterface` | Local web viewer | `create_app`, `launch` |

## Cross-cutting conventions

### Units

- distance and height: meters
- internal trigonometric calculations: radians
- public azimuth, elevation, and fixed-angle values: degrees
- lon/lat input and output: degrees

Some historical plotting interfaces still accept kilometers because that was
the historical public behavior.

### `range_mode`

Several APIs accept:

- `range_mode="aligned"`: historical aligned workflow
- `range_mode="native"`: native long-range reflectivity when available

Recommended rule:

- use `aligned` for backward-compatible plots and products
- use `native` when low-level reflectivity coverage matters

Velocity fields usually remain aligned because their valid range is often
shorter than reflectivity.

### Optional dependencies

Base install covers:

- readers
- `PRD`
- geometry
- interpolation
- NetCDF-style export

Full install adds:

- plotting
- map plotting
- QC
- web viewer
- Py-ART and xradar interop

Install the full stack with:

```bash
pip install "pycwr[full]"
```

## `pycwr.io`

### Recommended reader: `read_auto`

```python
read_auto(
    filename,
    station_lon=None,
    station_lat=None,
    station_alt=None,
    effective_earth_radius=None,
)
```

Purpose:

- detect the radar family automatically
- parse the file into a `PRD`
- preserve the package's compatible geometry and sweep layout

Arguments:

- `filename`: radar file path
- `station_lon`, `station_lat`, `station_alt`: optional site override values
- `effective_earth_radius`: optional beam-geometry radius in meters

Returns:

- `pycwr.core.NRadar.PRD`

Typical usage:

```python
from pycwr.io import read_auto

radar = read_auto("your_radar_file.bin.bz2")
print(radar.summary())
```

### Format-specific readers

```python
read_WSR98D(...)
read_SAB(...)
read_CC(...)
read_SC(...)
read_PA(...)
```

Use these when you already know the file family.
They return the same `PRD` object type as `read_auto`.

### Writers

Writer functions:

```python
write_wsr98d(prd, filename, **kwargs)
write_nexrad_level2_msg31(prd, filename, **kwargs)
write_nexrad_level2_msg1(prd, filename, **kwargs)
```

Preferred object-style export helpers:

- `radar.to_wsr98d(...)`
- `radar.to_nexrad_level2_msg31(...)`
- `radar.to_nexrad_level2_msg1(...)`

### Supported families

The public readers target:

- `WSR98D`
- `SAB`
- `CC`
- `SC`
- `PA`
- selected NEXRAD Level II workflows where enabled

If the input is not recognized, `read_auto` raises a format error instead of
silently guessing.

## `pycwr.core.NRadar.PRD`

`PRD` is the central object in `pycwr`. Almost all user workflows start by
reading a file into a `PRD`.

### Core attributes

- `fields`: list of sweep-level `xarray.Dataset`
- `scan_info`: volume metadata as an `xarray.Dataset`
- `extended_fields`: sidecar storage for native-range fields
- `product`: computed product dataset
- `sitename`: site name
- `nsweeps`: number of sweeps
- `nrays`: number of rays
- `effective_earth_radius`: geometry radius used by this volume

### Sweep data model

`fields[sweep]` is usually the best place to inspect one sweep.

Typical coordinates:

- `time`
- `range`
- `azimuth`
- `elevation`
- `x`, `y`, `z`
- `lon`, `lat`

Typical variables:

- `dBZ`
- `V`
- `W`
- `ZDR`
- `CC`
- `PhiDP`
- `KDP`
- corrected fields such as `Zc`, `ZDRc`, `PhiDPc`, `KDPc` when QC has run

### Inspection helpers

```python
radar.summary()
radar.available_fields(sweep=None, range_mode="aligned")
radar.sweep_summary()
radar.get_sweep_field(sweep, field_name, range_mode="aligned", sort_by_azimuth=False)
radar.get_native_sweep_field(sweep, field_name)
radar.has_extended_field(sweep, field_name)
radar.ordered_az(inplace=False)
```

#### `summary()`

Returns a lightweight `dict` with:

- site information
- scan type
- number of sweeps and rays
- field list
- per-sweep summaries

Use this first when opening an unfamiliar file.

#### `available_fields(sweep=None, range_mode="aligned")`

Returns:

- all visible field names for the full volume if `sweep is None`
- field names for one sweep if `sweep` is specified

With `range_mode="native"`, native sidecar fields are included when present.

#### `sweep_summary()`

Returns one summary row per sweep. Common keys include:

- `sweep`
- `fixed_angle`
- `rays`
- `aligned_fields`
- `native_fields`
- `aligned_max_range_m`

#### `get_sweep_field(...)`

```python
radar.get_sweep_field(
    sweep,
    field_name,
    range_mode="aligned",
    sort_by_azimuth=False,
)
```

Returns one `xarray.DataArray`.

Use this when you want:

- one field only
- explicit `aligned` or `native` control
- optional azimuth sorting

#### `get_native_sweep_field(sweep, field_name)`

Returns the native-range field when an extended sidecar exists.
If no native sidecar exists, it falls back to the aligned field.

#### `has_extended_field(sweep, field_name)`

Returns `True` or `False`.
Useful when you want to branch between aligned and native workflows explicitly.

#### `ordered_az(inplace=False)`

Returns or applies an azimuth-sorted view.

- `inplace=False`: return a sorted view
- `inplace=True`: mutate the current object

### Reflectivity access: aligned vs native

For some low sweeps:

- `radar.fields[sweep]["dBZ"]` is the aligned field on the shared range grid
- `radar.get_native_sweep_field(sweep, "dBZ")` is the native long-range field

Recommended rule:

- use aligned access for historical products and compatibility checks
- use native access when the actual low-level reflectivity range matters

Example:

```python
aligned = radar.get_sweep_field(0, "dBZ", range_mode="aligned")
native = radar.get_sweep_field(0, "dBZ", range_mode="native")
```

### Product generation

Public product builders on `PRD` include:

```python
radar.add_product_CR_xy(XRange, YRange, range_mode="aligned")
radar.add_product_CAPPI_xy(XRange, YRange, level_height, range_mode="aligned")
radar.add_product_CAPPI_3d_xy(XRange, YRange, level_heights, range_mode="aligned")
radar.add_product_VIL_xy(XRange, YRange, level_heights, range_mode="aligned")
radar.add_product_ET_xy(XRange, YRange, level_heights, range_mode="aligned")
radar.add_product_CR_lonlat(XLon, YLat, range_mode="aligned")
radar.add_product_CAPPI_lonlat(XLon, YLat, level_height, range_mode="aligned")
radar.add_product_VIL_lonlat(XLon, YLat, level_heights, range_mode="aligned")
radar.add_product_ET_lonlat(XLon, YLat, level_heights, range_mode="aligned")
radar.add_product_VWP(sweeps=None, field_name=None, range_mode="aligned", **kwargs)
```

Units:

- Cartesian grids: meters
- geographic grids: degrees
- heights: meters

Results are written into `radar.product`.

### Section extraction

```python
radar.extract_section(
    start,
    end,
    field_name="dBZ",
    point_units="km",
    interpolation="linear",
    range_mode="aligned",
    sample_spacing=None,
)
radar.extract_section_lonlat(
    start_lonlat,
    end_lonlat,
    field_name="dBZ",
    interpolation="linear",
    range_mode="aligned",
    sample_spacing=None,
)
radar.get_RHI_data(...)
radar.get_vcs_data(...)
```

Use these for vertical sections, lon/lat sections, and RHI-style extraction.

### Wind retrieval on `PRD`

```python
radar.retrieve_vad(sweeps=None, field_name=None, range_mode="aligned", **kwargs)
radar.retrieve_vvp(sweep, field_name=None, range_mode="aligned", **kwargs)
radar.retrieve_vwp(sweeps=None, field_name=None, range_mode="aligned", **kwargs)
radar.add_product_VWP(sweeps=None, field_name=None, range_mode="aligned", **kwargs)
```

Behavior:

- `retrieve_vad(...)`: returns ring-wise `xarray.Dataset` results
- `retrieve_vvp(...)`: returns local wind analysis on one sweep
- `retrieve_vwp(...)`: returns a vertical wind profile `xarray.Dataset`
- `add_product_VWP(...)`: stores the profile in `radar.product` as `VWP_*`

### Export and interop

Common object methods:

- `radar.to_pyart_radar(...)`
- `radar.to_xradar(...)`
- `radar.to_wsr98d(...)`
- `radar.to_nexrad_level2_msg31(...)`
- `radar.to_nexrad_level2_msg1(...)`
- `radar.to_cfgridded_netcdf(...)`

These are the main public export paths for downstream interoperability.

## `pycwr.draw`

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

### Return type

The easy plotting functions return an `EasyPlotResult` with:

- `fig`: matplotlib figure
- `ax`: target axis or axes
- `artist`: main plotted artist

### Common plotting entry points

- `plot_ppi(radar, field="dBZ", sweep=0, ...)`
- `plot_ppi_map(radar, field="dBZ", sweep=0, ...)`
- `plot_rhi(radar, field="dBZ", azimuth=..., ...)`
- `plot_section(radar, start=..., end=..., field="dBZ", ...)`
- `plot_section_lonlat(radar, start_lonlat=..., end_lonlat=..., field="dBZ", ...)`
- `plot_vvp(radar, sweep=0, background_field="dBZ", ...)`
- `plot_wind_profile(profile_or_radar, ...)`

Recommended usage:

- use `plot_ppi` for quick Cartesian PPI inspection
- use `plot_ppi_map` when geographic context matters
- use `plot_section` or `plot_section_lonlat` for vertical analysis
- use `plot_vvp` for vector-field output from wind retrieval
- use `plot_wind_profile` for `VWP` profile products

## `pycwr.qc`

### Main QC workflow

```python
from pycwr.qc import apply_dualpol_qc, run_dualpol_qc
```

Typical usage:

```python
qc_radar = apply_dualpol_qc(radar, inplace=False, clear_air_mode="mask")
```

Purpose:

- despeckle or suppress non-meteorological signals
- generate corrected polarimetric fields
- emit mask fields used by later workflows

Typical output fields include:

- `Zc`
- `ZDRc`
- `PhiDPc`
- `KDPc`
- `QC_MASK`
- `CLEAR_AIR_MASK`

Use `inplace=False` when you want to preserve the raw input volume.

## `pycwr.retrieve`

This module contains two public families of retrieval tools:

- hydrometeor classification
- single-radar wind retrieval

### Hydrometeor classification

Public helpers:

```python
from pycwr.retrieve import (
    apply_hydrometeor_classification,
    classify_hydrometeors,
    interpolate_temperature_profile,
)
```

Typical object workflow:

```python
hcl_radar = radar.classify_hydrometeors(
    inplace=False,
    band="C",
    profile_height=[0.0, 2000.0, 4000.0, 8000.0],
    profile_temperature=[24.0, 12.0, 2.0, -16.0],
    confidence_field="HCL_CONF",
)
```

Use cases:

- classify hydrometeors directly from arrays
- interpolate a temperature profile to gate heights
- write `HCL` and confidence fields back into a `PRD`

### Wind retrieval

Public helpers:

```python
from pycwr.retrieve import retrieve_vad, retrieve_vvp, retrieve_vwp
```

The main algorithms are:

- `VAD`: harmonic fit on one or more sweeps
- `VVP`: local least-squares horizontal wind retrieval on one sweep
- `VWP`: robust vertical profile built from multiple VAD layers

Typical usage:

```python
vad = radar.retrieve_vad(sweeps=[0, 1, 2], max_range_km=40.0, gate_step=4)
vvp = radar.retrieve_vvp(0, max_range_km=20.0, az_num=91, bin_num=5)
vwp = radar.retrieve_vwp(sweeps=[0, 1, 2], max_range_km=40.0, height_step=500.0)
```

Outputs are `xarray.Dataset` objects containing variables such as:

- `u`
- `v`
- `wind_speed`
- `wind_direction`
- `fit_rmse`
- coverage or sample-count metrics depending on the algorithm

Important behavior:

- velocity-field selection prefers `Vc` when present, then falls back to `V`
- retrievals are designed to tolerate missing data and partial azimuth coverage
- `attrs` record the method, the actual velocity field used, and references

Method references included in the module:

- Browning and Wexler (1968), VAD
- Waldteufel and Corbin (1979), VVP

## `pycwr.interp`

Recommended high-level network entry points:

```python
from pycwr.interp import (
    parse_radar_time_from_filename,
    discover_radar_files,
    select_radar_files,
    build_latlon_grid,
    load_network_config,
    run_radar_network_3d,
    radar_network_3d_to_netcdf,
)
```

Typical use:

1. discover or select input radar files
2. build a lon/lat grid
3. run `run_radar_network_3d(...)`
4. optionally write NetCDF

Outputs commonly include:

- network `CR`
- network `CAPPI`
- 3D reflectivity volumes
- per-radar metadata and range summaries

## `pycwr.GraphicalInterface`

Public entry points:

```python
from pycwr.GraphicalInterface import create_app, launch
```

Typical use:

```python
app = create_app()
launch()
```

Or run the script:

```bash
python scripts/LaunchGUI.py
```

The viewer is designed for local use:

- loopback-only binding
- token-guarded API
- restricted file access

## Common workflow recipes

### Read, inspect, and plot

```python
from pycwr.io import read_auto
from pycwr.draw import plot_ppi

radar = read_auto("your_radar_file")
print(radar.summary())
plot_ppi(radar, field="dBZ", sweep=0, show=True)
```

### Use native low-level reflectivity

```python
native_dBZ = radar.get_sweep_field(0, "dBZ", range_mode="native")
```

### Run QC and plot corrected reflectivity

```python
qc_radar = radar.apply_dualpol_qc(inplace=False)
plot_ppi(qc_radar, field="Zc", sweep=0, show=True)
```

### Build a wind profile and plot it

```python
from pycwr.draw import plot_wind_profile

profile = radar.retrieve_vwp(sweeps=[0, 1, 2], max_range_km=40.0, height_step=500.0)
plot_wind_profile(profile, show=True)
```

### Run multi-radar compositing

```python
from pycwr.interp import run_radar_network_3d

network = run_radar_network_3d(...)
```

## What to read next

- [../README.md](../README.md): project overview
- [../test/README.md](../test/README.md): runnable examples
- [draw_quickstart.md](draw_quickstart.md): plotting entry points
- [web_viewer_quickstart.md](web_viewer_quickstart.md): local viewer setup
