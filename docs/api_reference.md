# pycwr API Reference

This page is a practical guide to the public interfaces most users need.
It is intentionally selective: the goal is to help you find the right API,
understand what it expects, and know what it returns.

## How To Read This Document

`pycwr` is easiest to learn in this order:

1. Read radar data with `pycwr.io`
2. Inspect the returned `PRD` object
3. Plot or extract sections from `PRD`
4. Compute products or run QC
5. Run hydrometeor classification if you need gate-level phase labels
6. Export or build multi-radar network products

If you are new to the package, start with `read_auto`, `PRD.summary()`,
`PRD.fields`, and `pycwr.draw.plot_ppi`.
For runnable examples grouped by feature, see [../test/README.md](../test/README.md).

## Package Layout

| Module | Main purpose | Typical entry points |
| --- | --- | --- |
| `pycwr.io` | Read and write radar base data | `read_auto`, `read_WSR98D`, `read_SAB`, `read_CC`, `read_SC`, `read_PA`, `write_wsr98d`, `write_nexrad_level2_msg31`, `write_nexrad_level2_msg1` |
| `pycwr.core` | Core radar object, geometry, products | `PRD`, `grid_3d_network_xy` |
| `pycwr.draw` | Plotting | `plot`, `plot_ppi`, `plot_ppi_map`, `plot_rhi`, `plot_section`, `plot_section_lonlat`, `plot_vvp`, `plot_wind_profile` |
| `pycwr.qc` | Dual-pol QC primitives and workflow | `apply_dualpol_qc`, `run_dualpol_qc` |
| `pycwr.retrieve` | Hydrometeor and wind retrieval helpers | `apply_hydrometeor_classification`, `classify_hydrometeors`, `interpolate_temperature_profile`, `retrieve_vad`, `retrieve_vvp`, `retrieve_vwp`, `VAD`, `VVP` |
| `pycwr.interp` | Multi-radar network compositing | `parse_radar_time_from_filename`, `discover_radar_files`, `select_radar_files`, `build_latlon_grid`, `load_network_config`, `run_radar_network_3d`, `radar_network_3d_to_netcdf` |
| `pycwr.GraphicalInterface` | Lightweight local web viewer | `create_app`, `launch` |

## Cross-Cutting Conventions

### Units

Core internal rules:

- Distance and height: meters
- Internal trigonometric calculations: radians
- Lon/lat input and output: degrees

Historical public compatibility is still preserved:

- Public `azimuth`, `elevation`, and `fixed_angle` remain in degrees
- Some legacy plotting interfaces still use kilometers where that was the historical behavior

### `range_mode`

Several APIs accept:

- `range_mode="aligned"`: historical aligned workflow
- `range_mode="native"`: native long-range reflectivity where available

Use `aligned` when you need backward-compatible behavior.
Use `native` when you need the full low-level reflectivity range.

### Optional dependencies

Base install:

- reading
- core geometry
- `PRD`
- interpolation
- NetCDF export

Optional full stack:

- plotting
- map plotting
- web viewer
- QC
- Py-ART / xradar interop

Install the full stack with:

```bash
pip install "pycwr[full]"
```

## `pycwr.io`

This is the first API most users touch.

### `read_auto`

```python
read_auto(
    filename,
    station_lon=None,
    station_lat=None,
    station_alt=None,
    effective_earth_radius=None,
)
```

Use this first. It detects the radar family automatically and returns a `PRD`.

Parameters:

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

Use these only when you already know the file family.

They share the same main arguments as `read_auto` and also return `PRD`.

### Writers

```python
write_wsr98d(prd, filename, **kwargs)
write_nexrad_level2_msg31(prd, filename, **kwargs)
write_nexrad_level2_msg1(prd, filename, **kwargs)
```

For most user code, the `PRD` methods are the preferred entry points:

- `radar.to_wsr98d(...)`
- `radar.to_nexrad_level2_msg31(...)`
- `radar.to_nexrad_level2_msg1(...)`

### Supported file families

Current public readers target:

- `WSR98D`
- `SAB`
- `CC`
- `SC`
- `PA`

If a file is not recognized, `read_auto` raises a format error instead of silently guessing.

## `pycwr.core.NRadar.PRD`

`PRD` is the central data object in `pycwr`.
You usually do not instantiate it directly; you get it from `pycwr.io`.

### Core attributes

- `fields`: list of sweep-level `xarray.Dataset`
- `scan_info`: volume metadata as an `xarray.Dataset`
- `extended_fields`: sidecar storage for native-range fields
- `product`: computed product dataset
- `sitename`: site name
- `nsweeps`: number of sweeps
- `nrays`: number of rays
- `effective_earth_radius`: beam-geometry radius used by this volume

### Sweep data model

`fields[sweep]` is usually the most useful place to inspect one sweep.

Typical coordinates and variables:

- coordinates: `time`, `range`, `azimuth`, `elevation`, `x`, `y`, `z`, `lon`, `lat`
- fields: `dBZ`, `V`, `W`, `ZDR`, `CC`, `PhiDP`, `KDP`, and derived fields when available

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

Recommended usage:

- `summary()`: start here when inspecting a new file
- `available_fields()`: list fields quickly
- `sweep_summary()`: understand each sweep without dumping full arrays
- `get_sweep_field()`: stable access helper for aligned or native mode
- `ordered_az()`: get an azimuth-sorted view for analysis or plotting

#### `summary()`

Returns a lightweight `dict` with:

- site information
- scan type
- number of sweeps and rays
- available fields
- per-sweep summaries

#### `available_fields(sweep=None, range_mode="aligned")`

Returns:

- all available field names for the full volume if `sweep is None`
- field names for one sweep if `sweep` is specified

With `range_mode="native"`, native sidecar fields are included when present.

#### `sweep_summary()`

Returns a list of `dict` rows, one per sweep.

Each row includes:

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

Use this when:

- you want one field only
- you want `aligned` vs `native` chosen explicitly
- you want optional azimuth sorting

#### `get_native_sweep_field(sweep, field_name)`

Returns the field on its native range grid when native sidecar data exists.
If no native sidecar exists, it falls back to the aligned field.

#### `has_extended_field(sweep, field_name)`

Returns `True` or `False`.
Useful when you want to branch between aligned and native workflows explicitly.

#### `ordered_az(inplace=False)`

Returns a view where azimuth becomes the sweep dimension order.

Behavior:

- `inplace=False`: return a sorted view
- `inplace=True`: mutate the current object

### Reflectivity access: aligned vs native

This is the most important data-access rule in the package.

For some low sweeps:

- `radar.fields[sweep]["dBZ"]` is the aligned field used by the historical workflow
- `radar.get_native_sweep_field(sweep, "dBZ")` is the native long-range reflectivity field

Recommended rule:

- use `fields` or `get_sweep_field(..., range_mode="aligned")` for backward-compatible plotting and products
- use `get_native_sweep_field()` or `get_sweep_field(..., range_mode="native")` when full low-level coverage matters

### Product generation

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
- Geographic grids: degrees
- Heights: meters

Results are written into `radar.product`.

Common naming:

- aligned products keep their historical names
- native-range products use `_native` suffixes internally

### Wind retrieval on a PRD volume

```python
radar.retrieve_vad(sweeps=None, field_name=None, range_mode="aligned", **kwargs)
radar.retrieve_vvp(sweep, field_name=None, range_mode="aligned", **kwargs)
radar.retrieve_vwp(sweeps=None, field_name=None, range_mode="aligned", **kwargs)
radar.add_product_VWP(sweeps=None, field_name=None, range_mode="aligned", **kwargs)
```

Use these when you want single-radar wind diagnostics without leaving the `PRD` workflow.

Behavior:

- `retrieve_vad(...)` returns a sweep-by-gate ring-wind analysis
- `retrieve_vvp(...)` returns a local horizontal-wind analysis on one sweep
- `retrieve_vwp(...)` aggregates VAD output into a vertical wind profile
- `add_product_VWP(...)` stores the profile in `radar.product` as `VWP_*` variables

### Section and RHI extraction

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
radar.extract_rhi(azimuth=None, field_name="dBZ", range_mode="aligned")
radar.get_RHI_data(az, field_name="dBZ", range_mode="aligned")
radar.get_vcs_data(start_point, end_point, field_name, range_mode="aligned")
```

Guidance:

- `extract_section`: Cartesian start/end points
- `extract_section_lonlat`: lon/lat start/end points
- `extract_rhi`: standardized RHI dataset
- `get_RHI_data`: older lower-level arrays
- `get_vcs_data`: older lower-level vertical-section arrays

Recommended new code should prefer:

- `extract_section`
- `extract_section_lonlat`
- `extract_rhi`

### QC, HID, and export

```python
radar.apply_dualpol_qc(
    inplace=False,
    band="C",
    use_existing_kdp=True,
    clear_air_mode="label",
)
radar.classify_hydrometeors(
    inplace=False,
    band="C",
    profile_height=None,
    profile_temperature=None,
    confidence_field=None,
)
radar.add_hydrometeor_classification(band="C", confidence_field=None)
radar.to_pyart_radar(range_mode="aligned", field_names=None, use_external=None, strict=False)
radar.to_xradar_sweeps(range_mode="aligned", field_names=None)
radar.to_xradar(range_mode="aligned", field_names=None, strict=True)
radar.to_wsr98d(filename, field_names=None, strict=True, overwrite=False)
radar.to_nexrad_level2_msg31(filename, field_names=None, strict=True, overwrite=False)
radar.to_nexrad_level2_msg1(filename, field_names=None, strict=True, overwrite=False)
```

#### `apply_dualpol_qc(...)`

Applies dual-pol QC to selected sweeps and returns:

- a copy if `inplace=False`
- the same object if `inplace=True`

New fields commonly added:

- `Zc`
- `PIA`
- `PhiDP_smooth`
- `PhiDP_texture`
- `KDPc`
- `METEO_MASK`
- `CLEAR_AIR_MASK`
- `QC_MASK`
- `ZDRc`
- `PIA_ZDR`

#### `classify_hydrometeors(...)`

Adds gate-level hydrometeor class labels to selected sweeps and returns:

- a copy if `inplace=False`
- the same object if `inplace=True`

Important parameters:

- `profile_height` / `profile_temperature`: optional 1-D environmental profile
- `temperature`: optional direct gate-temperature field
- `confidence_field`: optional max-membership confidence output
- `temperature_field`: optional field storing the temperature field actually used

Behavior:

- writes `HCL` into the aligned dual-pol sweep range
- automatically prefers corrected fields like `Zc`, `ZDRc`, and `KDPc` when present
- falls back to a reduced-variable no-profile mode when temperature is omitted

#### `add_hydrometeor_classification(...)`

Convenience in-place wrapper for `classify_hydrometeors(...)`.

#### `to_pyart_radar(...)`

Use this when you want a Py-ART-style radar object.

Important parameters:

- `range_mode`: `aligned` or `native`
- `field_names`: export only selected fields if needed
- `use_external`: prefer upstream Py-ART or force bundled compatibility object
- `strict=True`: raise a clear error if the optional dependency is missing

#### `to_xradar_sweeps(...)`

Returns sweep-level `xarray.Dataset` objects suitable for xradar-style workflows.

#### `to_xradar(...)`

Returns a DataTree-style object when the optional stack is available.

#### `to_wsr98d(...)`

Writes a WSR-98D base-data file that can be read back by `pycwr.io.read_WSR98D(...)`.

#### `to_nexrad_level2_msg31(...)`

Writes a Py-ART-readable NEXRAD Level II MSG31 archive.
Supported moments are `dBZ`, `V`, `W`, `ZDR`, `CC`, and `PhiDP`.

#### `to_nexrad_level2_msg1(...)`

Writes a Py-ART-readable NEXRAD Level II MSG1 archive.
This legacy format is intentionally limited to `dBZ`, `V`, and `W`.

## `pycwr.draw`

This module provides beginner-friendly high-level plotting helpers.

Recommended imports:

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
    EasyRadarPlotter,
)
```

### Return type

All high-level plotting helpers return:

```python
EasyPlotResult(fig=..., ax=..., artist=...)
```

This makes it easy to:

- inspect the `matplotlib` figure or axes
- add extra annotations
- save figures yourself if needed

### `plot_ppi`

```python
plot_ppi(
    radar,
    field="dBZ",
    sweep=0,
    ax=None,
    figsize=(8, 8),
    show=False,
    save=None,
    title=None,
    **kwargs,
)
```

Use this as the default quicklook API for one sweep.

### `plot_ppi_map`

```python
plot_ppi_map(
    radar,
    field="dBZ",
    sweep=0,
    ax=None,
    figsize=(8, 8),
    projection=None,
    data_crs=None,
    show=False,
    save=None,
    title=None,
    **kwargs,
)
```

Use this when you want the same sweep on a map.

Requires the optional full plotting stack.

### `plot_section`

```python
plot_section(
    radar,
    start=(-50.0, 0.0),
    end=(50.0, 0.0),
    field="dBZ",
    point_units="km",
    ...
)
```

Use Cartesian section endpoints.
`point_units` can be `km` or `m`.

### `plot_section_lonlat`

```python
plot_section_lonlat(
    radar,
    start_lonlat,
    end_lonlat,
    field="dBZ",
    ...
)
```

Use lon/lat section endpoints.

### `plot_rhi`

```python
plot_rhi(
    radar,
    azimuth,
    field="dBZ",
    ...
)
```

Use this to extract and plot an RHI-style vertical slice at a target azimuth.

### `plot`

```python
plot(radar, kind="ppi", field="dBZ", sweep=0, show=False, save=None, **kwargs)
```

Supported `kind` values:

- `ppi`
- `ppi_map`
- `section`
- `section_lonlat`
- `rhi`
- `wind_profile`
- `vvp`

### `EasyRadarPlotter`

Small wrapper class for users who prefer:

```python
plotter = EasyRadarPlotter(radar)
plotter.ppi(...)
plotter.section(...)
plotter.rhi(...)
plotter.wind_profile(...)
plotter.vvp(...)
```

## `pycwr.qc`

This module provides dual-polarization QC.

### High-level workflow

```python
apply_dualpol_qc(
    prd,
    sweeps=None,
    inplace=False,
    band="C",
    use_existing_kdp=True,
    clear_air_mode="label",
)
```

Use this when you already have a `PRD` object.

Requirements:

- `dBZ` must exist
- at least one of `KDP` or `PhiDP` must exist

Behavior:

- computes corrected reflectivity and related QC diagnostics
- writes new fields into the selected sweeps
- returns a new object unless `inplace=True`

### Low-level array workflow

```python
run_dualpol_qc(
    ref,
    zdr=None,
    phidp=None,
    kdp=None,
    rhohv=None,
    snr=None,
    dr=0.075,
    band="C",
    use_existing_kdp=True,
    clear_air_mode="label",
)
```

Use this if you already work with arrays outside a `PRD`.

Returns a `dict` with arrays such as:

- `ref_corrected`
- `zdr_corrected`
- `pia`
- `pia_zdr`
- `phidp_smooth`
- `kdp_used`
- `phidp_texture`
- `meteo_mask`
- `clear_air_mask`
- `qc_mask`

QC reference notes:

- NOAA WDTD dual-pol QC: https://vlab.noaa.gov/web/wdtd/-/dual-pol-quality-control
- Tang et al. 2020 MRMS QC updates: doi:10.1175/JTECH-D-19-0165.1
- 吴翀;双偏振雷达的资料质量分析,相态识別及组网应用[D];南京信息工程大学;2018年. In the local workflow context, it treats clear-air echo as a non-meteorological class to be separated during QC.

### Lower-level QC helpers

Frequently useful functions:

```python
smooth_phidp(...)
kdp_from_phidp(...)
phidp_texture(...)
build_meteo_mask(...)
build_clear_air_mask(...)
despeckle_mask(...)
correct_attenuation_kdp(...)
pia_from_kdp(...)
```

Use them when you want to build your own QC sequence instead of using the packaged workflow.

## `pycwr.retrieve`

This module provides hydrometeor classification and single-radar wind-retrieval helpers.

### Hydrometeor classification on a PRD volume

```python
apply_hydrometeor_classification(
    prd,
    sweeps=None,
    inplace=False,
    band="C",
    method="hybrid",
    temperature=None,
    profile_height=None,
    profile_temperature=None,
    confidence_field=None,
    temperature_field=None,
)
```

Use this when you already have a `PRD` and want an `HCL` field written back into the sweeps.

Requirements:

- `dBZ` must exist
- at least one of `ZDR`, `KDP`, `CC`, or `LDR` must exist

Behavior:

- writes `HCL` as a gate-level sweep field
- optionally writes confidence and interpolated temperature fields
- uses the packaged 10-class fuzzy HID table
- if temperature is omitted, falls back to a reduced-variable no-profile mode

### Array-level hydrometeor classification

```python
classify_hydrometeors(
    dBZ,
    ZDR=None,
    KDP=None,
    CC=None,
    LDR=None,
    T=None,
    method="hybrid",
    band="C",
    weights=None,
    return_scores=False,
    return_confidence=False,
)
```

Use this when you already have arrays outside a `PRD`.

Returns:

- class ids in `[1, 10]`
- optional fuzzy-score cube when `return_scores=True`
- optional winning-class confidence when `return_confidence=True`

### Interpolate an environmental profile to gate temperature

```python
interpolate_temperature_profile(
    gate_altitude,
    profile_height,
    profile_temperature,
    radar_altitude=0.0,
    height_reference="asl",
)
```

Use `dataset["z"]` from a PRD sweep as `gate_altitude`.

### Hydrometeor class ids

`pycwr.retrieve.available_hydrometeor_classes()` returns:

1. `Drizzle`
2. `Rain`
3. `Ice Crystals`
4. `Dry Aggregates Snow`
5. `Wet Snow`
6. `Vertical Ice`
7. `Low-Density Graupel`
8. `High-Density Graupel`
9. `Hail`
10. `Big Drops`

References used by the packaged classifier:

- Dolan et al. (2013), *Journal of Applied Meteorology and Climatology*
- Marzano et al. (2006), *Advances in Geosciences*

### Wind retrieval helpers

```python
retrieve_vad(
    prd,
    sweeps=None,
    field_name=None,
    range_mode="aligned",
    gate_step=1,
    max_range_km=None,
    **kwargs,
)
retrieve_vwp(
    prd,
    sweeps=None,
    field_name=None,
    range_mode="aligned",
    height_bins=None,
    height_step=250.0,
    max_height_m=None,
    gate_step=1,
    max_range_km=None,
    **kwargs,
)
retrieve_vvp(
    prd,
    sweep,
    field_name=None,
    range_mode="aligned",
    az_num=91,
    bin_num=9,
    azimuth_step=1,
    range_step=1,
    max_range_km=None,
    **kwargs,
)
```

Use these when you want wind retrieval results as datasets rather than writing them back into `radar.product`.

Returns:

- `retrieve_vad(...)`: sweep-by-gate VAD dataset
- `retrieve_vwp(...)`: 1-D vertical wind profile dataset
- `retrieve_vvp(...)`: local single-sweep VVP dataset

Legacy compatibility helpers are still exported:

- `vad` / `VAD`
- `vvp` / `VVP`

## `pycwr.interp`

This module provides multi-radar network interpolation and NetCDF export.

### File discovery and time selection

```python
parse_radar_time_from_filename(path)
discover_radar_files(radar_dirs, pattern="*.bin*")
select_radar_files(radar_dirs, target_time, tolerance_minutes=10, pattern="*.bin*")
```

#### `parse_radar_time_from_filename(path)`

Parses scan time from a radar filename.
Raises `ValueError` when the filename does not contain a recognizable timestamp.

#### `discover_radar_files(...)`

Returns a dictionary keyed by radar directory name with candidate files.

#### `select_radar_files(...)`

Returns a list of records, one per selected radar file, containing:

- `radar_id`
- `path`
- `scan_time`
- `delta_seconds`

### Grid construction

```python
build_latlon_grid(lon_min, lon_max, lat_min, lat_max, lon_res_deg, lat_res_deg)
```

Returns:

- `lon`
- `lat`
- `grid_lon`
- `grid_lat`

### End-to-end 3D network workflow

```python
run_radar_network_3d(
    target_time,
    config_path=None,
    radar_dirs=None,
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    lon_res_deg=None,
    lat_res_deg=None,
    level_heights=None,
    field_names=None,
    output_products=None,
    output_path=None,
    time_tolerance_minutes=None,
    composite_method=None,
    influence_radius_m=None,
    fillvalue=None,
    effective_earth_radius=None,
    field_range_mode=None,
    max_range_km=None,
    parallel=None,
    max_workers=None,
    blind_method=None,
    use_qc=None,
    qc_method=None,
    qc_band=None,
    qc_use_existing_kdp=None,
    qc_fallback=None,
    ...
)
```

Use this when you want a complete multi-radar network product.

Minimum inputs you must define:

- `target_time`
- `radar_dirs`
- lon/lat bounds and resolution
- `level_heights`

Returns:

- one `xarray.Dataset`

Important groups of arguments:

- data selection: `target_time`, `radar_dirs`, `time_tolerance_minutes`, `pattern`
- grid definition: `lon_min`, `lon_max`, `lat_min`, `lat_max`, `lon_res_deg`, `lat_res_deg`, `level_heights`
- compositing: `field_names`, `field_range_mode`, `composite_method`, `influence_radius_m`, `max_range_km`
- QC: `use_qc`, `qc_method`, `qc_band`, `qc_use_existing_kdp`, `qc_fallback`
- output: `output_products`, `output_path`
- execution: `parallel`, `max_workers`

### NetCDF export

```python
radar_network_3d_to_netcdf(dataset, output_path)
```

Writes the dataset and returns the output path as a string.

## `pycwr.GraphicalInterface`

This is the lightweight local web viewer API.

### `create_app`

```python
create_app(allowed_roots=None, auth_token=None)
```

Parameters:

- `allowed_roots`: directories the viewer may browse and open
- `auth_token`: fixed token for automation or tests; if omitted, a random token is generated

Returns:

- a configured Flask app

### `launch`

```python
launch(host="127.0.0.1", port=8787, open_browser=True)
```

Notes:

- only loopback hosts are allowed
- non-index API calls require a token
- file access is restricted to allowed roots

## Recommended Learning Order

1. `pycwr.io.read_auto`
2. `PRD.summary()` and `PRD.fields`
3. `PRD.get_sweep_field(...)`
4. `pycwr.draw.plot_ppi`
5. `PRD.add_product_CR_xy` / `PRD.add_product_CAPPI_xy`
6. `PRD.apply_dualpol_qc`
7. `pycwr.interp.run_radar_network_3d`

## Related Docs

- [docs/draw_quickstart.md](draw_quickstart.md)
- [docs/web_viewer_quickstart.md](web_viewer_quickstart.md)
- [README.md](../README.md)
