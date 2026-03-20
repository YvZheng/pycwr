# pycwr

A Python toolkit for operational Chinese weather radar data, covering reading, analysis, plotting, quality control, multi-radar compositing, and export.

- Current release: `1.0.0`
- [中文说明](README_CN.md)
- [API reference](docs/api_reference.md)
- [Test guide](test/README.md)
- [Draw quickstart](docs/draw_quickstart.md)
- [Web viewer quickstart](docs/web_viewer_quickstart.md)
- [Contributors](CONTRIBUTORS.txt)

## What's New In 1.0.0

This is a major cleanup-and-hardening release focused on compatibility, geometry correctness, and lower operational complexity.

Module highlights:

- `pycwr.io`: lighter default install while keeping the core radar readers available
- `pycwr.core`: stricter geometry path with round-trip validation for radar/cartesian/geographic transforms
- `pycwr.draw` and `GraphicalInterface`: optional dependencies isolated from the base import path
- `pycwr.qc`: clearer example-oriented test coverage for dual-pol and X-band workflows
- `pycwr.retrieve`: packaged hydrometeor classification with profile-aware and profile-free workflows
- `pycwr.interp`: better-documented network compositing entry points and sample-driven tests

Performance highlights:

- redundant geometry implementations were removed so the main path now stays on one high-precision formula chain
- vectorized Cartesian conversion hot paths were tightened to reduce intermediate arrays and repeated angle-grid construction
- packaging now prefers a smaller default dependency set, which improves cross-platform installation success and lowers environment friction

## What pycwr does

`pycwr` follows one simple main workflow:

1. Read a radar file with `pycwr.io.read_auto`
2. Get back a `PRD` volume object
3. Use that volume for plotting, product generation, QC, network interpolation, or export

Core capabilities:

- Read `WSR98D`, `SAB`, `CC`, `SC`, and `PA`
- Preserve the historical aligned workflow while exposing native long-range reflectivity when available
- Plot PPI, map PPI, RHI, and vertical sections
- Compute `CR`, `CAPPI`, and 3D network products
- Run dual-polarization QC
- Run hydrometeor classification (`HCL`) with or without a temperature profile
- Export to Py-ART, xradar, and NetCDF
- Launch a local web viewer

## Installation

Install the default dependencies and the base package:

```bash
python -m pip install -r requirements-core.txt
python -m pip install .
```

Install the full-featured environment:

```bash
python -m pip install -r requirements-full.txt
python -m pip install ".[full]"
```

Rebuild the Cython extension after editing `pycwr/core/RadarGridC.pyx`:

```bash
python setup.py build_ext --inplace
```

Build source and wheel distributions:

```bash
python scripts/build_dist.py
```

The helper prefers `python -m build` when available and falls back to setuptools otherwise.

## 5-Minute Quickstart

Read one radar file:

```python
from pycwr.io import read_auto

radar = read_auto("Z_RADR_I_Z9046_20260317065928_O_DOR_SAD_CAP_FMT.bin.bz2")
print(radar.summary())
```

Inspect the returned volume:

```python
print(radar.available_fields())
print(radar.sweep_summary()[0])
print(radar.fields[0]["dBZ"])
```

Make a minimal plot:

```python
from pycwr.draw import plot_ppi

plot_ppi(radar, field="dBZ", sweep=0, show=True)
```

Extract a vertical section:

```python
from pycwr.draw import plot_section

plot_section(radar, start=(-50, 0), end=(50, 0), field="dBZ", show=True)
```

## Core Concepts

### `PRD`

All readers return `pycwr.core.NRadar.PRD`.

The most important members are:

- `fields`: one `xarray.Dataset` per sweep
- `scan_info`: site and scan metadata
- `extended_fields`: sidecar storage for native long-range fields
- `product`: computed products such as `CR`, `CAPPI`, and 3D grids

The most useful methods for new users:

- `summary()`: compact summary of the full volume
- `available_fields(sweep=None)`: list available fields for the whole volume or one sweep
- `sweep_summary()`: one summary record per sweep
- `get_sweep_field(sweep, field_name, range_mode="aligned" | "native")`
- `ordered_az(inplace=False)`

### Aligned Reflectivity vs Native Reflectivity

For some low sweeps, `pycwr` intentionally keeps two reflectivity views:

- `radar.fields[sweep]["dBZ"]`: the aligned result used by historical workflows
- `radar.get_native_sweep_field(sweep, "dBZ")`: the native long-range reflectivity

Which one should you use?

- Use the aligned field when you need strict compatibility with older processing chains
- Use the native field when you need the full low-level reflectivity coverage

## Public Interfaces

### 1. Reading radar data

```python
from pycwr.io import (
    read_auto,
    read_WSR98D,
    read_SAB,
    read_CC,
    read_SC,
    read_PA,
    write_wsr98d,
    write_nexrad_level2_msg31,
    write_nexrad_level2_msg1,
)
```

All readers share the same public signature:

```python
read_xxx(
    filename,
    station_lon=None,
    station_lat=None,
    station_alt=None,
    effective_earth_radius=None,
)
```

### 2. Plotting

Recommended entry points:

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

These functions return an `EasyPlotResult` containing `fig`, `ax`, and `artist`.

### 3. Product generation

Common `PRD` methods:

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
- `get_RHI_data`
- `get_vcs_data`

### 4. Quality control

```python
from pycwr.qc import apply_dualpol_qc, run_dualpol_qc
```

You can also call QC directly on a `PRD` object:

```python
processed = radar.apply_dualpol_qc(inplace=False)
```

### 5. Multi-radar 3D workflow

```python
from pycwr.interp import run_radar_network_3d
```

This workflow builds 3D radar network products on a regular lon/lat grid and can write NetCDF directly.

### 6. Hydrometeor classification

```python
from pycwr.retrieve import classify_hydrometeors, interpolate_temperature_profile
```

You can classify one sweep from arrays, or write `HCL` directly back into a `PRD`.

### 7. Wind retrieval

```python
from pycwr.retrieve import retrieve_vad, retrieve_vvp, retrieve_vwp
```

You can retrieve:

- ring-wise VAD winds from one or more sweeps
- local VVP winds from one sweep
- a VWP profile either directly or via `PRD.add_product_VWP(...)`

### 8. Web viewer

Launch the local viewer from Python:

```python
from pycwr.GraphicalInterface import create_app, launch
```

Or use the helper script:

```bash
python scripts/LaunchGUI.py
```

Current viewer boundaries:

- It only listens on loopback addresses
- All API requests require a token
- File access is restricted to configured allowed roots

## Common Workflows

### Plot one PPI sweep

```python
from pycwr.io import read_auto
from pycwr.draw import plot_ppi

radar = read_auto("your_radar_file")
plot_ppi(radar, field="dBZ", sweep=0, show=True)
```

### Compute CR, CAPPI, VIL, and ET

```python
import numpy as np

x = np.arange(-150_000.0, 150_001.0, 1_000.0)
y = np.arange(-150_000.0, 150_001.0, 1_000.0)
radar.add_product_CR_xy(x, y)
radar.add_product_CAPPI_xy(x, y, 3000.0)
radar.add_product_VIL_xy(x, y, np.array([1000.0, 2000.0, 3000.0]))
radar.add_product_ET_xy(x, y, np.array([1000.0, 2000.0, 3000.0]))

print(radar.product)
```

### Run QC and use corrected fields

```python
qc_radar = radar.apply_dualpol_qc(inplace=False, clear_air_mode="mask")
plot_ppi(qc_radar, field="Zc", sweep=0, show=True)
```

Notes:

- `CLEAR_AIR_MASK` is written as a dedicated weak-echo diagnostic for likely clear-air / biological echoes.
- `clear_air_mode="label"` keeps legacy QC acceptance and only labels likely clear-air echoes.
- `clear_air_mode="mask"` removes likely clear-air echoes from `QC_MASK`, `Zc`, `KDPc`, and downstream QC products.
- This clear-air branch is intentionally conservative and is not a full hydrometeor classification algorithm.

### Run hydrometeor classification

With an environmental temperature profile:

```python
hcl_radar = radar.classify_hydrometeors(
    inplace=False,
    band="C",
    profile_height=[0.0, 2000.0, 4000.0, 8000.0, 12000.0],
    profile_temperature=[24.0, 12.0, 2.0, -16.0, -40.0],
    confidence_field="HCL_CONF",
    temperature_field="HCL_T",
)
plot_ppi(hcl_radar, field="HCL", sweep=0, show=True)
```

Without a temperature profile:

```python
radar.add_hydrometeor_classification(
    band="C",
    confidence_field="HCL_CONF",
)
```

Notes:

- `HCL` is added as a gate-level sweep field, not a Cartesian product.
- The classifier uses the aligned dual-pol range shared by `dBZ/ZDR/KDP/CC/LDR`.
- If `profile_height/profile_temperature` are omitted, `pycwr` falls back to a reduced-variable no-profile mode.
- Class ids follow the packaged 10-class warm-season fuzzy HID table: drizzle, rain, ice crystals, dry aggregates snow, wet snow, vertical ice, low-density graupel, high-density graupel, hail, big drops.

### Retrieve single-radar winds

```python
from pycwr.retrieve import retrieve_vad, retrieve_vvp
from pycwr.draw import plot_vvp, plot_wind_profile

vad = retrieve_vad(radar, sweeps=0, max_range_km=40.0, gate_step=4)
vvp = retrieve_vvp(radar, sweep=0, max_range_km=20.0, az_num=91, bin_num=9)
vwp = radar.retrieve_vwp(sweeps=[0, 1, 2], max_range_km=40.0, gate_step=4, height_step=500.0)
radar.add_product_VWP(sweeps=[0, 1, 2], max_range_km=40.0, gate_step=4, height_step=500.0)

plot_vvp(vvp, show=True)
plot_wind_profile(vwp, show=True)
```

Notes:

- `retrieve_vad(...)` returns a sweep-by-gate VAD dataset.
- `retrieve_vvp(...)` returns a single-sweep local VVP analysis.
- `retrieve_vwp(...)` aggregates VAD output into a vertical wind profile.
- `add_product_VWP(...)` stores the profile into `radar.product` as `VWP_*` variables.

### Build a 3D network mosaic

```python
from pycwr.interp import run_radar_network_3d

dataset = run_radar_network_3d(
    target_time="2026-03-17T07:00:00",
    radar_dirs=["./Z9046", "./Z9047"],
    lon_min=118.5,
    lon_max=119.5,
    lat_min=31.8,
    lat_max=32.8,
    lon_res_deg=0.01,
    lat_res_deg=0.01,
    level_heights=[1000.0, 3000.0, 5000.0],
    field_names=["dBZ"],
    output_path="./network.nc",
)
```

Notes:

- `dBZ`-class scalar fields are supported for 3D mosaics.
- Velocity fields are intentionally not supported yet in `run_radar_network_3d(...)`.
- `VIL` / `ET` can use a denser `product_level_heights` grid than the display `level_heights`.
- When `ET` is requested, the output dataset also includes `ET_TOPPED`, which flags columns whose highest sampled level still exceeds the echo-top threshold.
- `VIL` uses the Greene and Clark (1972) liquid-water relation with a 56 dBZ cap; the low-reflectivity cutoff remains configurable and should be treated as a workflow assumption.

Algorithm references used by `pycwr.interp.RadarInterp`:

- VIL: Greene, D. R., and R. A. Clark, 1972, Monthly Weather Review, doi:10.1175/1520-0493(1972)100<0548:VILWNA>2.3.CO;2
- VIL operational guidance: NOAA WDTD, https://vlab.noaa.gov/web/wdtd/-/vertically-integrated-liquid-vil-
- Echo tops: Lakshmanan et al., 2013, Weather and Forecasting, doi:10.1175/WAF-D-12-00084.1
- Echo tops operational guidance: NOAA WDTD, https://vlab.noaa.gov/web/wdtd/-/xx-dbz-echo-top-et-
- Enhanced Echo Tops metadata semantics: WSR-88D ROC ICD 2620003Y, https://www.roc.noaa.gov/public-documents/icds/2620003Y.pdf
- Radar mosaic overlap handling: 吴翀;双偏振雷达的资料质量分析,相态识別及组网应用[D];南京信息工程大学;2018年, Chapter 4; and Lakshmanan et al., 2006, Weather and Forecasting, doi:10.1175/WAF942.1

QC references used by `pycwr.qc`:

- Dual-Pol QC operational guidance: NOAA WDTD, https://vlab.noaa.gov/web/wdtd/-/dual-pol-quality-control
- MRMS QC updates: Tang et al., 2020, Journal of Atmospheric and Oceanic Technology, doi:10.1175/JTECH-D-19-0165.1
- Clear-air and non-meteorological echo treatment in the local radar workflow context: 吴翀;双偏振雷达的资料质量分析,相态识別及组网应用[D];南京信息工程大学;2018年

### Export to Py-ART / xradar

```python
pyart_radar = radar.to_pyart_radar()
sweeps = radar.to_xradar_sweeps()
tree = radar.to_xradar(strict=True)
```

### Export to WSR98D / NEXRAD Level II

```python
radar.to_wsr98d("./roundtrip_wsr98d.bin", overwrite=True)
radar.to_nexrad_level2_msg31("./export_msg31.ar2v", overwrite=True)
radar.to_nexrad_level2_msg1("./export_msg1.ar2v", overwrite=True)
```

Notes:

- `to_wsr98d(...)` writes a pycwr-readable WSR-98D base-data file
- `to_nexrad_level2_msg31(...)` writes a Py-ART-readable NEXRAD Level II MSG31 archive
- `to_nexrad_level2_msg1(...)` writes a Py-ART-readable NEXRAD Level II MSG1 archive
- export currently supports `ppi` volumes only
- MSG31 supports `dBZ/V/W/ZDR/CC/PhiDP`
- MSG1 is intentionally limited to `dBZ/V/W`

## Units

Internal core conventions:

- Distance and height: meters
- Internal trigonometric calculations: radians
- Lon/lat input and output: degrees

Historical public compatibility is still preserved:

- Public `azimuth`, `elevation`, and `fixed_angle` remain in degrees
- Some legacy plotting interfaces still use kilometers

## Tests and Sample Data

Run the automated test suite:

```bash
python3 -m unittest discover -s test -p 'test_*.py'
```

Recommended test groups:

- [test/README.md](test/README.md): full index of example and regression scripts
- `test_examples_public_api.py`: quickstart-style integration tests on real sample files
- `test_examples_qc.py`: QC primitives and pipeline tests
- `test_examples_hid.py`: hydrometeor classification examples with and without profiles
- `test_examples_interp.py`: multi-radar interpolation and compositing examples
- `test_regression_geometry.py`: geometry round-trip validation
- `test_regression_security.py`: path-boundary and malformed-input regression tests

If the expected sample files are not available in the current workspace, some sample-based tests are skipped automatically.
Slow real-sample 3D network tests are opt-in via `PYCWR_RUN_SLOW_SAMPLE_TESTS=1`.

## Documentation Map

- [docs/api_reference.md](docs/api_reference.md): public modules, functions, and common signatures
- [docs/draw_quickstart.md](docs/draw_quickstart.md): plotting quickstart
- [docs/web_viewer_quickstart.md](docs/web_viewer_quickstart.md): local viewer usage

## License

This repository is distributed under **PolyForm Noncommercial 1.0.0**.
Commercial use, paid redistribution, and profit-making deployment are not permitted by default.
Contact the project separately if a commercial license is required.
See [LICENSE.txt](LICENSE.txt).
