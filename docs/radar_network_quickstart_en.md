# PyCWR Radar Network Quickstart

This guide is for first-time users who want to build a multi-radar network product with `pycwr`.

The goal is simple:

- combine multiple radars into one 3D network volume
- derive common products such as `CR`, `VIL`, and `ET`
- save the result as NetCDF for later analysis or plotting

If you only want the shortest working example, jump to section 3.

---

## 1. What you will use

The recommended high-level entrypoints are:

```python
from pycwr.interp import run_radar_network_3d, radar_network_3d_to_netcdf
```

In most cases, you do not need to grid each radar manually and then combine them yourself.
`run_radar_network_3d(...)` already handles the full workflow.

It will:

1. find the nearest radar file for each radar directory
2. read the files and interpolate them onto one common lat/lon grid
3. compose a 3D network reflectivity volume
4. optionally derive `CR`, `VIL`, and `ET`
5. optionally write the output NetCDF file

---

## 2. Recommended directory layout

The best beginner-friendly layout is one directory per radar:

```text
/data/radar_network/
├── Z9001/
│   ├── Z_RADR_I_Z9001_20260317070000_O_DOR_....bin.bz2
│   └── Z_RADR_I_Z9001_20260317070600_O_DOR_....bin.bz2
├── Z9002/
│   ├── Z_RADR_I_Z9002_20260317065930_O_DOR_....bin.bz2
│   └── Z_RADR_I_Z9002_20260317070530_O_DOR_....bin.bz2
└── Z9003/
    └── Z_RADR_I_Z9003_20260317070010_O_DOR_....bin.bz2
```

Important notes:

- `radar_dirs` is a list of directories, not a list of files
- filenames must contain a parseable `YYYYMMDDHHMMSS` timestamp
- each radar directory may contain many files; `pycwr` picks the one nearest to `target_time`

---

## 3. Shortest working example

```python
from pycwr.interp import run_radar_network_3d

dataset = run_radar_network_3d(
    target_time="2026-03-17T07:00:00",
    radar_dirs=[
        "/data/radar_network/Z9001",
        "/data/radar_network/Z9002",
        "/data/radar_network/Z9003",
    ],
    lon_min=118.0,
    lon_max=122.0,
    lat_min=29.0,
    lat_max=32.0,
    lon_res_deg=0.01,
    lat_res_deg=0.01,
    level_heights=[1000.0, 3000.0, 5000.0, 7000.0],
    field_names=["dBZ"],
    output_products=["CR", "VIL", "ET"],
    parallel=False,
    output_path="./radar_network_20260317070000.nc",
)

print(dataset)
print(dataset.data_vars)
print(dataset.attrs)
```

This produces:

- a 3D network reflectivity volume: `dBZ`
- a 2D composite reflectivity product: `CR`
- a 2D vertically integrated liquid product: `VIL`
- a 2D echo-top product: `ET`
- a NetCDF file: `./radar_network_20260317070000.nc`

For a first run, keep `parallel=False` so debugging is easier.

---

## 4. What comes back

The return value is an `xarray.Dataset`.

Common data variables:

- `dBZ`: 3D network reflectivity, usually with dimensions `("z", "lat", "lon")`
- `CR`: 2D composite reflectivity
- `VIL`: 2D vertically integrated liquid
- `ET`: 2D echo-top height
- `coverage_count`: how many radars can really observe each grid point
- `blind_mask`: blind-zone mask

Common coordinates:

- `z`
- `lat`
- `lon`

Common attributes:

- `composite_method`
- `field_range_mode`
- `output_products`
- `qc_applied_by_radar`
- `plot_files`

Example:

```python
level0 = dataset["dBZ"].isel(z=0)
print(level0.shape)
```

You can also reopen the file with `xarray`:

```python
import xarray as xr

ds = xr.load_dataset("./radar_network_20260317070000.nc")
print(ds)
```

---

## 5. Key parameters that beginners should understand

### `target_time`

The target scan time. You can pass:

- a `datetime`
- an ISO string such as `"2026-03-17T07:00:00"`

`pycwr` will choose the nearest file in each radar directory.

### `radar_dirs`

A list of radar directories:

```python
radar_dirs=[
    "/data/radar_network/Z9001",
    "/data/radar_network/Z9002",
]
```

### `lon_min / lon_max / lat_min / lat_max`

The geographic extent of the network domain, in degrees.

### `lon_res_deg / lat_res_deg`

Grid spacing in degrees.

- `0.01` is a reasonable first choice
- smaller values create denser grids and higher memory usage

### `level_heights`

Vertical levels in meters:

```python
level_heights=[1000.0, 3000.0, 5000.0, 7000.0]
```

### `field_names`

Fields to mosaic:

```python
field_names=["dBZ"]
```

You can also request more than one field:

```python
field_names=["dBZ", "ZDR"]
```

If one radar does not carry a requested field, `pycwr` skips that radar for that field.
For the first run, start with `["dBZ"]`.

### `output_products`

Additional derived products, usually:

- `CR`
- `VIL`
- `ET`

### `output_path`

If provided, the dataset is written automatically to NetCDF.

---

## 6. A practical recommended template

```python
from pycwr.interp import run_radar_network_3d

dataset = run_radar_network_3d(
    target_time="2026-03-17T07:00:00",
    radar_dirs=[
        "/data/radar_network/Z9001",
        "/data/radar_network/Z9002",
        "/data/radar_network/Z9003",
    ],
    lon_min=118.0,
    lon_max=122.0,
    lat_min=29.0,
    lat_max=32.0,
    lon_res_deg=0.01,
    lat_res_deg=0.01,
    level_heights=[1000.0, 3000.0, 5000.0, 7000.0],
    field_names=["dBZ"],
    output_products=["CR", "VIL", "ET"],
    composite_method="quality_weighted",
    max_range_km=460.0,
    blind_method="hybrid",
    parallel=False,
    output_path="./network_20260317070000.nc",
)
```

This is a safe starting template for real use.

---

## 7. Check which files will be selected

Before you run the full workflow, you can inspect candidate files:

```python
from pycwr.interp import discover_radar_files, select_radar_files

radar_dirs = [
    "/data/radar_network/Z9001",
    "/data/radar_network/Z9002",
]

files = discover_radar_files(radar_dirs, pattern="*.bin*")
print(files)

selected = select_radar_files(
    radar_dirs,
    target_time="2026-03-17T07:00:00",
    tolerance_minutes=10,
    pattern="*.bin*",
)
print(selected)
```

This is useful when you need to debug:

- why one radar was not selected
- whether the filename timestamp is parseable
- whether the time tolerance is too small

---

## 8. Config-file workflow

`pycwr` can also read network options from:

- `.json`
- `.toml`
- `.yml`
- `.yaml`

For beginners, YAML is usually the easiest.

Example `network.yml`:

```yaml
network:
  radar_dirs:
    - /data/radar_network/Z9001
    - /data/radar_network/Z9002
    - /data/radar_network/Z9003
  field_names: [dBZ]
  field_range_mode:
    dBZ: native
  lon_min: 118.0
  lon_max: 122.0
  lat_min: 29.0
  lat_max: 32.0
  lon_res_deg: 0.01
  lat_res_deg: 0.01
  level_heights: [1000.0, 3000.0, 5000.0, 7000.0]
  output_products: [CR, VIL, ET]
  max_range_km: 460.0
  composite_method: quality_weighted
  blind_method: hybrid
  parallel: false
  output_dir: ./output
```

Then run:

```python
from pycwr.interp import run_radar_network_3d

dataset = run_radar_network_3d(
    target_time="2026-03-17T07:00:00",
    config_path="./network.yml",
)
```

This is useful when:

- you have many parameters
- you want to rerun the same domain many times
- you want to hand the setup to another user

---

## 9. Optional QC before compositing

You can ask `pycwr` to run QC before compositing:

```python
dataset = run_radar_network_3d(
    target_time="2026-03-17T07:00:00",
    radar_dirs=[
        "/data/radar_network/Z9001",
        "/data/radar_network/Z9002",
    ],
    lon_min=118.0,
    lon_max=122.0,
    lat_min=29.0,
    lat_max=32.0,
    lon_res_deg=0.01,
    lat_res_deg=0.01,
    level_heights=[1000.0, 3000.0],
    field_names=["dBZ"],
    use_qc=True,
    parallel=False,
)
```

Notes:

- if `use_qc=True`, reflectivity may be sourced from corrected fields such as `Zc`
- the dataset attributes record whether QC was actually applied by radar
- for a first run, it is usually better to get the base network workflow working first

---

## 10. Common beginner mistakes

### `No radar files matched the target time within tolerance`

Usually means:

- the target time is wrong
- the filename does not contain a parseable timestamp
- the time tolerance is too small
- `radar_dirs` points to the wrong directories

### The run is too slow or memory-heavy

Usually caused by:

- a very large lon/lat domain
- very fine `lon_res_deg` / `lat_res_deg`
- too many `level_heights`

Start small first.

### One field is missing in the result

Usually means:

- some radars do not carry that field
- the requested `field_names` value is wrong
- QC changed the actual source field used in the workflow

### Only `dBZ` exists, but no `CR/VIL/ET`

Because you did not request:

```python
output_products=["CR", "VIL", "ET"]
```

---

## 11. Recommended learning order

If this is your first time:

1. make sure the directory layout is one directory per radar
2. run `discover_radar_files` / `select_radar_files` first
3. start with `field_names=["dBZ"]`
4. start with a small domain and only a few height levels
5. keep `parallel=False`
6. add `CR/VIL/ET` after the base workflow works
7. enable `use_qc=True` only after that

---

## 12. Read next

- [Chinese API reference](api_reference_cn.md)
- [English API reference](api_reference.md)
- [Draw quickstart](draw_quickstart.md)
- [README](../README.md)
- [README_CN](../README_CN.md)
- [Read the Docs](https://pycwr.readthedocs.io/en/latest/)

If you plan to hand this workflow to complete beginners, it is worth preparing one extra internal example with your own station directories and your own standard network domain.
