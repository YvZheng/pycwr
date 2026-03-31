# PyCWR 雷达组网快速上手

这份文档是给第一次使用 `pycwr` 做多雷达组网的人准备的。
目标很简单：

- 先把多部雷达数据拼成一个 3D 组网结果
- 再输出常见产品，比如 `CR`、`VIL`、`ET`
- 最后把结果保存成 NetCDF，方便后续分析或出图

如果你只想先跑通一遍，直接看下面的“最短可跑示例”。

---

## 1. 先知道你会用到什么

`pycwr` 当前推荐的高层组网入口是：

```python
from pycwr.interp import run_radar_network_3d, radar_network_3d_to_netcdf
```

你通常不需要自己逐部雷达做 3D 插值再手工合成，直接调用 `run_radar_network_3d(...)` 就行。

它会帮你做这些事情：

1. 从每个雷达目录里找接近目标时次的文件
2. 读取雷达文件并插值到统一经纬度网格
3. 合成 3D 组网反射率体
4. 按需生成 `CR`、`VIL`、`ET`
5. 按需写出 NetCDF

---

## 2. 你的数据目录应该怎么放

最推荐的目录结构是“每部雷达一个目录”：

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

有几个重要要求：

- `radar_dirs` 传的是“目录列表”，不是文件列表
- 文件名里要能解析出时间，格式里要包含 `YYYYMMDDHHMMSS`
- 每个目录里可以放多个时次，`pycwr` 会自动挑最接近 `target_time` 的文件

---

## 3. 最短可跑示例

下面这个例子适合第一次跑：

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

这个例子会输出：

- 一个 3D 组网反射率体：`dBZ`
- 一个二维组合反射率：`CR`
- 一个二维垂直积分液态含水量：`VIL`
- 一个二维回波顶高：`ET`
- 一个 NetCDF 文件：`./radar_network_20260317070000.nc`

第一次跑建议先把 `parallel=False` 写死，方便排查问题。

---

## 4. 结果里会有什么

返回值是一个 `xarray.Dataset`。

最常见的数据变量包括：

- `dBZ`：三维组网反射率，维度通常是 `("z", "lat", "lon")`
- `CR`：二维组合反射率
- `VIL`：二维垂直积分液态含水量
- `ET`：二维回波顶高
- `coverage_count`：每个格点有多少部雷达真实可观测
- `blind_mask`：盲区掩码

最常见的坐标包括：

- `z`
- `lat`
- `lon`

最常见的属性包括：

- `composite_method`
- `field_range_mode`
- `output_products`
- `qc_applied_by_radar`
- `plot_files`

你可以这样看一个高度层：

```python
level0 = dataset["dBZ"].isel(z=0)
print(level0.shape)
```

也可以直接用 `xarray` 打开刚写出的结果：

```python
import xarray as xr

ds = xr.load_dataset("./radar_network_20260317070000.nc")
print(ds)
```

---

## 5. 最重要的参数怎么理解

下面这些参数，是小白第一次最容易搞混的。

### `target_time`

目标时次。可以传：

- `datetime`
- ISO 时间字符串，比如 `"2026-03-17T07:00:00"`

`pycwr` 会从每个雷达目录中挑最接近这个时次的文件。

### `radar_dirs`

每部雷达的目录列表。

例如：

```python
radar_dirs=[
    "/data/radar_network/Z9001",
    "/data/radar_network/Z9002",
]
```

### `lon_min / lon_max / lat_min / lat_max`

定义组网区域范围，单位是经纬度。

例如：

- `lon_min=118.0`
- `lon_max=122.0`
- `lat_min=29.0`
- `lat_max=32.0`

### `lon_res_deg / lat_res_deg`

网格分辨率，单位是度。

例如：

- `0.01` 大约是 1 km 量级
- 数值越小，格点越密，计算越慢，内存占用越大

第一次建议先用：

```python
lon_res_deg=0.01
lat_res_deg=0.01
```

### `level_heights`

垂直层高度，单位是米。

例如：

```python
level_heights=[1000.0, 3000.0, 5000.0, 7000.0]
```

这是组网 3D 体的高度层。

### `field_names`

要组网的场名列表。

最常见的写法：

```python
field_names=["dBZ"]
```

如果你还想同时组网别的场，可以写：

```python
field_names=["dBZ", "ZDR"]
```

注意：

- 某一部雷达如果缺少某个场，`pycwr` 会跳过这部雷达在这个场上的贡献
- 所以第一次建议先只做 `dBZ`

### `output_products`

要附带计算的产品。

当前常见值：

- `CR`
- `VIL`
- `ET`

例如：

```python
output_products=["CR", "VIL", "ET"]
```

### `output_path`

如果你传了这个参数，函数跑完会自动把结果写成 NetCDF。

例如：

```python
output_path="./network.nc"
```

如果你不传，就只返回内存里的 `xarray.Dataset`。

---

## 6. 一个更实用的推荐模板

如果你不是做测试，而是准备开始真正使用，建议从下面这个模板抄：

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

这是一个相对稳妥的起点。

---

## 7. 如果想先看看每部雷达会选中哪个文件

可以先用这些辅助函数检查：

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

这一步非常适合排查：

- 为什么某部雷达没被选上
- 文件名时间是不是不规范
- 时间容差是不是太小

---

## 8. 想用配置文件也可以

`pycwr` 支持从配置文件读组网参数，支持这些格式：

- `.json`
- `.toml`
- `.yml`
- `.yaml`

最适合小白的是 `YAML`。

例如，新建一个 `network.yml`：

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

然后在 Python 里这样跑：

```python
from pycwr.interp import run_radar_network_3d

dataset = run_radar_network_3d(
    target_time="2026-03-17T07:00:00",
    config_path="./network.yml",
)
```

这种方式适合：

- 参数很多，不想全写在函数里
- 要反复跑不同时间但同一组网配置
- 要交给别人重复使用

---

## 9. 想在组网前先做 QC 怎么办

可以。

`run_radar_network_3d(...)` 支持在组网前对单雷达先做 QC：

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

需要注意：

- 开启 `use_qc=True` 后，反射率可能会自动改用订正后的 `Zc`
- 结果属性里会记录每部雷达是否成功应用了 QC
- 如果你只是第一次试跑，建议先不要开 QC，先保证最基础组网能跑通

---

## 10. 组网后怎么快速看结果

最简单的方法：

```python
print(dataset)
print(list(dataset.data_vars))
```

看二维产品：

```python
cr = dataset["CR"]
vil = dataset["VIL"]
et = dataset["ET"]
```

看某个高度层的 3D 反射率：

```python
dbz_3km = dataset["dBZ"].sel(z=3000.0)
print(dbz_3km.values.shape)
```

---

## 11. 小白最常见的报错原因

### 1）`No radar files matched the target time within tolerance`

说明没有找到满足目标时次的文件。常见原因：

- `target_time` 写错了
- 文件名里没有正确时间
- `time_tolerance_minutes` 太小
- `radar_dirs` 指到了错误目录

### 2）组网区域太大，运行很慢或内存很高

常见原因：

- `lon_min/lon_max/lat_min/lat_max` 范围过大
- `lon_res_deg/lat_res_deg` 太细
- `level_heights` 层数太多

第一次建议小范围、低层数先试。

### 3）某个字段没有结果

常见原因：

- 某些雷达没有这个字段
- 你写的 `field_names` 不存在
- 开启 QC 后，实际参与组网的是订正字段

### 4）结果里只有 `dBZ`，没有 `CR/VIL/ET`

因为你没有显式传：

```python
output_products=["CR", "VIL", "ET"]
```

---

## 12. 给新手的最终建议

如果你是第一次上手，请按下面顺序来：

1. 先确认目录结构是“每部雷达一个目录”
2. 先用 `discover_radar_files` / `select_radar_files` 检查选文件是否正常
3. 第一次只做 `field_names=["dBZ"]`
4. 第一次只做少量高度层
5. 第一次先 `parallel=False`
6. 跑通以后再加 `CR/VIL/ET`
7. 最后再考虑 `use_qc=True`

---

## 13. 你下一步还可以看什么

- [api_reference_cn.md](api_reference_cn.md)：完整接口参考
- [draw_quickstart.md](draw_quickstart.md)：绘图怎么用
- [README_CN.md](../README_CN.md)：项目整体说明
- `test/test_examples_interp.py`：项目里的真实组网示例和回归测试

如果你想把这份文档交给完全没写过 Python 的同事，建议你再配一份“你们自己的目录路径示例”和“你们常用的网格范围模板”，这样上手会更快。
