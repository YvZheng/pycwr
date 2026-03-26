# pycwr 接口参考

这份文档面向 `pycwr 1.0.3` 的公开接口。
它不是源码逐行索引，而是告诉你：

- 应该从哪个模块进入
- 每个接口大致接收什么参数
- 返回的对象是什么
- aligned 和 native 反射率工作流该怎么区分

如果你想看按功能整理的可运行示例，请继续看
[../test/README.md](../test/README.md)。

## 模块布局

| 模块 | 主要用途 | 推荐入口 |
| --- | --- | --- |
| `pycwr.io` | 读取和写出雷达基数据 | `read_auto`, `read_WSR98D`, `read_SAB`, `read_CC`, `read_SC`, `read_PA` |
| `pycwr.core` | 核心体扫对象、几何和导出辅助 | `PRD`, `radar.summary()`, `radar.get_sweep_field()` |
| `pycwr.draw` | 绘图和快速出图 | `plot_ppi`, `plot_ppi_map`, `plot_rhi`, `plot_section`, `plot_vvp`, `plot_wind_profile` |
| `pycwr.qc` | 双偏振质量控制 | `apply_dualpol_qc`, `run_dualpol_qc` |
| `pycwr.retrieve` | 水凝物分类和风场反演 | `classify_hydrometeors`, `retrieve_vad`, `retrieve_vvp`, `retrieve_vwp` |
| `pycwr.interp` | 多雷达组网插值 | `run_radar_network_3d`, `radar_network_3d_to_netcdf` |
| `pycwr.GraphicalInterface` | 本地 Web viewer | `create_app`, `launch` |

## 通用约定

### 单位

- 距离和高度：米
- 内部三角函数计算：弧度
- 对外 `azimuth`、`elevation`、`fixed_angle`：度
- 经纬度输入输出：度

部分历史绘图接口仍然接受千米，这是为了兼容旧行为。

### `range_mode`

很多接口都接受：

- `range_mode="aligned"`：历史兼容的对齐工作流
- `range_mode="native"`：有原生长距离反射率时走原生距离库

建议：

- 做历史兼容产品和对比时，用 `aligned`
- 需要低层真实反射率覆盖范围时，用 `native`

速度场通常仍保持 aligned，因为速度有效距离往往短于反射率。

### 可选依赖

基础安装包含：

- reader
- `PRD`
- 几何
- 插值
- NetCDF 风格导出

full 安装额外包含：

- 绘图
- 地图绘图
- QC
- Web viewer
- Py-ART 和 xradar 互操作

安装 full 版本：

```bash
pip install "pycwr[full]"
```

## `pycwr.io`

### 推荐入口：`read_auto`

```python
read_auto(
    filename,
    station_lon=None,
    station_lat=None,
    station_alt=None,
    effective_earth_radius=None,
)
```

作用：

- 自动识别雷达文件家族
- 解析成 `PRD`
- 保持项目当前兼容的 sweep 布局和几何主链

参数：

- `filename`：雷达文件路径
- `station_lon`、`station_lat`、`station_alt`：可选站点覆盖值
- `effective_earth_radius`：可选有效地球半径，单位米

返回：

- `pycwr.core.NRadar.PRD`

典型用法：

```python
from pycwr.io import read_auto

radar = read_auto("your_radar_file.bin.bz2")
print(radar.summary())
```

### 指定格式 reader

```python
read_WSR98D(...)
read_SAB(...)
read_CC(...)
read_SC(...)
read_PA(...)
```

只有在你已经明确知道文件格式时，才建议直接调用这些 reader。
它们返回的仍然是同一个 `PRD` 对象类型。

### 写出接口

函数式 writer：

```python
write_wsr98d(prd, filename, **kwargs)
write_nexrad_level2_msg31(prd, filename, **kwargs)
write_nexrad_level2_msg1(prd, filename, **kwargs)
```

更推荐的对象式导出：

- `radar.to_wsr98d(...)`
- `radar.to_nexrad_level2_msg31(...)`
- `radar.to_nexrad_level2_msg1(...)`

### 当前支持的雷达家族

公开 reader 主要面向：

- `WSR98D`
- `SAB`
- `CC`
- `SC`
- `PA`
- 启用时的部分 NEXRAD Level II 工作流

如果格式无法识别，`read_auto` 会直接报格式错误，不会静默猜测。

## `pycwr.core.NRadar.PRD`

`PRD` 是整个项目的中心对象。几乎所有用户工作流都会先得到一个 `PRD`。

### 核心属性

- `fields`：每层一个 `xarray.Dataset`
- `scan_info`：体扫元数据 `xarray.Dataset`
- `extended_fields`：原生距离 sidecar
- `product`：产品结果集
- `sitename`：站点名
- `nsweeps`：层数
- `nrays`：总射线数
- `effective_earth_radius`：几何计算使用的有效地球半径

### sweep 数据模型

`fields[sweep]` 往往是查看单层数据最直接的入口。

常见坐标：

- `time`
- `range`
- `azimuth`
- `elevation`
- `x`、`y`、`z`
- `lon`、`lat`

常见变量：

- `dBZ`
- `V`
- `W`
- `ZDR`
- `CC`
- `PhiDP`
- `KDP`
- 运行 QC 后生成的 `Zc`、`ZDRc`、`PhiDPc`、`KDPc`

### 体扫查看接口

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

返回一个轻量 `dict`，通常包含：

- 站点信息
- 扫描类型
- 层数和射线数
- 场列表
- 每层摘要

第一次打开陌生数据时，建议先看这个。

#### `available_fields(sweep=None, range_mode="aligned")`

返回：

- `sweep is None` 时，返回整个体扫可见场名
- 指定 `sweep` 时，返回某一层可见场名

如果 `range_mode="native"`，会把 sidecar 中可见的原生场也纳入结果。

#### `sweep_summary()`

按层返回摘要列表。常见字段包括：

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

返回一个 `xarray.DataArray`。

适合：

- 只取一个场
- 显式控制 `aligned` / `native`
- 顺手做按方位排序

#### `get_native_sweep_field(sweep, field_name)`

如果某层存在原生 sidecar，就返回原生距离库上的场；
如果没有，则自动回退到 aligned。

#### `has_extended_field(sweep, field_name)`

返回 `True` / `False`。
适合显式区分 aligned 与 native 路径。

#### `ordered_az(inplace=False)`

返回或应用“按 azimuth 排序后的视图”。

- `inplace=False`：返回排序视图
- `inplace=True`：直接修改当前对象

### 反射率访问：aligned 与 native

对部分低层：

- `radar.fields[sweep]["dBZ"]` 是对齐版
- `radar.get_native_sweep_field(sweep, "dBZ")` 是原生长距离反射率

建议规则：

- 做历史兼容产品、图和对比时，用 aligned
- 做低层完整反射率分析时，用 native

示例：

```python
aligned = radar.get_sweep_field(0, "dBZ", range_mode="aligned")
native = radar.get_sweep_field(0, "dBZ", range_mode="native")
```

### 产品计算

`PRD` 上的公开产品接口包括：

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

单位：

- 笛卡尔网格：米
- 经纬度网格：度
- 高度：米

结果会写回 `radar.product`。

### 剖面提取

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

适合做垂直剖面、经纬度剖面和 RHI 风格提取。

### 在 `PRD` 上做风场反演

```python
radar.retrieve_vad(sweeps=None, field_name=None, range_mode="aligned", **kwargs)
radar.retrieve_vvp(sweep, field_name=None, range_mode="aligned", **kwargs)
radar.retrieve_vwp(sweeps=None, field_name=None, range_mode="aligned", **kwargs)
radar.add_product_VWP(sweeps=None, field_name=None, range_mode="aligned", **kwargs)
```

行为：

- `retrieve_vad(...)`：返回环状风反演 `xarray.Dataset`
- `retrieve_vvp(...)`：返回单层局地水平风反演结果
- `retrieve_vwp(...)`：返回垂直风廓线 `xarray.Dataset`
- `add_product_VWP(...)`：把廓线写入 `radar.product`，变量名是 `VWP_*`

### 导出与互操作

常用对象方法：

- `radar.to_pyart_radar(...)`
- `radar.to_xradar(...)`
- `radar.to_wsr98d(...)`
- `radar.to_nexrad_level2_msg31(...)`
- `radar.to_nexrad_level2_msg1(...)`
- `radar.to_cfgridded_netcdf(...)`

这些是下游互操作的主要公开出口。

## `pycwr.draw`

推荐使用的公开绘图接口：

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

### 返回值

easy 绘图函数统一返回 `EasyPlotResult`，其中包含：

- `fig`：matplotlib figure
- `ax`：目标轴或轴数组
- `artist`：主绘图对象

### 常用绘图入口

- `plot_ppi(radar, field="dBZ", sweep=0, ...)`
- `plot_ppi_map(radar, field="dBZ", sweep=0, ...)`
- `plot_rhi(radar, field="dBZ", azimuth=..., ...)`
- `plot_section(radar, start=..., end=..., field="dBZ", ...)`
- `plot_section_lonlat(radar, start_lonlat=..., end_lonlat=..., field="dBZ", ...)`
- `plot_vvp(radar, sweep=0, background_field="dBZ", ...)`
- `plot_wind_profile(profile_or_radar, ...)`

建议：

- 快速看 PPI，用 `plot_ppi`
- 需要地理背景，用 `plot_ppi_map`
- 做垂直分析，用 `plot_section` 或 `plot_section_lonlat`
- 看 VVP 矢量场，用 `plot_vvp`
- 看 `VWP` 结果，用 `plot_wind_profile`

## `pycwr.qc`

### 主 QC 工作流

```python
from pycwr.qc import apply_dualpol_qc, run_dualpol_qc
```

典型用法：

```python
qc_radar = apply_dualpol_qc(radar, inplace=False, clear_air_mode="mask")
```

作用：

- 去除或抑制非气象目标
- 生成订正后的极化场
- 输出后续流程可复用的掩码场

常见输出场：

- `Zc`
- `ZDRc`
- `PhiDPc`
- `KDPc`
- `QC_MASK`
- `CLEAR_AIR_MASK`

如果你想保留原始体扫，建议用 `inplace=False`。

## `pycwr.retrieve`

这个模块包含两类公开反演能力：

- 水凝物分类
- 单雷达风场反演

### 水凝物分类

公开入口：

```python
from pycwr.retrieve import (
    apply_hydrometeor_classification,
    classify_hydrometeors,
    interpolate_temperature_profile,
)
```

典型对象式工作流：

```python
hcl_radar = radar.classify_hydrometeors(
    inplace=False,
    band="C",
    profile_height=[0.0, 2000.0, 4000.0, 8000.0],
    profile_temperature=[24.0, 12.0, 2.0, -16.0],
    confidence_field="HCL_CONF",
)
```

适用场景：

- 直接对数组做分类
- 把温度廓线插值到 gate 高度
- 把 `HCL` 和置信度场写回 `PRD`

### 风场反演

公开入口：

```python
from pycwr.retrieve import retrieve_vad, retrieve_vvp, retrieve_vwp
```

三条主算法：

- `VAD`：一个或多个 sweep 的谐波拟合
- `VVP`：单层 sweep 的局地最小二乘水平风反演
- `VWP`：由多层 VAD 聚合的稳健垂直风廓线

典型用法：

```python
vad = radar.retrieve_vad(sweeps=[0, 1, 2], max_range_km=40.0, gate_step=4)
vvp = radar.retrieve_vvp(0, max_range_km=20.0, az_num=91, bin_num=5)
vwp = radar.retrieve_vwp(sweeps=[0, 1, 2], max_range_km=40.0, height_step=500.0)
```

返回值都是 `xarray.Dataset`，常见变量包括：

- `u`
- `v`
- `wind_speed`
- `wind_direction`
- `fit_rmse`
- 不同算法对应的覆盖率或样本数指标

重要行为：

- 速度场选择优先 `Vc`，没有时回退到 `V`
- 反演会尽量容忍缺测和不完整方位覆盖
- `attrs` 中会记录方法、实际使用的速度场和参考文献

模块中写入的主要参考文献：

- Browning and Wexler (1968), VAD
- Waldteufel and Corbin (1979), VVP

## `pycwr.interp`

推荐的高层组网入口：

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

典型流程：

1. 发现或筛选输入雷达文件
2. 构建经纬度网格
3. 执行 `run_radar_network_3d(...)`
4. 按需写 NetCDF

常见输出包括：

- 组网 `CR`
- 组网 `CAPPI`
- 3D 反射率体
- 每部雷达的元数据和距离摘要

## `pycwr.GraphicalInterface`

公开入口：

```python
from pycwr.GraphicalInterface import create_app, launch
```

典型使用：

```python
app = create_app()
launch()
```

或者直接运行：

```bash
python scripts/LaunchGUI.py
```

viewer 的设计边界：

- 只监听本机回环地址
- API 受 token 保护
- 文件访问受目录限制

## 常见工作流配方

### 读数、查看、出一张图

```python
from pycwr.io import read_auto
from pycwr.draw import plot_ppi

radar = read_auto("your_radar_file")
print(radar.summary())
plot_ppi(radar, field="dBZ", sweep=0, show=True)
```

### 使用低层原生反射率

```python
native_dBZ = radar.get_sweep_field(0, "dBZ", range_mode="native")
```

### 做 QC 后画订正反射率

```python
qc_radar = radar.apply_dualpol_qc(inplace=False)
plot_ppi(qc_radar, field="Zc", sweep=0, show=True)
```

### 生成风廓线并绘图

```python
from pycwr.draw import plot_wind_profile

profile = radar.retrieve_vwp(sweeps=[0, 1, 2], max_range_km=40.0, height_step=500.0)
plot_wind_profile(profile, show=True)
```

### 跑多雷达组网

```python
from pycwr.interp import run_radar_network_3d

network = run_radar_network_3d(...)
```

## 接下来读什么

- [../README_CN.md](../README_CN.md)：项目总览
- [../test/README.md](../test/README.md)：可运行示例
- [draw_quickstart.md](draw_quickstart.md)：绘图入口
- [web_viewer_quickstart.md](web_viewer_quickstart.md)：本地 viewer 使用说明
