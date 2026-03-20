# pycwr 接口参考

这份文档面向“第一次真正上手 `pycwr` 的人”。
它不是源码逐行索引，而是告诉你：

- 该从哪个模块进入
- 关键函数的参数是什么
- 返回值长什么样
- 什么场景该用哪个接口

## 如何阅读这份文档

建议按下面顺序学习：

1. 先用 `pycwr.io` 读取数据
2. 看 `PRD` 对象
3. 先画图、看剖面
4. 再做产品计算和 QC
5. 需要相态标签时再做水凝物分类
6. 最后再做导出和组网

如果你第一次碰这个项目，最推荐的起点是：

- `read_auto`
- `PRD.summary()`
- `PRD.fields`
- `pycwr.draw.plot_ppi`

如果你想直接看“按功能整理的可运行示例”，请继续看 [../test/README.md](../test/README.md)。

## 模块总览

| 模块 | 作用 | 常用入口 |
| --- | --- | --- |
| `pycwr.io` | 读取和写出雷达基数据 | `read_auto`, `read_WSR98D`, `read_SAB`, `read_CC`, `read_SC`, `read_PA`, `write_wsr98d`, `write_nexrad_level2_msg31`, `write_nexrad_level2_msg1` |
| `pycwr.core` | 核心雷达对象、几何和产品计算 | `PRD`, `grid_3d_network_xy` |
| `pycwr.draw` | 绘图 | `plot`, `plot_ppi`, `plot_ppi_map`, `plot_rhi`, `plot_section`, `plot_section_lonlat`, `plot_vvp`, `plot_wind_profile` |
| `pycwr.qc` | 双偏振质控 | `apply_dualpol_qc`, `run_dualpol_qc` |
| `pycwr.retrieve` | 水凝物分类和风场反演辅助接口 | `apply_hydrometeor_classification`, `classify_hydrometeors`, `interpolate_temperature_profile`, `retrieve_vad`, `retrieve_vvp`, `retrieve_vwp`, `VAD`, `VVP` |
| `pycwr.interp` | 多雷达组网插值 | `parse_radar_time_from_filename`, `discover_radar_files`, `select_radar_files`, `build_latlon_grid`, `load_network_config`, `run_radar_network_3d`, `radar_network_3d_to_netcdf` |
| `pycwr.GraphicalInterface` | 本地轻量 Web viewer | `create_app`, `launch` |

## 通用约定

### 单位

核心内部统一约定：

- 距离和高度：米
- 内部三角函数计算：弧度
- 经纬度输入输出：度

历史兼容仍然保留：

- 对外 `azimuth`、`elevation`、`fixed_angle` 仍然是度
- 一些旧绘图接口仍然沿用千米

### `range_mode`

很多接口都接受：

- `range_mode="aligned"`：历史兼容的对齐工作流
- `range_mode="native"`：原生长距离反射率工作流

建议：

- 要与旧流程保持一致，用 `aligned`
- 要保留低层完整反射率范围，用 `native`

### 可选依赖

默认安装提供：

- 读数
- 核心几何
- `PRD`
- 插值
- NetCDF 导出

全功能可选依赖提供：

- 绘图
- 地图绘图
- Web viewer
- QC
- Py-ART / xradar 互操作

安装全功能版本：

```bash
pip install "pycwr[full]"
```

## `pycwr.io`

这是绝大多数用户的第一个入口。

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

最推荐的入口。自动识别雷达格式并返回 `PRD`。

参数：

- `filename`：雷达文件路径
- `station_lon`、`station_lat`、`station_alt`：站点信息覆盖值
- `effective_earth_radius`：波束几何使用的有效地球半径，单位米

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

只有你已经明确知道文件格式时，才建议直接调用这些 reader。

它们和 `read_auto` 共用同一组主要参数，也同样返回 `PRD`。

### 写出接口

```python
write_wsr98d(prd, filename, **kwargs)
write_nexrad_level2_msg31(prd, filename, **kwargs)
write_nexrad_level2_msg1(prd, filename, **kwargs)
```

对大多数用户代码，更推荐使用 `PRD` 自己的方法：

- `radar.to_wsr98d(...)`
- `radar.to_nexrad_level2_msg31(...)`
- `radar.to_nexrad_level2_msg1(...)`

### 当前支持的雷达家族

对外支持：

- `WSR98D`
- `SAB`
- `CC`
- `SC`
- `PA`

如果格式无法识别，`read_auto` 会明确报错，而不是静默猜测。

## `pycwr.core.NRadar.PRD`

`PRD` 是整个项目的中心数据对象。
通常你不会手工实例化它，而是从 `pycwr.io` 的 reader 拿到它。

### 核心属性

- `fields`：每层一个 sweep 级 `xarray.Dataset`
- `scan_info`：体扫元数据 `xarray.Dataset`
- `extended_fields`：原生长距离 sidecar
- `product`：产品结果集
- `sitename`：站点名
- `nsweeps`：层数
- `nrays`：射线数
- `effective_earth_radius`：当前体扫使用的有效地球半径

### sweep 数据模型

`fields[sweep]` 往往是最常用的入口。

常见坐标和变量：

- 坐标：`time`、`range`、`azimuth`、`elevation`、`x`、`y`、`z`、`lon`、`lat`
- 场：`dBZ`、`V`、`W`、`ZDR`、`CC`、`PhiDP`、`KDP`，以及后处理生成的场

### 查看体扫信息

```python
radar.summary()
radar.available_fields(sweep=None, range_mode="aligned")
radar.sweep_summary()
radar.get_sweep_field(sweep, field_name, range_mode="aligned", sort_by_azimuth=False)
radar.get_native_sweep_field(sweep, field_name)
radar.has_extended_field(sweep, field_name)
radar.ordered_az(inplace=False)
```

推荐理解方式：

- `summary()`：第一次看一个体扫，先用它
- `available_fields()`：快速列出场
- `sweep_summary()`：按层看摘要
- `get_sweep_field()`：稳定地取一个场
- `ordered_az()`：拿到按方位排序后的视图

#### `summary()`

返回一个轻量 `dict`，通常包含：

- 站点信息
- 扫描类型
- 层数和射线数
- 可用场
- 每层摘要

#### `available_fields(sweep=None, range_mode="aligned")`

返回：

- `sweep is None` 时，返回整个体扫所有可见场名
- 指定 `sweep` 时，返回该层可见场名

如果 `range_mode="native"`，会把原生 sidecar 中的场也纳入结果。

#### `sweep_summary()`

返回一个列表，每一项是一个 `dict`，通常包含：

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

- 只想取一个场
- 想显式指定 `aligned` 或 `native`
- 想顺手按方位角排序

#### `get_native_sweep_field(sweep, field_name)`

如果存在原生长距离 sidecar，就返回原生距离库上的场；
如果没有，则自动回退到对齐场。

#### `has_extended_field(sweep, field_name)`

返回 `True` / `False`。
适合在代码里显式区分 aligned/native 路径。

#### `ordered_az(inplace=False)`

返回或应用“按 azimuth 排序后的视图”。

行为：

- `inplace=False`：返回一个排序视图
- `inplace=True`：直接修改当前对象

### 反射率访问：aligned 与 native

这是整个项目里最重要的读取规则。

对于部分低层：

- `radar.fields[sweep]["dBZ"]` 是历史兼容的对齐版
- `radar.get_native_sweep_field(sweep, "dBZ")` 是原生长距离反射率

建议：

- 做历史兼容绘图、产品和旧流程时，用 aligned
- 做低层完整反射率分析时，用 native

### 产品计算

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

结果写入 `radar.product`。

命名规则：

- aligned 产品保持原有名字
- native 产品通常带 `_native` 后缀

### 在 `PRD` 上做风场反演

```python
radar.retrieve_vad(sweeps=None, field_name=None, range_mode="aligned", **kwargs)
radar.retrieve_vvp(sweep, field_name=None, range_mode="aligned", **kwargs)
radar.retrieve_vwp(sweeps=None, field_name=None, range_mode="aligned", **kwargs)
radar.add_product_VWP(sweeps=None, field_name=None, range_mode="aligned", **kwargs)
```

适合不离开 `PRD` 工作流，直接做单雷达风诊断的场景。

行为：

- `retrieve_vad(...)` 返回按 sweep 和 gate 组织的环状风场反演
- `retrieve_vvp(...)` 返回单层 sweep 上的局地水平风分析
- `retrieve_vwp(...)` 把多层 VAD 结果聚合成垂直风廓线
- `add_product_VWP(...)` 把廓线写入 `radar.product`，变量名为 `VWP_*`

### 剖面与 RHI 提取

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

使用建议：

- `extract_section`：笛卡尔起点终点
- `extract_section_lonlat`：经纬度起点终点
- `extract_rhi`：标准化 RHI dataset
- `get_RHI_data`：旧式低层数组接口
- `get_vcs_data`：旧式低层垂直剖面接口

新代码更推荐：

- `extract_section`
- `extract_section_lonlat`
- `extract_rhi`

### QC、水凝物分类与导出

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

对选定层应用双偏振 QC，返回：

- `inplace=False` 时返回拷贝
- `inplace=True` 时返回原对象

常见新增场：

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

对选定层做 gate 级水凝物分类，返回：

- `inplace=False` 时返回拷贝
- `inplace=True` 时返回原对象

重要参数：

- `profile_height` / `profile_temperature`：可选的 1-D 环境温度廓线
- `temperature`：可选的直接 gate 温度场
- `confidence_field`：可选的最大隶属度置信度输出
- `temperature_field`：可选的“实际用于分类的温度场”输出

行为：

- 把 `HCL` 写回 sweep 场
- 如果 sweep 里已经有 `Zc`、`ZDRc`、`KDPc`，会优先使用订正场
- 不提供温度时，会自动退化到 reduced-variable 无廓线模式

#### `add_hydrometeor_classification(...)`

`classify_hydrometeors(...)` 的原地便捷封装。

#### `to_pyart_radar(...)`

当你要接 Py-ART 工作流时，用这个接口。

重要参数：

- `range_mode`：`aligned` 或 `native`
- `field_names`：只导出指定场
- `use_external`：优先用上游 Py-ART 或强制仓库内兼容对象
- `strict=True`：可选依赖缺失时明确报错

#### `to_xradar_sweeps(...)`

返回 sweep 级 `xarray.Dataset`，适合 xradar 风格工作流。

#### `to_xradar(...)`

在可选依赖满足时，返回 DataTree 风格对象。

#### `to_wsr98d(...)`

写出 WSR98D 基数据文件，并可被 `pycwr.io.read_WSR98D(...)` 回读。

#### `to_nexrad_level2_msg31(...)`

写出可被 Py-ART 读取的 NEXRAD Level II MSG31 文件。
当前支持 `dBZ`、`V`、`W`、`ZDR`、`CC`、`PhiDP`。

#### `to_nexrad_level2_msg1(...)`

写出可被 Py-ART 读取的 NEXRAD Level II MSG1 文件。
这个旧格式当前只支持 `dBZ`、`V`、`W`。

## `pycwr.draw`

这个模块提供对新手更友好的高层绘图接口。

推荐导入：

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

### 返回值

所有高层绘图接口都返回：

```python
EasyPlotResult(fig=..., ax=..., artist=...)
```

这样你可以：

- 继续操作 `matplotlib` Figure / Axes
- 追加注释
- 自己控制保存逻辑

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

最推荐的单层 quicklook 接口。

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

想把同一层画到地图上时，用这个接口。

需要全功能绘图依赖。

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

使用笛卡尔坐标定义剖面端点。
`point_units` 允许 `km` 或 `m`。

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

使用经纬度定义剖面端点。

### `plot_rhi`

```python
plot_rhi(
    radar,
    azimuth,
    field="dBZ",
    ...
)
```

按指定方位角提取并绘制 RHI 风格剖面。

### `plot`

```python
plot(radar, kind="ppi", field="dBZ", sweep=0, show=False, save=None, **kwargs)
```

支持的 `kind`：

- `ppi`
- `ppi_map`
- `section`
- `section_lonlat`
- `rhi`
- `wind_profile`
- `vvp`

### `EasyRadarPlotter`

适合喜欢这种风格的用户：

```python
plotter = EasyRadarPlotter(radar)
plotter.ppi(...)
plotter.section(...)
plotter.rhi(...)
plotter.wind_profile(...)
plotter.vvp(...)
```

## `pycwr.qc`

这个模块提供双偏振质控能力。

### 高层工作流

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

适合已经有 `PRD` 对象时直接调用。

要求：

- 必须有 `dBZ`
- 至少有 `KDP` 或 `PhiDP` 之一

行为：

- 计算订正后的反射率和相关 QC 诊断量
- 把新场写回所选层
- 默认返回新对象，除非 `inplace=True`

### 底层数组工作流

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

适合你已经在 `numpy` 数组层工作，不想先封成 `PRD` 的情况。

返回一个 `dict`，常见键包括：

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

主要参考文献：

- NOAA WDTD 双偏振 QC： https://vlab.noaa.gov/web/wdtd/-/dual-pol-quality-control
- Tang et al. 2020 MRMS QC 更新：doi:10.1175/JTECH-D-19-0165.1
- 吴翀;双偏振雷达的资料质量分析,相态识別及组网应用[D];南京信息工程大学;2018年。在本地业务背景下把晴空回波作为非气象回波类别处理。

### 常用底层 QC 函数

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

当你想自定义 QC 流程而不是直接用封装工作流时，使用这些函数。

## `pycwr.retrieve`

这个模块提供水凝物分类和单雷达风场反演接口。

### 在 `PRD` 上做水凝物分类

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

适合已经有 `PRD`，想把 `HCL` 直接写回各层 sweep 的场景。

要求：

- 必须有 `dBZ`
- 至少有 `ZDR`、`KDP`、`CC`、`LDR` 之一

行为：

- 把 `HCL` 写成 gate 级 sweep 字段
- 可选写出 confidence 和插值后的温度场
- 使用内置的 10 类 fuzzy HID 参数表
- 不提供温度时，会自动退化到 reduced-variable 无廓线模式

### 数组级水凝物分类

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

适合你已经在 `numpy` 数组层工作，不想先封成 `PRD` 的情况。

返回：

- `[1, 10]` 范围内的类别编号
- `return_scores=True` 时额外返回模糊得分立方体
- `return_confidence=True` 时额外返回胜出类别的置信度

### 把环境廓线插值到 gate 温度场

```python
interpolate_temperature_profile(
    gate_altitude,
    profile_height,
    profile_temperature,
    radar_altitude=0.0,
    height_reference="asl",
)
```

典型用法是把 `PRD` sweep 里的 `dataset["z"]` 作为 `gate_altitude`。

### 水凝物类别编号

`pycwr.retrieve.available_hydrometeor_classes()` 返回：

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

内置分类器参考：

- Dolan et al. (2013), *Journal of Applied Meteorology and Climatology*
- Marzano et al. (2006), *Advances in Geosciences*

### 风场反演 helper

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

适合想直接拿结果 dataset，而不是先写回 `radar.product` 的场景。

返回：

- `retrieve_vad(...)`：按 sweep 和 gate 组织的 VAD 数据集
- `retrieve_vwp(...)`：1-D 垂直风廓线数据集
- `retrieve_vvp(...)`：单层 sweep 上的局地 VVP 数据集

历史兼容 helper 仍然保留：

- `vad` / `VAD`
- `vvp` / `VVP`

## `pycwr.interp`

这个模块提供多雷达组网插值和 NetCDF 输出。

### 文件发现与时间选择

```python
parse_radar_time_from_filename(path)
discover_radar_files(radar_dirs, pattern="*.bin*")
select_radar_files(radar_dirs, target_time, tolerance_minutes=10, pattern="*.bin*")
```

#### `parse_radar_time_from_filename(path)`

从文件名里解析扫描时间。
文件名里没有可识别时间时会抛 `ValueError`。

#### `discover_radar_files(...)`

返回一个字典，键是雷达目录名，值是匹配到的文件列表。

#### `select_radar_files(...)`

返回一个列表，每个元素通常包含：

- `radar_id`
- `path`
- `scan_time`
- `delta_seconds`

### 目标网格构建

```python
build_latlon_grid(lon_min, lon_max, lat_min, lat_max, lon_res_deg, lat_res_deg)
```

返回：

- `lon`
- `lat`
- `grid_lon`
- `grid_lat`

### 端到端 3D 组网工作流

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

当你要生成完整多雷达组网产品时，用这个接口。

最少必须定义：

- `target_time`
- `radar_dirs`
- 经纬度范围和分辨率
- `level_heights`

返回：

- 一个 `xarray.Dataset`

参数可以粗分为几类：

- 数据选择：`target_time`、`radar_dirs`、`time_tolerance_minutes`、`pattern`
- 网格定义：`lon_min`、`lon_max`、`lat_min`、`lat_max`、`lon_res_deg`、`lat_res_deg`、`level_heights`
- 合成控制：`field_names`、`field_range_mode`、`composite_method`、`influence_radius_m`、`max_range_km`
- QC：`use_qc`、`qc_method`、`qc_band`、`qc_use_existing_kdp`、`qc_fallback`
- 输出：`output_products`、`output_path`
- 执行控制：`parallel`、`max_workers`

衍生产品和算法口径：

- `CR`：对单站组合反射率先求值，再在组网侧做 `max` 合成。
- `VIL`：使用 Greene and Clark (1972) 的液态水经验关系 `3.44e-6 * Z^(4/7)`，并按业务常见做法对反射率做 `56 dBZ` 封顶；低反射率截止值是 `pycwr` 的可配置工作流参数。
- `ET`：使用阈值回波顶高思路，寻找最高超过阈值的层，并在阈值穿越处做线性插值。
- `ET_TOPPED`：当最高可用层仍超过 `ET` 阈值时，对应格点记为 `1`，用于区分“真实顶高已解析”和“仅到最高采样层”的情况。
- `product_level_heights`：可用于 `VIL` / `ET`，通常建议比显示层 `level_heights` 更密。

当前 3D 组网参考文献与说明：

- 吴翀;双偏振雷达的资料质量分析,相态识別及组网应用[D];南京信息工程大学;2018年，第 4 章给出了单站格点化、共同覆盖区指数权重法和质量权重法的业务背景。
- `pycwr` 的 `quality_weighted` 是在论文所述“质量权重法”基础上做的工程扩展，额外叠加了距离衰减、波束展宽、仰角代表性、QC 掩膜和可选遮挡权重。
- 如果需要与某一业务系统逐式对齐，请不要把 `pycwr` 当前实现直接当成某一篇论文或某一套国家级组网系统的逐项复现。

建议参考：

- Greene, D. R., and R. A. Clark, 1972, *Monthly Weather Review*, doi:10.1175/1520-0493(1972)100<0548:VILWNA>2.3.CO;2
- Lakshmanan et al., 2013, *Weather and Forecasting*, doi:10.1175/WAF-D-12-00084.1
- Lakshmanan et al., 2006, *Weather and Forecasting*, doi:10.1175/WAF942.1
- NOAA WDTD `VIL` / `ET` 产品页
- WSR-88D ROC ICD 2620003Y

### NetCDF 输出

```python
radar_network_3d_to_netcdf(dataset, output_path)
```

写出数据并返回输出路径字符串。

## `pycwr.GraphicalInterface`

这是轻量本地 Web viewer 的 API。

### `create_app`

```python
create_app(allowed_roots=None, auth_token=None)
```

参数：

- `allowed_roots`：viewer 可浏览、可打开的目录
- `auth_token`：固定 token，适合自动化或测试；不传则自动生成随机 token

返回：

- 配置好的 Flask app

### `launch`

```python
launch(host="127.0.0.1", port=8787, open_browser=True)
```

说明：

- 只允许回环地址
- 非首页 API 请求都要求 token
- 文件访问严格限制在允许目录内

## 推荐学习顺序

1. `pycwr.io.read_auto`
2. `PRD.summary()` 和 `PRD.fields`
3. `PRD.get_sweep_field(...)`
4. `pycwr.draw.plot_ppi`
5. `PRD.add_product_CR_xy` / `PRD.add_product_CAPPI_xy` / `PRD.add_product_VIL_xy` / `PRD.add_product_ET_xy`
6. `PRD.apply_dualpol_qc`
7. `pycwr.interp.run_radar_network_3d`

## 相关文档

- [docs/draw_quickstart.md](draw_quickstart.md)
- [docs/web_viewer_quickstart.md](web_viewer_quickstart.md)
- [README_CN.md](../README_CN.md)
