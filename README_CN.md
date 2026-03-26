# pycwr

`pycwr` 是一个面向中国天气雷达业务流程的 Python 工具库，覆盖雷达基数据读取、
几何计算、绘图、质量控制、水凝物分类、单雷达风场反演、多雷达组网插值和导出。

- 当前版本：`1.0.3`
- [English](README.md)
- [接口参考](docs/api_reference_cn.md)
- [测试索引](test/README.md)
- [绘图快速上手](docs/draw_quickstart.md)
- [Web viewer 快速上手](docs/web_viewer_quickstart.md)
- [开发人员](CONTRIBUTORS_CN.txt)

## 为什么是 1.0.3

`1.0.3` 延续了第一条稳定发布线，目标仍然是“可发布、可集成、可维护”。

重点变化：

- 常见中国天气雷达格式 reader 更稳定，兼容性更清楚
- `PRD` 数据模型更适合做后续处理和二次开发
- 低层反射率同时保留 aligned 和 native 两条访问路径
- 默认依赖更轻，绘图、QC、web viewer、互操作放到 full 依赖
- 几何链路和回归检查更完整
- 单雷达风场反演 `VAD`、`VVP`、`VWP` 已纳入公开接口
- 组网合成和参考风格组合反射率流程更清楚

## 安装

基础安装：

```bash
python -m pip install -r requirements-core.txt
python -m pip install .
```

全功能安装：

```bash
python -m pip install -r requirements-full.txt
python -m pip install ".[full]"
```

说明：

- `pycwr 1.0.3` 要求 Python `>=3.9`
- 基础安装足够支持 reader、`PRD`、几何、插值和 NetCDF 风格导出
- 全功能安装建议用于绘图、地图绘图、QC、Py-ART/xradar 互操作和 web viewer
- 上游 `arm_pyart` 和 `xradar` 当前要求 Python `>=3.10`，因此在 Python
  `3.9` 上，全功能安装仍可覆盖绘图、QC 和 web viewer，但不包含这两类
  可选互操作依赖
- `1.0.3` 中 `pandas` 已限制为 `<3`，优先保证发布稳定性

修改 `pycwr/core/RadarGridC.pyx` 后重编译：

```bash
python setup.py build_ext --inplace
```

构建发布产物：

```bash
python -m build
```

## 5 分钟上手

读取一个雷达文件并查看体扫摘要：

```python
from pycwr.io import read_auto

radar = read_auto("Z_RADR_I_Z9046_20260317065928_O_DOR_SAD_CAP_FMT.bin.bz2")
print(radar.summary())
print(radar.available_fields())
print(radar.sweep_summary()[0])
```

取某一层某一个场：

```python
dBZ0 = radar.get_sweep_field(0, "dBZ")
velocity0 = radar.get_sweep_field(0, "V")
```

画一张 PPI：

```python
from pycwr.draw import plot_ppi

plot_ppi(radar, field="dBZ", sweep=0, show=True)
```

提取垂直剖面：

```python
from pycwr.draw import plot_section

plot_section(radar, start=(-50, 0), end=(50, 0), field="dBZ", show=True)
```

生成一个简单产品：

```python
import numpy as np

x = np.arange(-150_000.0, 150_001.0, 1_000.0)
y = np.arange(-150_000.0, 150_001.0, 1_000.0)
radar.add_product_CR_xy(x, y)
print(radar.product)
```

## API 地图

| 模块 | 用途 | 推荐起点 |
| --- | --- | --- |
| `pycwr.io` | 读取和写出雷达基数据 | `read_auto`, `read_WSR98D`, `read_SAB`, `read_CC`, `read_SC`, `read_PA` |
| `pycwr.core` | 核心体扫对象、几何和导出辅助 | `PRD`, `radar.summary()`, `radar.get_sweep_field()` |
| `pycwr.draw` | 绘图和快速出图 | `plot_ppi`, `plot_ppi_map`, `plot_rhi`, `plot_section`, `plot_vvp`, `plot_wind_profile` |
| `pycwr.qc` | 双偏振质控 | `apply_dualpol_qc`, `run_dualpol_qc` |
| `pycwr.retrieve` | 水凝物分类和风场反演 | `classify_hydrometeors`, `retrieve_vad`, `retrieve_vvp`, `retrieve_vwp` |
| `pycwr.interp` | 多雷达组网插值 | `run_radar_network_3d` |
| `pycwr.GraphicalInterface` | 本地 Web viewer | `create_app`, `launch` |

## 核心对象模型

所有 reader 最终都返回 `pycwr.core.NRadar.PRD`。

`PRD` 最重要的组成有：

- `fields`：每层一个 `xarray.Dataset`
- `scan_info`：站点和扫描元数据
- `extended_fields`：aligned/native 距离不一致时的原生 sidecar
- `product`：产品结果集

最常用的查看接口：

- `summary()`：体扫整体摘要
- `available_fields(sweep=None, range_mode="aligned")`
- `sweep_summary()`
- `get_sweep_field(sweep, field_name, range_mode="aligned", sort_by_azimuth=False)`
- `get_native_sweep_field(sweep, field_name)`
- `ordered_az(inplace=False)`

### aligned 与 native 反射率

部分低层反射率可能同时存在两条访问路径：

- aligned：和历史流程一致的共享距离库
- native：反射率原生距离库，没有被速度距离库截短

建议：

- 需要和旧业务链严格一致，用 `range_mode="aligned"`
- 需要低层完整反射率探测范围，用 `range_mode="native"`

示例：

```python
aligned = radar.get_sweep_field(0, "dBZ", range_mode="aligned")
native = radar.get_sweep_field(0, "dBZ", range_mode="native")
```

## 常见工作流

### 绘图

推荐的公开绘图接口：

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

这些函数统一返回 `EasyPlotResult`，里面包含 `fig`、`ax`、`artist`。

### 产品计算

常用 `PRD` 产品接口：

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

示例：

```python
radar.add_product_CAPPI_xy(x, y, 3000.0)
radar.add_product_VIL_xy(x, y, [1000.0, 2000.0, 3000.0])
```

### 质量控制

```python
from pycwr.qc import apply_dualpol_qc

qc_radar = apply_dualpol_qc(radar, inplace=False, clear_air_mode="mask")
```

常见订正场包括 `Zc`、`ZDRc`、`PhiDPc`、`KDPc`，以及在开启时输出的
`QC_MASK`、`CLEAR_AIR_MASK` 等标记场。

### 水凝物分类

```python
hcl_radar = radar.classify_hydrometeors(
    inplace=False,
    band="C",
    profile_height=[0.0, 2000.0, 4000.0, 8000.0],
    profile_temperature=[24.0, 12.0, 2.0, -16.0],
    confidence_field="HCL_CONF",
)
```

如果你已经有数组，也可以直接调用
`pycwr.retrieve.classify_hydrometeors(...)`。

### 风场反演

`pycwr` 现在内置三条单雷达风诊断工作流：

- `retrieve_vad`：单层或多层的环状谐波拟合
- `retrieve_vvp`：单层局地最小二乘水平风反演
- `retrieve_vwp`：由多层 VAD 聚合出来的垂直风廓线

示例：

```python
vad = radar.retrieve_vad(sweeps=[0, 1, 2], max_range_km=40.0, gate_step=4)
vvp = radar.retrieve_vvp(0, max_range_km=20.0, az_num=91, bin_num=5)
vwp = radar.retrieve_vwp(sweeps=[0, 1, 2], max_range_km=40.0, height_step=500.0)
```

如果要把风廓线写回产品集：

```python
radar.add_product_VWP(sweeps=[0, 1, 2], max_range_km=40.0, height_step=500.0)
```

### 导出与互操作

常用导出接口：

- `radar.to_wsr98d(...)`
- `radar.to_nexrad_level2_msg31(...)`
- `radar.to_nexrad_level2_msg1(...)`
- `radar.to_pyart_radar(...)`
- `radar.to_xradar(...)`
- `radar.to_cfgridded_netcdf(...)`

适合把 `pycwr` 数据继续送进 Py-ART、xradar 或其他 NetCDF 流程。

### 多雷达组网

```python
from pycwr.interp import run_radar_network_3d
```

这是规则经纬度网格 3D 组网产品的推荐高层入口，也支持直接写 NetCDF。

### Web viewer

```python
from pycwr.GraphicalInterface import create_app, launch
```

或者直接用脚本：

```bash
python scripts/LaunchGUI.py
```

viewer 设计上只允许本机访问，并要求 token 才能调用 API。

## 下一步读什么

- [接口参考](docs/api_reference_cn.md)：更细的公开 API 说明
- [测试索引](test/README.md)：按功能整理的可运行示例
- [绘图快速上手](docs/draw_quickstart.md)：绘图入口和示例
- [Web viewer 快速上手](docs/web_viewer_quickstart.md)：本地 viewer 使用方式

## 给发布用户的说明

`1.0.3` 最需要明确的行为规则有这几条：

- 所有 reader 统一返回稳定的 `PRD` 对象
- 低层反射率可以显式选择 aligned 或 native 距离库
- QC 和水凝物分类可以把订正场直接写回 `PRD`
- 单雷达风场反演已经属于公开接口的一部分
- 多雷达组网已经有明确的高层工作流入口
- 依赖拆分成 base 和 full 两组，安装摩擦更低
