# pycwr

`pycwr` 是一个面向中国天气雷达业务流程的 Python 工具库，覆盖雷达基数据读取、
几何计算、绘图、质量控制、水凝物分类、单雷达风场反演、多雷达组网插值和导出。

- 当前版本：`1.0.8`
- [English](README.md)
- [接口参考](docs/api_reference_cn.md)
- [雷达组网快速上手](docs/radar_network_quickstart.md)
- [Read the Docs 在线文档](https://pycwr.readthedocs.io/en/latest/)
- [测试索引](test/README.md)
- [绘图快速上手](docs/draw_quickstart.md)

## 为什么是 1.0.8

`1.0.8` 延续了第一条稳定发布线，目标仍然是“可发布、可集成、可维护”。
这次版本在保持 `PA` reader 行为不变的前提下，重点整理了相控阵解码代码结构，并保留默认 colormap 兼容性更新。

重点变化：

- 常见中国天气雷达格式 reader 更稳定，兼容性更清楚
- `PRD` 数据模型更适合做后续处理和二次开发
- 低层反射率同时保留 aligned 和 native 两条访问路径
- 默认依赖更轻，绘图、QC、web viewer、互操作放到 full 依赖
- 几何链路和回归检查更完整
- 单雷达风场反演 `VAD`、`VVP`、`VWP` 和规则格点三维水平风场都已纳入公开接口
- 组网合成和参考风格组合反射率流程更清楚
- `read_auto` 和 `read_PA` 现在能正确识别并构建支持范围内的 `PA` 样本
- `plot_ppi`、`plot_ppi_map` 以及兼容旧代码的 `Graph` / `GraphMap`
  现在都会自动注册 `CN_ref`、`CN_vel` 这类 `pycwr` colormap 名称
- `PA` reader 内部现在把 moment 解码、field padding、距离库构造拆成了更清楚的 helper，后续维护更稳

## 安装

从 PyPI 直接安装：

```bash
python -m pip install pycwr
```

如果需要完整功能依赖：

```bash
python -m pip install "pycwr[full]"
```

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

- `pycwr 1.0.8` 要求 Python `>=3.9`
- 对普通用户来说，优先推荐直接使用 `python -m pip install pycwr`
- 基础安装足够支持 reader、`PRD`、几何、插值和 NetCDF 风格导出
- 全功能安装建议用于绘图、地图绘图、QC、Py-ART/xradar 互操作和 web viewer
- 上游 `arm_pyart` 和 `xradar` 当前要求 Python `>=3.10`，因此在 Python
  `3.9` 上，全功能安装仍可覆盖绘图、QC 和 web viewer，但不包含这两类
  可选互操作依赖
- `1.0.8` 中 `pandas` 已限制为 `<3`，优先保证发布稳定性
- 如果你在本地开发、调试或需要重编译 Cython 扩展，再使用源码安装方式

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

画一张带地图底图的 PPI：

```python
from pycwr.draw import plot_ppi_map

plot_ppi_map(radar, field="dBZ", sweep=0, show=True)
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
| `pycwr.retrieve` | 水凝物分类和风场反演 | `classify_hydrometeors`, `retrieve_vad`, `retrieve_vvp`, `retrieve_vwp`, `retrieve_wind_volume_xy`, `retrieve_wind_volume_lonlat` |
| `pycwr.interp` | 多雷达组网插值 | `run_radar_network_3d` |
| `pycwr.GraphicalInterface` | 本地 Web viewer | `create_app`, `launch` |

## 推荐上手顺序

如果你是第一次接触 `pycwr`，建议按下面顺序理解：

1. 先用 `read_auto()` 读一个真实雷达文件，确认 `PRD` 长什么样。
2. 再用 `summary()`、`available_fields()`、`sweep_summary()` 看清楚每层有哪些参量、多少径向、多少 gate。
3. 然后决定你要的是：
   - 看图：走 `pycwr.draw`
   - 做产品：走 `PRD.add_product_*`
   - 做订正和分类：走 `pycwr.qc` 与 `pycwr.retrieve`
   - 做风场：走 `retrieve_vad / retrieve_vvp / retrieve_vwp`
   - 做组网：走 `run_radar_network_3d`
4. 真正落到业务前，再确认 `aligned/native` 距离库、坐标单位、网格分辨率和输出格式。

## 核心对象模型

所有 reader 最终都返回 `pycwr.core.NRadar.PRD`。

`PRD` 最重要的组成有：

- `fields`：每层一个 `xarray.Dataset`
- `scan_info`：站点和扫描元数据
- `extended_fields`：aligned/native 距离不一致时的原生 sidecar
- `product`：产品结果集

可以把它简单理解成：

- `scan_info` 负责“这一整部体扫是什么”，比如站点经纬高、扫到几层、各层仰角、Nyquist 速度、无模糊距离等
- `fields[i]` 负责“第 `i` 层真正有哪些矩阵场”，比如 `dBZ`、`V`、`W`、`ZDR`、`CC` 等
- `extended_fields` 负责“当某些字段原生距离更长，但历史共享距离库更短时，原生版本放在哪里”
- `product` 负责“基于这部体扫进一步算出来的网格产品”

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

这个差异最常见在低仰角层：

- 速度和谱宽等多普勒变量因为业务体制、PRT 或编码方式，距离可能比较短
- 反射率实际还能探测得更远
- 老流程为了让不同变量共享一套 gate，常常把反射率一起裁短
- `pycwr` 现在保留了两套访问口径，方便你既做历史兼容，也做真实反射率分析

建议：

- 需要和旧业务链严格一致，用 `range_mode="aligned"`
- 需要低层完整反射率探测范围，用 `range_mode="native"`
- 如果你是在核对“最大探测距离”或做组合反射率范围分析，不要只看 aligned，要先确认 native

快速示例：

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

#### 默认 colormap 策略

当你不传 `cmap`，也就是保持 `cmap=None` 时，`pycwr.draw` 会根据字段自动选择默认色标。

- `dBZ` 默认使用反射率色标，例如 `CN_ref`
- `V` 默认使用速度色标，例如 `CN_vel`
- `HCL` 默认使用离散的水凝物分类色标
- 未知字段默认回退到 `viridis`

这套默认策略在 `plot_ppi`、`plot_ppi_map`、`plot_rhi`、`plot_section`，
以及兼容旧代码的 `Graph` / `GraphMap` 类里都是一致的。

推荐写法：

```python
from pycwr.draw import plot_ppi, plot_ppi_map

plot_ppi(radar, field="dBZ", sweep=0, show=True)
plot_ppi_map(radar, field="dBZ", sweep=0, show=True)
```

如果你想显式指定色标：

```python
plot_ppi(radar, field="dBZ", sweep=0, cmap="CN_ref", cmap_bins=16, show=True)
plot_ppi_map(radar, field="V", sweep=0, cmap="CN_vel", min_max=(-27.0, 27.0), show=True)
```

旧类接口也可以继续这样写：

```python
from pycwr.draw.RadarPlot import GraphMap
from cartopy import crs as ccrs

display = GraphMap(radar, ccrs.PlateCarree())
display.plot_ppi_map(ax, 0, "dBZ", cmap="CN_ref")
```

说明：

- `CN_ref`、`CN_vel` 这类 `pycwr` colormap 名称现在会在绘图时自动注册，不需要再先“预热”一次封装绘图函数。
- 如果当前雷达文件没有某个字段，绘图会明确报字段不存在，而不是先报 colormap 错误。比如有些文件只有 `dBZ`，没有 `V`。
- 如果你要自己做更底层的 matplotlib 定制，也可以显式 `import pycwr.draw.colormap` 来注册这些色标名。

一般建议：

- 先用 `plot_ppi`、`plot_rhi`、`plot_section` 快速看数值和结构
- 需要地图底图时再用 `plot_ppi_map`
- 需要风场展示时用 `plot_vvp` 和 `plot_wind_profile`
- 做正式业务图时，优先先确认数值和网格对不对，再去微调 colormap、extent、annotation

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

快速示例：

```python
radar.add_product_CAPPI_xy(x, y, 3000.0)
radar.add_product_VIL_xy(x, y, [1000.0, 2000.0, 3000.0])
```

产品接口的理解方式可以统一成一句话：

- 先准备目标网格
- 再指定产品类型
- 必要时指定高度层、体积层或输出字段
- 最终结果写回 `radar.product`

实际使用时建议统一几件事：

- `x/y/z` 内部都按米处理
- 经纬度网格和投影网格不要混用
- 同一个流程里，网格分辨率、范围和缺测值口径尽量固定
- 做跨时次对比时，不要一张图 1 km、一张图 2 km 地混着算

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

`pycwr` 现在内置四条单雷达风诊断工作流：

- `retrieve_vad`：单层或多层的环状谐波拟合
- `retrieve_vvp`：单层局地最小二乘水平风反演
- `retrieve_vwp`：由多层 VAD 聚合出来的垂直风廓线
- `retrieve_wind_volume_xy` / `retrieve_wind_volume_lonlat`：把多个 sweep 的 VVP 结果重建到固定高度层后的规则格点三维水平风场

示例：

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

如果要把风廓线写回产品集：

```python
radar.add_product_VWP(sweeps=[0, 1, 2], max_range_km=40.0, height_step=500.0)
```

如果要把规则格点三维水平风场写回产品集：

```python
radar.add_product_WIND_VOLUME_xy(
    XRange=np.arange(-20_000.0, 20_001.0, 10_000.0),
    YRange=np.arange(-20_000.0, 20_001.0, 10_000.0),
    level_heights=np.array([500.0, 1000.0, 1500.0]),
    sweeps=[0, 1, 2],
)
```

#### 三维风场反演到底是什么

`retrieve_wind_volume_xy` 和 `retrieve_wind_volume_lonlat` 返回的是单雷达、
固定高度层、规则格点的水平风体。

应该这样理解：

- 这是规则格点上的 `u/v` 水平风场
- 每个水平格点对应一根风场廓线
- 它适合做业务诊断、产品分析和固定高度层风场展示

不应该这样理解：

- 它不是完整的 `u/v/w` 三维动力风场
- 它不是双多普勒或多雷达变分三维风场
- 它也不等价于对流动力研究里那种垂直速度反演结果

#### 三维风场反演的计算流程

这条链路的主要步骤是：

1. 对每个候选 sweep 选择速度场，优先使用订正速度。
2. 在每个 sweep 上做局地 `VVP` 反演，得到单层水平风诊断结果。
3. 如果 `sweeps=None` 或 `sweeps="auto"`，先对各层结果做质量汇总，再自动挑选更有代表性的 sweep。
4. 把各层有效 `VVP` 结果通过反距离加权映射到目标水平格点。
5. 再把这些 sweep 结果重建到固定高度层，优先做局地垂直聚合，必要时做有限间隔的线性插值。
6. 同时输出残差、样本数、支持层数、质量评分等诊断量。

这样做的目的不是让结果“看起来更平滑”，而是为了让真实业务数据在缺测多、
距离不一致、低质量 sweep 混入时仍然尽量稳。

#### 三维风场反演推荐接口

直接使用：

- `radar.retrieve_wind_volume_xy(...)`：输出 `x/y/z` 规则格点三维水平风场
- `radar.retrieve_wind_volume_lonlat(...)`：输出规则经纬网格三维水平风场
- `radar.add_product_WIND_VOLUME_xy(...)`
- `radar.add_product_WIND_VOLUME_lonlat(...)`

推荐的业务化使用方式：

- `sweeps=None`：自动选层，优先推荐
- `sweeps="all"`：把所有 sweep 都纳入，不做自动筛选
- `sweeps=[0, 1, 2]`：显式指定层号，只在你非常明确要哪些层时使用
- `max_range_km`：一定要给一个合理上限，不要把远距离低质量速度盲目并进去
- `workers`：在多核机器上可显著加速

一个更完整的业务示例：

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

#### 三维风场反演关键参数怎么理解

- `XRange`、`YRange`、`level_heights`：输出目标格点，单位都是米
- `XLon`、`YLat`：经纬网格模式下的目标网格，单位是度
- `sweeps`：`None` 或 `"auto"` 表示自动选层，`"all"` 表示全部层都参与，列表表示手工指定层
- `max_range_km`：最重要的稳定性控制参数之一
- `az_num`、`bin_num`：局地 `VVP` 窗口大小
- `azimuth_step`、`range_step`：单层 `VVP` 输出抽样步长
- `horizontal_radius_m`：把 sweep 结果映射到目标格点时的水平搜索半径
- `max_horizontal_radius_m`：样本过稀时允许扩张的最大搜索半径
- `vertical_tolerance_m`：固定高度层聚合时的局地垂直容差
- `max_vertical_gap_m`：允许做垂直插值的最大层间距
- `workers`：基于进程的并行计算参数

更实用的调参建议：

- 速度缺测多、噪声大时，`az_num` 宁可大一点，不要过小
- `bin_num` 不宜太大，否则容易过度平滑
- 如果结果断得很厉害，先检查 `max_range_km` 和搜索半径，不要先怀疑公式错了
- 如果某些 sweep 明显质量差，优先用 `sweeps=None` 让它自动筛掉

#### 三维风场反演返回什么

返回的数据集里，核心变量包括：

- `u`、`v`：东西向、南北向水平风分量
- `wind_speed`、`wind_direction`
- `fit_rmse`：局地拟合残差
- `valid_count`：有效样本数
- `source_sweep_count`：参与该体素重建的 sweep 数
- `neighbor_count`：水平邻域支持数
- `effective_vertical_support`：垂直方向上实际参与支持的层数
- `sweep_spread_m`：参与支持层的垂直离散程度
- `quality_score`：`0-100` 的启发式质量评分
- `quality_flag`：`0/1/2/3`，分别表示无效、低可信、可用、高可信

属性里最值得看的是：

- `selection_mode`
- `requested_sweeps`
- `selected_sweeps`
- `rejected_sweeps`
- `rejected_sweep_reasons`
- `workers_used`
- `parallel_mode`

真实业务里，这些元数据往往比图像本身更能解释“为什么这里有空洞”。

#### 三维风场结果怎么判断质量

如果真实个例里出现大片缺测，不要立刻判定算法失效。

最常见的原因有：

- 无降水或弱回波，本来就没有可靠的径向速度
- 方位覆盖不足，局地拟合不稳
- 请求的高度层超出了当前 sweep 的有效几何支撑范围
- 搜索半径太小，或者距离上限设得过于激进

排查时优先看：

- `quality_score`
- `quality_flag`
- `fit_rmse`
- `valid_count`
- `source_sweep_count`
- `selected_sweeps` 和 `rejected_sweep_reasons`

#### 风场这几条产品线怎么区分

- `VAD` 更适合做分层平均风和风廓线基础层
- `VVP` 更适合看某一层的局地水平风结构
- `VWP` 是把多个 VAD 层做稳健聚合后得到的业务化风廓线
- `WIND_VOLUME` 是把多个 sweep 的 `VVP` 结果重建到规则高度层之后得到的单雷达三维水平风场

#### 三维风场反演与并行

`workers` 使用的是基于进程的并行，不是 Python 线程池。

经验上：

- 小网格时提速可能不明显
- 中等及以上规模网格，`workers=2` 或 `workers=4` 往往就有明显收益
- 并行数不是越大越好，进程调度和数据搬运开销会限制继续扩展

做正式验收时，建议至少对一个代表性个例同时跑一遍：

- `workers=1`
- `workers=2` 或 `workers=4`

然后确认：

- `u/v` 数值一致
- `quality_flag` 一致
- 总耗时确实下降

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

如果你要做单层组合反射率或单层 CAPPI，建议先固定四个要素：

- 经度范围
- 纬度范围
- 网格分辨率
- 高度层

这样后面跨时次、跨雷达对比才有意义。组网结果是否“看起来一致”，不要只看图，要同时核对：

- 参与站点列表
- 每站该时次哪些 sweep 和字段参与了计算
- 反射率用的是 aligned 还是 native 距离库
- 网格边界和分辨率是否完全一致

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
- [雷达组网快速上手](docs/radar_network_quickstart.md)：面向小白的组网步骤和示例
- [测试索引](test/README.md)：按功能整理的可运行示例
- [绘图快速上手](docs/draw_quickstart.md)：绘图入口和示例
- [Read the Docs 在线文档](https://pycwr.readthedocs.io/en/latest/)：持续更新的完整手册
- Web viewer：直接运行 `python scripts/LaunchGUI.py`

## 给发布用户的说明

`1.0.8` 最需要明确的行为规则有这几条：

- 所有 reader 统一返回稳定的 `PRD` 对象
- 低层反射率可以显式选择 aligned 或 native 距离库
- QC 和水凝物分类可以把订正场直接写回 `PRD`
- 单雷达风场反演已经属于公开接口的一部分，包含规则格点三维水平风场
- 多雷达组网已经有明确的高层工作流入口
- 依赖拆分成 base 和 full 两组，安装摩擦更低
- 符合支持标记的 `PA` 样本会走 `read_PA`，不再被误判成 `WSR98D`
- 默认绘图现在会在 easy API 和旧类接口里统一解析 `pycwr` colormap 名称

如果你是准备把 `pycwr` 接进自己的业务脚本，最实用的建议只有三条：

1. 先锁定你要用的输入格式和样本基线，不要一边接入一边猜格式。
2. 先锁定数值口径，再做绘图和界面，不要只凭图像判断算法对不对。
3. 对外发布前，把你实际会用到的 reader、产品、绘图和导出路径各跑一遍真实样本回归。
