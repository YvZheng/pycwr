# pycwr

面向中国天气雷达业务数据的 Python 工具库，覆盖读取、分析、绘图、质控、组网插值和导出。

- 当前版本：`1.0.0`
- [English](README.md)
- [接口参考](docs/api_reference_cn.md)
- [测试索引](test/README.md)
- [绘图快速上手](docs/draw_quickstart.md)
- [Web viewer 快速上手](docs/web_viewer_quickstart.md)
- [开发人员](CONTRIBUTORS_CN.txt)

## 1.0.0 重磅更新

这是一次以兼容性、几何正确性和依赖减负为核心的大版本整理。

模块重点：

- `pycwr.io`：默认安装更轻，但核心雷达 reader 保持可用
- `pycwr.core`：雷达坐标、笛卡尔坐标、地理坐标主链统一到更严格的几何实现，并补上闭环校验
- `pycwr.draw` 与 `GraphicalInterface`：可选依赖从基础导入链拆开，默认环境更稳
- `pycwr.qc`：双偏振和 X 波段流程补成更适合参考的示例测试
- `pycwr.retrieve`：补齐了带温度廓线和无廓线两种水凝物分类接口
- `pycwr.interp`：组网插值入口、样例测试和文档更加完整

性能重点：

- 去掉了重复几何实现，主计算链只保留高精度公式
- 向量化坐标转换减少了中间数组和重复角度网格构造
- 默认依赖更精简，跨平台安装成功率更高，环境准备更简单

## pycwr 能做什么

`pycwr` 的主线很简单：

1. 用 `pycwr.io.read_auto` 读取雷达文件
2. 拿到 `PRD` 对象
3. 对这个体扫做绘图、产品计算、QC、组网或导出

核心能力：

- 读取 `WSR98D`、`SAB`、`CC`、`SC`、`PA`
- 保持历史对齐工作流，同时暴露原生长距离反射率
- 绘制 PPI、地图 PPI、RHI、垂直剖面
- 计算 `CR`、`CAPPI` 和 3D 组网产品
- 运行双偏振质控
- 运行水凝物分类（`HCL`），支持带廓线和无廓线两种模式
- 导出到 Py-ART、xradar、NetCDF
- 启动本地 Web viewer

## 安装

安装默认依赖和基础包：

```bash
python -m pip install -r requirements-core.txt
python -m pip install .
```

安装全功能版本：

```bash
python -m pip install -r requirements-full.txt
python -m pip install ".[full]"
```

修改 `pycwr/core/RadarGridC.pyx` 后重编译：

```bash
python setup.py build_ext --inplace
```

## 5 分钟上手

读取一个雷达文件：

```python
from pycwr.io import read_auto

radar = read_auto("Z_RADR_I_Z9046_20260317065928_O_DOR_SAD_CAP_FMT.bin.bz2")
print(radar.summary())
```

查看这个体扫：

```python
print(radar.available_fields())
print(radar.sweep_summary()[0])
print(radar.fields[0]["dBZ"])
```

画一张最简单的图：

```python
from pycwr.draw import plot_ppi

plot_ppi(radar, field="dBZ", sweep=0, show=True)
```

提取一个垂直剖面：

```python
from pycwr.draw import plot_section

plot_section(radar, start=(-50, 0), end=(50, 0), field="dBZ", show=True)
```

## 核心概念

### `PRD`

所有 reader 最终都返回 `pycwr.core.NRadar.PRD`。

最重要的成员：

- `fields`：每层一个 `xarray.Dataset`
- `scan_info`：站点和扫描元数据
- `extended_fields`：原生长距离场的 sidecar
- `product`：产品结果，如 `CR`、`CAPPI`、3D 网格

对新手最有用的方法：

- `summary()`：整个体扫的紧凑摘要
- `available_fields(sweep=None)`：查看所有层或某一层有哪些场
- `sweep_summary()`：每层的摘要信息
- `get_sweep_field(sweep, field_name, range_mode="aligned" | "native")`
- `ordered_az(inplace=False)`

### 对齐反射率与原生反射率

对部分低层，`pycwr` 会同时保留两种 `dBZ` 视图：

- `radar.fields[sweep]["dBZ"]`：历史兼容的对齐结果
- `radar.get_native_sweep_field(sweep, "dBZ")`：原生长距离反射率

什么时候用哪个：

- 要与旧流程保持一致，用对齐版
- 要看低层完整反射率探测范围，用原生版

## 公开接口

### 1. 读取

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

统一签名：

```python
read_xxx(
    filename,
    station_lon=None,
    station_lat=None,
    station_alt=None,
    effective_earth_radius=None,
)
```

### 2. 绘图

推荐优先使用：

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

这些函数返回 `EasyPlotResult`，里面有 `fig`、`ax`、`artist`。

### 3. 产品计算

常用 `PRD` 方法：

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

### 4. 质控

```python
from pycwr.qc import apply_dualpol_qc, run_dualpol_qc
```

也可以直接在 `PRD` 上调用：

```python
processed = radar.apply_dualpol_qc(inplace=False)
```

### 5. 组网 3D 工作流

```python
from pycwr.interp import run_radar_network_3d
```

用于构建规则经纬度网格上的 3D 雷达组网产品，并可直接写 NetCDF。

### 6. 水凝物分类

```python
from pycwr.retrieve import classify_hydrometeors, interpolate_temperature_profile
```

既可以直接对数组做分类，也可以把 `HCL` 直接写回 `PRD`。

### 7. 风场反演

```python
from pycwr.retrieve import retrieve_vad, retrieve_vvp, retrieve_vwp
```

可以做：

- 一个或多个 sweep 的 VAD 环状反演
- 单层 sweep 的局地 VVP 反演
- 直接生成或写回 `PRD` 的 VWP 垂直风廓线

### 8. Web viewer

用 Python 入口启动：

```python
from pycwr.GraphicalInterface import create_app, launch
```

或者使用脚本：

```bash
python scripts/LaunchGUI.py
```

当前 viewer 的边界：

- 只允许监听本机回环地址
- 所有 API 请求都要求 token
- 文件访问限制在允许目录内

## 常见工作流

### 画一层 PPI

```python
from pycwr.io import read_auto
from pycwr.draw import plot_ppi

radar = read_auto("your_radar_file")
plot_ppi(radar, field="dBZ", sweep=0, show=True)
```

### 计算 CR、CAPPI、VIL 和 ET

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

### 做 QC 并使用订正场

```python
qc_radar = radar.apply_dualpol_qc(inplace=False, clear_air_mode="mask")
plot_ppi(qc_radar, field="Zc", sweep=0, show=True)
```

补充说明：

- `CLEAR_AIR_MASK` 会作为独立场输出，用于标记“可能的晴空 / 生物弱回波”。
- `clear_air_mode="label"` 只做标记，保持原有 QC 接受逻辑。
- `clear_air_mode="mask"` 会把这类弱晴空回波从 `QC_MASK`、`Zc`、`KDPc` 以及后续 QC 产品里剔除。
- 这条晴空回波分支是保守诊断，不等价于完整的相态识别 / HCA 分类器。

### 做水凝物分类

带环境温度廓线：

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

不提供温度廓线：

```python
radar.add_hydrometeor_classification(
    band="C",
    confidence_field="HCL_CONF",
)
```

说明：

- `HCL` 是 gate 级别的 sweep 字段，不是笛卡尔格点产品。
- 分类使用 `dBZ/ZDR/KDP/CC/LDR` 共享的对齐距离范围。
- 如果不提供 `profile_height/profile_temperature`，`pycwr` 会自动退化到不带温度项的 reduced-variable 模式。
- 类别编号固定为 10 类：毛毛雨、雨、冰晶、干雪聚合体、湿雪、竖直冰晶、低密度霰、高密度霰、冰雹、大水滴。

### 做单雷达风场反演

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

说明：

- `retrieve_vad(...)` 返回按 sweep 和 gate 组织的 VAD 结果。
- `retrieve_vvp(...)` 返回单层 sweep 上的局地 VVP 分析。
- `retrieve_vwp(...)` 会把多层 VAD 结果聚合成垂直风廓线。
- `add_product_VWP(...)` 会把结果写入 `radar.product`，变量名为 `VWP_*`。

### 计算 3D 组网

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

补充说明：

- 目前 3D 组网支持 `dBZ` 一类标量场。
- `run_radar_network_3d(...)` 暂不支持速度场 3D 组网。
- `VIL` / `ET` 可以通过 `product_level_heights` 使用比显示层 `level_heights` 更密的垂直层。
- 请求 `ET` 时，输出数据集会同时带出 `ET_TOPPED`，用于标记“最高可用层仍超过阈值”的柱状格点。
- `VIL` 采用 Greene and Clark (1972) 的液态水经验关系，并按业务常见做法对反射率做 `56 dBZ` 封顶；低反射率截止值仍然是可配置的工程参数。

`pycwr.interp.RadarInterp` 当前采用和引用的主要算法文献：

- `VIL`：Greene, D. R., and R. A. Clark, 1972, *Monthly Weather Review*, doi:10.1175/1520-0493(1972)100<0548:VILWNA>2.3.CO;2
- `VIL` 业务说明：NOAA WDTD, https://vlab.noaa.gov/web/wdtd/-/vertically-integrated-liquid-vil-
- `ET / ETOP`：Lakshmanan et al., 2013, *Weather and Forecasting*, doi:10.1175/WAF-D-12-00084.1
- `ET / ETOP` 业务说明：NOAA WDTD, https://vlab.noaa.gov/web/wdtd/-/xx-dbz-echo-top-et-
- `Enhanced Echo Tops` 元数据语义：WSR-88D ROC ICD 2620003Y, https://www.roc.noaa.gov/public-documents/icds/2620003Y.pdf
- 雷达组网重叠区处理：吴翀;双偏振雷达的资料质量分析,相态识別及组网应用[D];南京信息工程大学;2018年，第 4 章；以及 Lakshmanan et al., 2006, *Weather and Forecasting*, doi:10.1175/WAF942.1

`pycwr.qc` 当前采用和引用的主要 QC 文献：

- 双偏振 QC 业务说明：NOAA WDTD, https://vlab.noaa.gov/web/wdtd/-/dual-pol-quality-control
- MRMS QC 更新：Tang et al., 2020, *Journal of Atmospheric and Oceanic Technology*, doi:10.1175/JTECH-D-19-0165.1
- 本地雷达业务背景下的晴空 / 非气象回波处理：吴翀;双偏振雷达的资料质量分析,相态识別及组网应用[D];南京信息工程大学;2018年

### 导出到 Py-ART / xradar

```python
pyart_radar = radar.to_pyart_radar()
sweeps = radar.to_xradar_sweeps()
tree = radar.to_xradar(strict=True)
```

### 导出到 WSR98D / NEXRAD Level II

```python
radar.to_wsr98d("./roundtrip_wsr98d.bin", overwrite=True)
radar.to_nexrad_level2_msg31("./export_msg31.ar2v", overwrite=True)
radar.to_nexrad_level2_msg1("./export_msg1.ar2v", overwrite=True)
```

说明：

- `to_wsr98d(...)` 会写出可被 `pycwr` 自己回读的 WSR98D 基数据文件
- `to_nexrad_level2_msg31(...)` 会写出可被 Py-ART 读取的 NEXRAD Level II MSG31 文件
- `to_nexrad_level2_msg1(...)` 会写出可被 Py-ART 读取的 NEXRAD Level II MSG1 文件
- 当前导出只支持 `ppi` 体扫
- MSG31 支持 `dBZ/V/W/ZDR/CC/PhiDP`
- MSG1 当前只支持 `dBZ/V/W`

## 单位约定

核心内部约定：

- 距离和高度：米
- 内部三角函数计算：弧度
- 经纬度输入输出：度

历史兼容仍然保留：

- `azimuth`、`elevation`、`fixed_angle` 仍然是度
- 部分旧绘图接口仍然沿用千米

## 测试与样例数据

运行自动化测试：

```bash
python3 -m unittest discover -s test -p 'test_*.py'
```

推荐关注的测试：

- [test/README.md](test/README.md)：示例脚本和回归脚本总索引
- `test_examples_public_api.py`：面向真实 sample 的快速上手集成测试
- `test_examples_qc.py`：QC 原语与流程
- `test_examples_hid.py`：水凝物分类示例，覆盖带廓线和无廓线两种调用
- `test_examples_interp.py`：多雷达组网插值示例
- `test_regression_geometry.py`：几何闭环校验
- `test_regression_security.py`：路径边界与恶意输入回归

如果当前 workspace 没有对应样例文件，部分 sample test 会自动跳过。
真实 3D 组网 slow test 默认不跑，需要显式设置 `PYCWR_RUN_SLOW_SAMPLE_TESTS=1`。

## 文档导航

- [docs/api_reference_cn.md](docs/api_reference_cn.md)：接口参考
- [docs/draw_quickstart.md](docs/draw_quickstart.md)：绘图快速上手
- [docs/web_viewer_quickstart.md](docs/web_viewer_quickstart.md)：本地 viewer 使用说明

## 许可证

本仓库采用 **PolyForm Noncommercial 1.0.0**。
默认禁止商业用途、收费分发和盈利部署。如需商业授权，请单独联系。
见 [LICENSE.txt](LICENSE.txt)。
