<div align=center>

![avatar](./pictures/NJIAS.png)

</div>

# 中国业务天气雷达开源库

- [English](README.md)
- [开发人员](CONTRIBUTORS_CN.txt)

文档
--------------
PyCWR库的文档正在readthedocs网站不断更新，具体可以查看[PyCWR文档](https://pycwr.readthedocs.io/en/latest/)。


安装pycwr库
----------
### 可以使用conda和pip来安装pycwr库

如果你没安装 cartopy, 推荐使用conda来安装cartopy:

```
conda install -c conda-forge cartopy
```
然后, 可以使用以下pip命令安装:
```
pip install pycwr
```

### 当然你也可以通过源码安装:

```
git clone https://github.com/YvZheng/pycwr.git
cd pycwr
python setup.py install    
```

读取雷达数据 PRD类或者Py-ART的Radar类
----------
```
from pycwr.io import read_auto
file = r"./Z_RADR_I_Z9898_20190828192401_O_DOR_SAD_CAP_FMT.bin.bz2"
PRD = read_auto(file)
PyartRadar = PRD.ToPyartRadar()
```
PRD类的数据结构如下:

<img src="./pictures/PRD_class.png" width="50%" height="50%"/>

可视化PPI图像并叠加地图
----------
```
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pycwr.draw.RadarPlot import Graph, GraphMap
ax = plt.axes(projection=ccrs.PlateCarree())
graph = GraphMap(PRD, ccrs.PlateCarree())
graph.plot_ppi_map(ax, 0, "dBZ", cmap="pyart_NWSRef")
ax.set_title("example of PPI with map", fontsize=16)
plt.show()
```
![avatar](pictures/graph_map.png)

可视化PPI图像
----------
```
fig, ax = plt.subplots()
graph = Graph(PRD)
graph.plot_ppi(ax, 0, "dBZ", cmap="pyart_NWSRef")
graph.add_rings(ax, [0, 50, 100, 150, 200, 250, 300])
ax.set_title("example of PPI", fontsize=16)
ax.set_xlabel("Distance From Radar In East (km)", fontsize=14)
ax.set_ylabel("Distance From Radar In North (km)", fontsize=14)
```
![avatar](pictures/graph.png)


根据起始点经纬度画垂直剖面
----------
```
fig, ax = plt.subplots()
graph = GraphMap(PRD, ccrs.PlateCarree())
graph.plot_vcs_map(ax, (120.8, 27.8), (122.9, 26.8), "dBZ", cmap="pyart_NWSRef")
ax.set_ylim([0,15])
ax.set_ylabel("Height (km)", fontsize=14)
ax.set_xlabel("Latitude, Longitude", fontsize=14)
ax.set_title("VCS exmaple", fontsize=16)
plt.show()
```

![avatar](pictures/vcs_map.png)

根据起始点笛卡尔坐标画垂直剖面
----------
```
fig, ax = plt.subplots()
graph = Graph(PRD)
graph.plot_vcs(ax, (0,0), (150, 0), "dBZ", cmap="pyart_NWSRef")
ax.set_ylim([0,15])
ax.set_ylabel("Height (km)", fontsize=14)
ax.set_xlabel("Distance From Section Start (Uints:km)", fontsize=14)
ax.set_title("VCS exmaple", fontsize=16)
plt.show()
```

![avatar](pictures/vcs.png)

启动图形化界面
----------

```
 python scripts/LaunchGUI.py
```

主窗口如下图所示:

<img src="./pictures/pycwr.png" width="50%" height="50%"/>

更多个例参见[pycwr例子](./notebooks/pycwr_example.ipynb)


项目开发计划
----------

- [x] WSR98D、CINRAD/SA/SB/CB、CINRAD/CC/CCJ、CINRAD/SC/CD 支持
- [ ] Cfradial 读取支持
- [x] 写入 Cfradial 支持
- [x] 自动识别雷达并获取经纬度信息（SA/SB/CB）
- [x] 自动识别雷达数据格式类型
- [x] 转换为 Pyart 雷达对象
- [x] 图形界面支持
- [x] 雷达垂直剖面支持
- [x] 插值算法支持
- [x] PPI 绘制支持，叠加地图支持
- [ ] RHI 绘制支持
- [ ] 多雷达反演算法支持
- [x] 雷达产品算法支持
- [ ] 多普勒雷达/双极化雷达质量控制算法
- [ ] 双极化雷达的 DSD 算法支持
- [ ] 多普勒雷达风场反演支持
- [ ] 雷达定量降水估计算法支持
- [ ] 雷达外推算法支持
- [ ] 雷达定量降水预报算法支持