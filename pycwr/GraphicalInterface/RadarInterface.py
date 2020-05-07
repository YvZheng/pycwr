# -*- coding: utf-8 -*-

"""
Module implementing MainWindow.
"""
import os
from .RadarUI import Ui_MainWindow
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5 import QtWidgets
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QMainWindow
from ..io import read_auto
from ..io.util import radar_format
from ..draw.SingleRadarPlot import RadarGraph
from ..draw.SingleRadarPlotMap import RadarGraphMap
from ..configure.location_config import last_open_dir
from glob import glob
import json
import sys
from .station_info import Ui_Dialog
from ..draw.VerticalSectionPlot import VerticalSection
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

field_name = ["dBZ", "V", "W", "ZDR", "KDP", "CC"]


class LineBuilder:
    def __init__(self, fig, ax, radar_data, product, map_bool):
        self.ax = ax
        self.xs = []
        self.ys = []
        self.fig = fig
        self.map = map_bool
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self)
        self.cursor = self.fig.canvas.mpl_connect('motion_notify_event', self.mouse_move)
        self.radar_dat = radar_data
        self.product = product

    def __call__(self, event):
        if len(self.xs) < 2:
            self.xs.append(event.xdata)
            self.ys.append(event.ydata)
            if len(self.xs) == 1:
                self.start = self.ax.scatter(event.xdata, event.ydata, color="r", marker="+", s=60,
                                             zorder=len(self.xs)+10)
            else:
                self.end = self.ax.scatter(event.xdata, event.ydata, color="r", marker="+", s=60,
                                             zorder=len(self.xs) + 10)
                self.rline = self.ax.plot(self.xs, self.ys, color="r", linewidth=1, zorder=13)
                cv = FigureCanvas(Figure(figsize=(8, 6)))
                ax = cv.figure.add_axes([0.1, 0.3, 0.8, 0.6])
                cax = cv.figure.add_axes([0.1, 0.1, 0.8, 0.06])
                if not self.map:
                    VerticalSection.GUI_section(cv.figure, ax, cax, self.radar_dat, [self.xs[0]*1000, self.ys[0]*1000],\
                                            [self.xs[1]*1000, self.ys[1]*1000], field_name[self.product])
                else:
                    VerticalSection.GUI_section_map(cv.figure, ax, cax, self.radar_dat,
                                                [self.xs[0], self.ys[0]], \
                                                [self.xs[1], self.ys[1]], field_name[self.product])
                cv.show()
            self.fig.canvas.draw()
        else:
            self.rline[0].remove()
            self.start.remove()
            self.end.remove()
            self.xs = []
            self.ys = []
            self.xs.append(event.xdata)
            self.ys.append(event.ydata)
            self.start = self.ax.scatter(event.xdata, event.ydata, color="r", marker="+", s=60,
                                         zorder=len(self.xs) + 10)
            self.fig.canvas.draw()

    def mouse_move(self, event):
        try:
            self.move_line[0].remove()
        except Exception:
            pass
        if len(self.xs) == 1:
            self.move_line = self.ax.plot([self.xs[0], event.xdata], [self.ys[0], event.ydata], color="r",
                                          linewidth=1, linestyle="--", zorder=100)
        self.fig.canvas.draw()


class Dialog(QDialog, Ui_Dialog):
    """
    Class documentation goes here.
    """

    def __init__(self, parent=None):
        """
        Constructor
        
        @param parent reference to the parent widget
        @type QWidget
        """
        super(Dialog, self).__init__(parent)
        self.setupUi(self)

    @pyqtSlot()
    def on_pushButton_clicked(self):
        """
        Slot documentation goes here.
        """
        self.lon = float(self.lineEdit.text())
        self.lat = float(self.lineEdit_2.text())
        self.height = float(self.lineEdit_3.text())
        self.close()

    @pyqtSlot()
    def on_pushButton_2_clicked(self):
        """
        Slot documentation goes here.
        """
        self.close()

    @pyqtSlot()
    def on_toolButton_clicked(self):
        """
        Slot documentation goes here.
        """
        lon, LonTrue = QInputDialog.getDouble(self, r"经度", "雷达站点经度（单位：度）", 131.3, -180, 180)
        if LonTrue:
            self.lineEdit.setText(str(lon))

    @pyqtSlot()
    def on_toolButton_2_clicked(self):
        """
        Slot documentation goes here.
        """
        # TODO: not implemented yet
        lat, LatTrue = QInputDialog.getDouble(self, r"纬度", "雷达站点纬度（单位：度）", 23, -90, 90)
        if LatTrue:
            self.lineEdit.setText(str(lat))

    @pyqtSlot()
    def on_toolButton_3_clicked(self):
        """
        Slot documentation goes here.
        """
        # TODO: not implemented yet
        height, HeightTrue = QInputDialog.getDouble(self, r"高度", "雷达站点高度（单位：米）", 57, -2000, 5000)
        if HeightTrue:
            self.lineEdit.setText(str(height))


class MainWindow(QMainWindow, Ui_MainWindow):
    """
    Class documentation goes here.
    """

    def __init__(self, parent=None):
        """
        Constructor
        
        @param parent reference to the parent widget
        @type QWidget
        """
        super(MainWindow, self).__init__(parent)
        self.setupUi(self)
        self.lastOpenDir = self.open_last_opendir()
        self.radar_dat = None
        self.dualPOL = False
        self.openbasename = None
        self.files = None
        self.radar_type = None

        self.org_lat = 131.3
        self.org_lon = 23
        self.org_height = 57

    def open_last_opendir(self):
        """打开上次关闭文件的位置"""
        with open(last_open_dir, "r") as f:
            dir_dict = json.load(f)
        return dir_dict["lastOpenDir"]

    def write_last_opendir(self, filedir):
        """将打开的位置写入json文件中"""
        with open(last_open_dir, "w") as f:
            json.dump({"lastOpenDir":filedir}, f)

    @pyqtSlot()
    def on_actionvertical_changed(self):
        """垂直剖面的绘制"""
        if self.actionvertical.isChecked():
            try:
                self.linebuilder = LineBuilder(self.fig, self.ax, self.radar_dat, self.find_var_in_groupBox(),\
                                               self.actionwithmap.isChecked())
                self.clickevent = True
            except AttributeError:
                pass
        else:
            self.fig.canvas.mpl_disconnect(self.linebuilder.cid)
            self.fig.canvas.mpl_disconnect(self.linebuilder.cursor)
            self.linebuilder.rline[0].remove()
            self.linebuilder.start.remove()
            self.linebuilder.end.remove()
            self.fig.canvas.draw()

    @pyqtSlot()
    def on_actionwithmap_changed(self):
        """
        Slot documentation goes here.
        """
        pass

    @pyqtSlot()
    def on_actioncontinuous_changed(self):
        """
        Slot documentation goes here.
        """

        pass

    def Read_radar(self, filename):
        if radar_format(filename) is not None:
            NRadar = read_auto(filename)
            self.org_lat = NRadar.scan_info.latitude.values
            self.org_lon = NRadar.scan_info.longitude.values
            self.org_height = NRadar.scan_info.altitude.values
            if "KDP" in NRadar.fields[0].keys():
                self.open_dual()
            else:
                self.close_non_dual()
            return NRadar
        else:
            QMessageBox.warning(self, "数据错误警告", "非SA/SB/CA/CB/98D/CC/CCJ/SC/CD数据",
                                QMessageBox.Yes)
            return 0

    def close_non_dual(self):
        """关闭非双偏振雷达变量"""
        self.radioButton_13.hide()
        self.radioButton_14.hide()
        self.radioButton_15.hide()

    def open_dual(self):
        """关闭非双偏振雷达变量"""
        self.radioButton_13.show()
        self.radioButton_14.show()
        self.radioButton_15.show()

    def setSelected(self, filename):
        """将选中数据高亮"""
        basename = os.path.basename(filename)
        self.openbasename = basename
        items = self.listWidget.findItems(basename, Qt.MatchExactly)
        if len(items) > 0:
            for item in items:
                self.listWidget.setCurrentItem(item)

    def import_basedat(self, direc):
        """查找文件夹中的所有雷达文件名，并以list返回"""
        self.lastOpenDir = direc
        self.write_last_opendir(direc)
        extensions = ["*.*A", "*.*V", "*.bz2", "*.bin",
                      "*.AR2", "*.gz", ".GZ"]
        files = []
        for iextend in extensions:
            file = glob(os.path.join(direc, iextend))
            files.extend(file)
        return [os.path.basename(ifile) for ifile in files]

    def add_listwidget(self, files):
        """将files添加到listWidget"""
        self.listWidget.clear()
        for item in files:
            self.listWidget.addItem(item)

    @pyqtSlot(QListWidgetItem)
    def on_listWidget_itemDoubleClicked(self, item):
        """
        Slot documentation goes here.
        
        @param item DESCRIPTION
        @type QListWidgetItem
        """
        filename = self.lastOpenDir + os.sep + item.text()
        self.radar_dat = self.Read_radar(filename)
        if self.radar_dat != 0:
            self.setSelected(filename)
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_actionopen_2_triggered(self):
        """
        Slot documentation goes here.
        """
        if self.lastOpenDir and os.path.exists(self.lastOpenDir):
            defaultOpenDirPath = self.lastOpenDir
        else:
            defaultOpenDirPath = '.'
        filename = QFileDialog.getOpenFileName(self, "打开一个雷达基数据", defaultOpenDirPath,
                                               "天气雷达基数据(*bin *bz2 *A *V *BIN *BZ2 *AR2 *GZ *gz)")
        ReadFile = filename[0]
        if ReadFile.strip() == "":
            return
        PathDir = os.path.dirname(ReadFile)
        self.files = self.import_basedat(PathDir)
        self.add_listwidget(self.files)
        self.radar_dat = self.Read_radar(ReadFile)
        if self.radar_dat != 0:
            self.setSelected(ReadFile)
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_actionopendir_2_triggered(self):
        """
        Slot documentation goes here.
        """
        if self.lastOpenDir and os.path.exists(self.lastOpenDir):
            defaultOpenDirPath = self.lastOpenDir
        else:
            defaultOpenDirPath = '.'

        self.targetDirPath = QFileDialog.getExistingDirectory(self, "打开新一代天气雷达数据文件夹",
                                                              defaultOpenDirPath,
                                                              QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks)

        if self.targetDirPath.strip() == '':
            return
        self.files = self.import_basedat(self.targetDirPath)
        self.add_listwidget(self.files)

    @pyqtSlot()
    def on_actionquit_2_triggered(self):
        """
        Slot documentation goes here.
        """
        sys.exit(0)

    @pyqtSlot()
    def on_actionstation_triggered(self):
        """
        Slot documentation goes here.
        """
        self.my_info = Dialog()
        self.my_info.lineEdit.setText(str(self.org_lon))
        self.my_info.lineEdit_2.setText(str(self.org_lat))
        self.my_info.lineEdit_3.setText(str(self.org_height))

        self.my_info.lat = self.org_lat
        self.my_info.lon = self.org_lon
        self.my_info.height = self.org_height

        self.my_info.exec_()
        self.org_lat = self.my_info.lat
        self.org_lon = self.my_info.lon
        self.org_height = self.my_info.height

    def find_checked_radiobutton(self, radiobuttons):
        ''' find the checked radiobutton '''
        for items in radiobuttons:
            if items.isChecked():
                checked_radiobutton = items.text()
                return checked_radiobutton

    def find_level_in_groupBox(self):
        """查找仰角"""
        level = self.find_checked_radiobutton(self.groupBox.findChildren(QtWidgets.QRadioButton))
        levels = ["第1层", "第2层", "第3层",
                  "第4层", "第5层", "第6层",
                  "第7层", "第8层", "第9层"]
        for i in range(9):
            if level == levels[i]:
                return i
        return 0

    def find_var_in_groupBox(self):
        """查找变量"""
        var = self.find_checked_radiobutton(self.groupBox_2.findChildren(QtWidgets.QRadioButton))
        vars = ["反射率因子", "径向速度", "谱宽", "差分反射率", "差分相位比", "相关系数"]
        for i in range(6):
            if var == vars[i]:
                return i
        return 0

    def plot_graph_PPI(self, radar, level, product, map, continuously):
        self.MplWidget.canvas.update()
        self.MplWidget.canvas.flush_events()
        try:
            self.fig.clf()
            self.ax.clear()
            self.cax.clear()
        except AttributeError:
            pass

        if not map:
            self.fig, self.ax, self.cax = self.MplWidget.canvas.get_fig_ax()
            self.ax.set_facecolor((0.95, 0.95, 0.95))
            self.pm = RadarGraph.GUI_plot(radar, self.fig, self.ax, self.cax, level, \
                                field_name[product], continuously=continuously)
        else:
            self.fig, self.ax, self.cax = self.MplWidget.canvas.get_fig_ax_map()
            self.ax.set_facecolor((0.95, 0.95, 0.95))
            self.pm = RadarGraphMap.GUI_plot(radar, self.fig, self.ax, self.cax, level, \
                                   field_name[product], continuously=continuously)
        self.ax.tick_params(axis="y", which="both", direction='in')
        self.ax.tick_params(axis="x", which="both", direction='in')
        self.MplWidget.canvas.draw()

        if self.actionvertical.isChecked(): #尝试重新绑定
            try:
                self.fig.canvas.mpl_disconnect(self.linebuilder.cid)
                self.fig.canvas.mpl_disconnect(self.linebuilder.cursor)
                self.linebuilder = LineBuilder(self.fig, self.ax, self.radar_dat, self.find_var_in_groupBox(), \
                                               self.actionwithmap.isChecked())
                self.clickevent = True
            except AttributeError:
                pass

    @pyqtSlot()
    def on_pushButton_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.files is not None:
            items = self.listWidget.findItems(self.openbasename, Qt.MatchExactly)
            row = self.listWidget.row(items[0])
            nrows = len(self.files)
            res_row = row - 1
            if res_row < 0:
                res_row = nrows - 1
            self.radar_dat = self.Read_radar(self.lastOpenDir + os.sep + self.files[res_row])
            if self.radar_dat != 0:
                self.setSelected(self.lastOpenDir + os.sep + self.files[res_row])
                self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                    self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_pushButton_2_clicked(self):
        """
        Slot documentation goes here.
        动画播放部分
        """
        if self.files is not None:
            items = self.listWidget.findItems(self.openbasename, Qt.MatchExactly)
            row = self.listWidget.row(items[0])
            nrows = len(self.files)
            for irow in range(row, nrows):
                self.radar_dat = self.Read_radar(os.path.join(self.lastOpenDir, self.files[irow]))
                if self.radar_dat != 0:
                    self.setSelected(self.lastOpenDir + os.sep + self.files[irow])
                    self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                        self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_pushButton_3_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.files is not None:
            items = self.listWidget.findItems(self.openbasename, Qt.MatchExactly)
            row = self.listWidget.row(items[0])
            nrows = len(self.files)
            res_row = row + 1
            if res_row == nrows:
                res_row = 0
            self.radar_dat = self.Read_radar(self.lastOpenDir + os.sep + self.files[res_row])
            if self.radar_dat != 0:
                self.setSelected(self.lastOpenDir + os.sep + self.files[res_row])
                self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                    self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_15_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_12_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_14_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_10_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_13_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_11_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_2_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_4_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_5_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_3_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_1_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_7_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_8_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_6_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())

    @pyqtSlot()
    def on_radioButton_9_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(), self.actioncontinuous.isChecked())


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    ui = MainWindow()
    ui.show()
    sys.exit(app.exec_())
