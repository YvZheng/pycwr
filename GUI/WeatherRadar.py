# -*- coding: utf-8 -*-

"""
Module implementing MainWindow.
"""
import os
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5 import QtWidgets
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QMainWindow
from Ui_WeatherRadar import Ui_MainWindow
from SCFileRead import IO
from CINRADPlot import CINRAD
from CPolPlot import CPol
from CPolRead import CIO
from glob import glob
import sys
import bz2
import gzip

from Ui_info import Ui_Dialog


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
        lon,  LonTrue = QInputDialog.getDouble(self, r"经度", "雷达站点经度（单位：度）",  131.3, -180, 180)
        if LonTrue:
            self.lineEdit.setText(str(lon))
    
    @pyqtSlot()
    def on_toolButton_2_clicked(self):
        """
        Slot documentation goes here.
        """
        # TODO: not implemented yet
        lat,  LatTrue = QInputDialog.getDouble(self, r"纬度", "雷达站点纬度（单位：度）",  23, -90, 90)
        if LatTrue:
            self.lineEdit.setText(str(lat))
    
    @pyqtSlot()
    def on_toolButton_3_clicked(self):
        """
        Slot documentation goes here.
        """
        # TODO: not implemented yet
        height,  HeightTrue = QInputDialog.getDouble(self, r"高度", "雷达站点高度（单位：米）",  57, -2000, 5000)
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
        self.lastOpenDir = False 
        self.radar_dat = None
        self.dualPOL = False
        self.openbasename = None
        self.files = None
        self.radar_type = None
        
        self.org_lat = 23
        self.org_lon = 131.3
        self.org_height = 57

    @pyqtSlot()
    def on_actionwithmap_changed(self):
        """
        Slot documentation goes here.
        """

        print(self.actionwithmap.isChecked())
    
    @pyqtSlot()
    def on_actioncontinuous_changed(self):
        """
        Slot documentation goes here.
        """

        print(self.actioncontinuous.isChecked())
        
    def radar_format(self, filename):
        if hasattr(filename, 'read'):
            return filename

        fh = open(filename, 'rb')
        magic = fh.read(3)
        fh.close()

        if magic.startswith(b'\x1f\x8b'):
            f = gzip.GzipFile(filename, 'rb')
        elif magic.startswith(b'BZh'):
            f = bz2.BZ2File(filename, 'rb')
        else:
            f = open(filename, 'rb')

        cpol_id = f.read(4)
        sc_id = f.read(24)
        f.close()
        if cpol_id == b'RSTM':
            self.open_dual()
            return 1
        elif sc_id[10:12] == b'\x01\x00':
            self.close_non_dual()
            return 0
        else:
            return 2

    def Read_radar(self,  filename):
        self.radar_type = self.radar_format(filename)
        if self.radar_type == 0:
            return IO.read_sc(filename,  self.org_lat,  self.org_lon,  self.org_height)
        elif self.radar_type == 1:
            return CIO.read(filename)
        else:
            QMessageBox.warning(self,"数据错误警告",  "非SA/SB/CA/CB/POL数据", 
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
        
    def setSelected(self,  filename):
        """将选中数据高亮"""
        basename = os.path.basename(filename)
        self.openbasename  = basename
        items = self.listWidget.findItems(basename,Qt.MatchExactly)
        if len(items) > 0:
            for item in items:
                self.listWidget.setCurrentItem(item)
    
    def import_basedat(self,  dir):
        """查找文件夹中的所有雷达文件名，并以list返回"""
        self.lastOpenDir = dir
        extensions = ["*.*A",  "*.*V", "*.bz2",  "*.bin",
                      "*.AR2", "*.gz",".GZ"]
        files = []
        for iextend in extensions:
            file = glob(os.path.join(dir,  iextend))
            files.extend(file)
        return [os.path.basename(ifile) for ifile in files]
        
    def add_listwidget(self,  files):
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
                                self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
            
    @pyqtSlot()
    def on_actionopen_2_triggered(self):
        """
        Slot documentation goes here.
        """
        if self.lastOpenDir and os.path.exists(self.lastOpenDir):
            defaultOpenDirPath = self.lastOpenDir
        else:
            defaultOpenDirPath = '.'
        filename = QFileDialog.getOpenFileName(self,"打开一个雷达基数据", defaultOpenDirPath, 
                        "天气雷达基数据(*bin *bz2 *A *V *BIN *BZ2 *AR2 *GZ *gz)")
        ReadFile = filename[0]
        if ReadFile.strip() == "":
            return
        PathDir = os.path.dirname(ReadFile)
        print("Path", PathDir)
        self.files = self.import_basedat(PathDir)
        self.add_listwidget(self.files)
        self.radar_dat = self.Read_radar(ReadFile)
        if self.radar_dat != 0:
            self.setSelected(ReadFile)
            self.plot_graph_PPI(self.radar_dat,  self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
        print(self.radioButton_2.isChecked())
        print(ReadFile)
    
    @pyqtSlot()
    def on_actionopendir_2_triggered(self):
        """
        Slot documentation goes here.
        """
        if self.lastOpenDir and os.path.exists(self.lastOpenDir):
            defaultOpenDirPath = self.lastOpenDir
        else:
            defaultOpenDirPath = '.'
    
        self.targetDirPath = QFileDialog.getExistingDirectory(self,"打开新一代天气雷达数据文件夹", 
                                    defaultOpenDirPath, QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks)
                                    
        if self.targetDirPath.strip()=='':
            return
        print("Path", self.targetDirPath)
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
        vars = ["反射率因子", "径向速度", "谱宽", "差分反射率",  "差分相位比", "相关系数"]
        for i in range(6):
            if var == vars[i]:
                return i
        return 0
        
    def plot_graph_PPI(self,  radar,  level, product,  map,  nws):
        self.MplWidget.canvas.update()
        self.MplWidget.canvas.flush_events()
        fig,  ax,  cax = self.MplWidget.canvas.get_fig_ax()
        #fig.subplots_adjust(top = 0.95, bottom = 0.055, right = 1, left = 0.03, hspace = 0, wspace = 0)
        #ax.margins(0, 0)
        ax.clear()
        cax.clear()
        ax.set_facecolor((0.95, 0.95, 0.95))
        if self.radar_type == 0:
            if map:
                CINRAD.PPI_MAP(fig, ax, cax, radar, level, product, nws=nws)
            else:
                CINRAD.PPI(fig,  ax,  cax,  radar,  level,  product, nws=nws)
        else:
            if map:
                CPol.PPI_MAP(fig,  ax,  cax,  radar,  level,  product, nws=nws)
            else:
                CPol.PPI(fig, ax, cax, radar, level, product, nws=nws)
        self.MplWidget.canvas.draw()
        
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
                self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(),  self.find_var_in_groupBox(),
                                self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_pushButton_2_clicked(self):
        """
        Slot documentation goes here.
        """
        pass
    
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
                                self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
            
    @pyqtSlot()
    def on_radioButton_15_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_radioButton_12_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_radioButton_14_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_radioButton_10_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_radioButton_13_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_radioButton_11_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_radioButton_2_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_radioButton_4_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_radioButton_5_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_radioButton_3_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_radioButton_1_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_radioButton_7_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_radioButton_8_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_radioButton_6_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())
    
    @pyqtSlot()
    def on_radioButton_9_clicked(self):
        """
        Slot documentation goes here.
        """
        if self.radar_dat is not None:
            self.plot_graph_PPI(self.radar_dat, self.find_level_in_groupBox(), self.find_var_in_groupBox(),
                            self.actionwithmap.isChecked(),self.actioncontinuous.isChecked())

if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    ui = MainWindow()
    ui.show()
    sys.exit(app.exec_())
    

