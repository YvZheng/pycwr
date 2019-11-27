# -*- coding: utf-8 -*-
###after install pycwr lib
import warnings
warnings.filterwarnings('ignore')
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pycwr.GraphicalInterface.RadarInterface import MainWindow
from PyQt5 import QtWidgets
app = QtWidgets.QApplication(sys.argv)
ui = MainWindow()
ui.show()
sys.exit(app.exec_())