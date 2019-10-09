import matplotlib
import sys
sys.path.append("../")
# Ensure using PyQt5 backend
matplotlib.use('QT5Agg')
from NuistRadar.GraphicalInterface.RadarInterface import MainWindow
import sys
from PyQt5 import QtWidgets
app = QtWidgets.QApplication(sys.argv)
ui = MainWindow()
ui.show()
sys.exit(app.exec_())