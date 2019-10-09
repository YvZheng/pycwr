from NuistRadar.GraphicalInterface.RadarInterface import MainWindow
import sys
from PyQt5 import QtWidgets

app = QtWidgets.QApplication(sys.argv)
ui = MainWindow()
ui.show()
sys.exit(app.exec_())