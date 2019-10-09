import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from NuistRadar.GraphicalInterface.RadarInterface import MainWindow
from PyQt5 import QtWidgets
app = QtWidgets.QApplication(sys.argv)
ui = MainWindow()
ui.show()
sys.exit(app.exec_())