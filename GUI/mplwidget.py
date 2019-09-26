# Imports
import matplotlib
# Ensure using PyQt5 backend
matplotlib.use('QT5Agg')
from PyQt5 import QtWidgets
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as Canvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

# Matplotlib canvas class to create figure
class MplCanvas(Canvas):
    def __init__(self):
        self.fig = Figure()
        Canvas.__init__(self, self.fig)
        Canvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        Canvas.updateGeometry(self)
        
    def get_fig_ax(self):
        self.ax = self.fig.add_axes([0.035, 0.06, 0.88, 0.88])
        self.cax = self.fig.add_axes([0.9, 0.06, 0.028, 0.88])
        return self.fig,  self.ax,  self.cax
        
# Matplotlib widget
class MplWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        QtWidgets.QWidget.__init__(self, parent)   # Inherit from QWidget
        self.canvas = MplCanvas()                  # Create canvas object
        self.ntb = NavigationToolbar(self.canvas,  parent)
        self.vbl = QtWidgets.QVBoxLayout()         # Set box for plotting
        self.vbl.addWidget(self.canvas)
        self.vbl.addWidget(self.ntb)
        self.setLayout(self.vbl)
