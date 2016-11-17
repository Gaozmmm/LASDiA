import sys
import matplotlib
# matplotlib.use("Qt5Agg")
from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QVBoxLayout, QSizePolicy, QMessageBox, QWidget
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class MplCanvas(FigureCanvas):
    
    def __init__(self):
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding,QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


class MplWidget(QWidget):
    
    def __init__(self, parent = None):
        QWidget.__init__(self, parent)
        self.canvas = MplCanvas()
        self.vbl = QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)

# class MplWidget(FigureCanvas):
    # """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    # def __init__(self, parent=None, width=5, height=4, dpi=100):
        # fig = Figure(figsize=(width, height), dpi=dpi)
        # self.axes = fig.add_subplot(111)
        # # We want the axes cleared every time plot() is called
        # self.axes.hold(False)
        
        # self.compute_initial_figure()
        
        # FigureCanvas.__init__(self, fig)
        # self.setParent(parent)
        
        # FigureCanvas.setSizePolicy(self,
                # QSizePolicy.Expanding,
                # QSizePolicy.Expanding)
        # FigureCanvas.updateGeometry(self)
    
    # def compute_initial_figure(self):
        # pass
