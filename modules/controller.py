# The MIT License (MIT)

# Copyright (c) 2015-2016 European Synchrotron Radiation Facility

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""LASDiA Controller
"""

from __future__ import (absolute_import, division, print_function, unicode_literals)
import six

from PyQt5 import QtCore, QtGui, QtWidgets

import sys
import os

import matplotlib
matplotlib.use("Qt4Agg")
import matplotlib.pyplot as plt
from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QVBoxLayout, QSizePolicy, QMessageBox, QWidget
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
# from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
# import numpy as np

from modules import Formalism
from modules import Geometry
from modules import IgorFunctions
from modules import KaplowMethod
from modules import MainFunctions
from modules import Minimization
from modules import Optimization
from modules import Utility
from modules import UtilityAnalysis

from modules import LASDiAGUI

class LASDiA(QtWidgets.QMainWindow, LASDiAGUI.Ui_LASDiAGUI):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        #self.toolbar = NavigationToolbar(self.canvas, self)
        self.ui = LASDiAGUI.Ui_LASDiAGUI()
        self.ui.setupUi(self)

        # Set the variable
        self.Q = None
        self.I_Q = None
        self.Qbkg = None
        self.I_Qbkg = None
        self.S_Q = None
        self.elementList = None
        self.Ztot = None
        self.fe_Q = None
        self.r = None
        self.F_r = None
        self.Sinf = None
        self.i_Q = None
        self.rho0 = None
        self.Iincoh_Q = None
        self.J_Q = None
        
        # Set the buttons
        self.ui.importData.clicked.connect(self.import_data)
        # self.ui.LoadBkg.clicked.connect(self.load_bkg)
        # self.ui.CalcSQ.clicked.connect(self.calcSQ)
        # self.ui.Calcgr.clicked.connect(self.calcgr)
        # self.ui.Optimization.clicked.connect(self.calcOptimization)
        #self.ui.Minimization.clicked.connect(self.calcMinimization)
        
    #---------------------------------------------------------  

    def import_data(self):
        '''load and plot the file'''
        # Open data file
        path = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', '.')
        
        self.Q, self.I_Q = Utility.read_file(path)

        self.ui.matplotlibwidget.canvas.ax.clear()
        self.ui.matplotlibwidget.canvas.ax.plot(self.Q, self.I_Q, label='Data')
        self.ui.matplotlibwidget.canvas.draw()

    #---------------------------------------------------------
        
    # def load_bkg(self):
        # '''load and plot the file'''
        # # Open bkg file
        # path = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '.')
        
        # # Modify the variables as numpy array
        # self.Qbkg, self.I_Qbkg = read_file(path)

        # self.ui.RawData.canvas.ax.plot(self.Qbkg, self.I_Qbkg, 'g--', label='Bkg')
        # self.ui.RawData.canvas.draw()

    # #---------------------------------------------------------

    # def calcSQ(self):
        # element = str(self.ui.Elements.currentText())
        # multiplicity = self.ui.Multiplicity.value()
        # self.elementList = {element:multiplicity}
    
        # self.fe_Q, self.Ztot = calc_eeff(self.elementList, self.Q)
    
        # #numAtoms =  sc.N_A
        # N = 1
        # if self.ui.sValue.value() == 0.00:
            # s = 1
        # else:
            # s = self.ui.sValue.value()
        
        # if self.ui.rho0Value.value() == 0.00:
            # self.rho0 = 25.0584
        # else:
            # self.rho0 = self.ui.rho0Value.value()

        # self.Iincoh_Q = calc_Iincoh(self.elementList, self.Q)            
        # self.J_Q = calc_JQ(self.Iincoh_Q, self.Ztot, self.fe_Q)
        # self.Sinf = calc_Sinf(self.elementList, self.fe_Q, self.Q, self.Ztot)

        # Isample_Q = self.I_Q - s * self.I_Qbkg
        
        # alpha = calc_alpha(self.J_Q, self.Sinf, self.Q, Isample_Q, self.fe_Q, self.Ztot, self.rho0)

        # Icoh_Q = calc_Icoh(N, alpha, Isample_Q, self.Iincoh_Q)
    
        # self.S_Q = calc_SQ(N, Icoh_Q, self.Ztot, self.fe_Q, self.Q)

        # qmax = np.amax(self.Q)
        # Qa = np.arange(qmax,100,self.Q[1]-self.Q[0])
        # newQ = np.concatenate([self.Q, Qa])
        # SQa = np.zeros(Qa.size)
        # SQa.fill(1)
        # newSQ = np.concatenate([self.S_Q, SQa])
        
        # self.ui.SQ.canvas.ax.clear()
        # self.ui.SQ.canvas.ax.plot(newQ, self.S_Q)
        # self.ui.SQ.canvas.draw()

    # #---------------------------------------------------------

    # def calcgr(self):
        # self.i_Q = calc_iQ(self.S_Q, self.Sinf)
        # self.r = calc_spectrum(self.i_Q)
        # self.F_r = calc_Fr(self.r, self.Q, self.i_Q)
        # #self.rho0 = calc_rho0(self.Q, self.i_Q)
        # #g_r = calc_gr(r, F_r, rho0)

        # # self.ui.gr.canvas.ax.clear()
        # # self.ui.gr.canvas.ax.plot(r, g_r)
        # # self.ui.gr.canvas.draw()

        # self.ui.Fr.canvas.ax.clear()
        # self.ui.Fr.canvas.ax.plot(self.r, self.F_r)
        # self.ui.Fr.canvas.draw()
    
    # #---------------------------------------------------------

    # def calcOptimization(self):
        # iteration = self.ui.Iteration.value()
        # r_cutoff = self.ui.rcutoff.value()

        # Fintra_r = calc_Fintra()
        # optF_r = calc_optimize_Fr(iteration, self.F_r, Fintra_r, self.rho0, self.i_Q, self.Q, self.Sinf, self.J_Q, self.r, r_cutoff)

        # #self.ui.Fr.canvas.ax.clear()
        # self.ui.Fr.canvas.ax.plot(self.r, optF_r)
        # self.ui.Fr.canvas.draw()

        # # SQ_F = calc_SQ_F(optF_r, self.r, self.Q, self.Sinf)

        # # self.ui.SQ.canvas.ax.plot(self.Q, SQ_F)
        # # self.ui.SQ.canvas.draw()
            
    # #---------------------------------------------------------

# #    def calcMinimization(self):
        
    
    # #---------------------------------------------------------
    
def main():
    # app = QtWidgets.QApplication(sys.argv)
    # LASDiA = LASDiAGUI.QtWidgets.QMainWindow()
    # ui = LASDiAGUI.Ui_LASDiAGUI()
    # ui.setupUi(LASDiA)
    # LASDiA.show()
    # sys.exit(app.exec_())
    
    app = QtWidgets.QApplication(sys.argv)
    ui = LASDiA()
    ui.show()
    sys.exit(app.exec_())
    

if __name__ == '__main__':
    main()
