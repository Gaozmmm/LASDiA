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
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QVBoxLayout, QSizePolicy, QMessageBox, QWidget
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

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
        self.Ibkg_Q = None
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
        self.ui.importBkg.clicked.connect(self.import_bkg)
        self.ui.calcSQ.clicked.connect(self.SQ)
        # self.ui.Calcgr.clicked.connect(self.calcgr)
        # self.ui.Optimization.clicked.connect(self.calcOptimization)
        #self.ui.Minimization.clicked.connect(self.calcMinimization)
        
    #---------------------------------------------------------  

    def import_data(self):
        '''load and plot the data file'''
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load Data File", r"C:\Users\devoto\work\ID27\data\cea_files\Ar", "Data File(*chi *xy)")
        
        self.Q, self.I_Q = Utility.read_file(path)

        self.ui.rawDataPlot.canvas.ax.clear()
        self.ui.rawDataPlot.canvas.ax.plot(self.Q, self.I_Q, label='Data')
        self.ui.rawDataPlot.canvas.draw()

    #---------------------------------------------------------
        
    def import_bkg(self):
        '''load and plot the bkg file'''
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load Bkg File", r"C:\Users\devoto\work\ID27\data\cea_files\Ar", "Data File(*chi *xy)")
        
        self.Qbkg, self.Ibkg_Q = Utility.read_file(path)

        self.ui.rawDataPlot.canvas.ax.plot(self.Qbkg, self.Ibkg_Q, 'g--', label='Bkg')
        self.ui.rawDataPlot.canvas.draw()

    # #---------------------------------------------------------

    def SQ(self):
        elementList = Utility.molToelemList("Ar")
        elementParameters = Utility.read_parameters(elementList, "./elementParameters.txt")
        
        self.Q, self.I_Q, self.Qbkg, self.Ibkg_Q = UtilityAnalysis.check_data_length(self.Q, \
            self.I_Q, self.Qbkg, self.Ibkg_Q, self.ui.minQ.value(), self.ui.maxQ.value())
        
        self.fe_Q, self.Ztot = MainFunctions.calc_eeff(elementList, self.Q, elementParameters)
        self.Iincoh_Q = MainFunctions.calc_Iincoh(elementList, self.Q, elementParameters)
        self.J_Q = MainFunctions.calc_JQ(self.Iincoh_Q, self.Ztot, self.fe_Q)
        self.Sinf, Sinf_Q = MainFunctions.calc_Sinf(elementList, self.fe_Q, \
            self.Q, self.Ztot, elementParameters)
        
        QmaxIntegrate = 90.0

        Isample_Q = MainFunctions.calc_IsampleQ(self.I_Q, self.ui.sf_value.value(), self.Ibkg_Q)
        alpha = MainFunctions.calc_alpha(self.J_Q[self.Q<=QmaxIntegrate], self.Sinf, \
            self.Q[self.Q<=QmaxIntegrate], Isample_Q[self.Q<=QmaxIntegrate], \
        self.fe_Q[self.Q<=QmaxIntegrate], self.Ztot, self.ui.rho0_value.value())
        Icoh_Q = MainFunctions.calc_Icoh(alpha, Isample_Q, self.Iincoh_Q)

        self.S_Q = MainFunctions.calc_SQ(Icoh_Q, self.Ztot, self.fe_Q, self.Sinf, self.Q, self.ui.minQ.value(), \
            QmaxIntegrate, self.ui.maxQ.value())
        S_Qsmoothed = UtilityAnalysis.calc_SQsmoothing(self.Q, self.S_Q, self.Sinf, self.ui.smooth_factor.value(), \
            self.ui.minQ.value(), QmaxIntegrate, self.ui.maxQ.value())
        S_QsmoothedDamp = UtilityAnalysis.calc_SQdamp(S_Qsmoothed, self.Q, self.Sinf, \
            QmaxIntegrate, self.ui.damping_factor.value())
        
        self.ui.factorPlot.canvas.ax.clear()
        self.ui.factorPlot.canvas.ax.plot(self.Q, S_QsmoothedDamp)
        self.ui.factorPlot.canvas.draw()

    #---------------------------------------------------------

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
