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

import sys
import os

import numpy as np

import matplotlib
matplotlib.use("Qt5Agg")
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow
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
        
        self.ui = LASDiAGUI.Ui_LASDiAGUI()
        self.ui.setupUi(self)
        
        # Set the variable
        self.Q = None
        self.I_Q = None
        self.Qbkg = None
        self.Ibkg_Q = None
        self.SsmoothDamp_Q = None
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
        self.ui.calcFr.clicked.connect(self.Fr)
        self.ui.optimize.clicked.connect(self.Optimization)
        #self.ui.Minimization.clicked.connect(self.calcMinimization)
        
    #---------------------------------------------------------  

    def import_data(self):
        """Function to load and plot the data file"""
        
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load Data File", \
            r"C:\Users\devoto\work\ID27\data\cea_files\Ar", "Data File(*chi *xy)")
        
        self.Q, self.I_Q = Utility.read_file(path)
        
        self.ui.fileSampleName.setPlainText(path)
        
        self.ui.rawDataPlot.canvas.ax.plot(self.Q, self.I_Q, label="Data")
        self.ui.rawDataPlot.canvas.ax.legend()
        self.ui.rawDataPlot.canvas.draw()

    #---------------------------------------------------------
        
    def import_bkg(self):
        """Function to load and plot the bkg file"""
        
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load Bkg File", \
            r"C:\Users\devoto\work\ID27\data\cea_files\Ar", "Data File(*chi *xy)")
        
        self.Qbkg, self.Ibkg_Q = Utility.read_file(path)
        
        self.ui.fileBkgName.setPlainText(path)
        
        self.ui.rawDataPlot.canvas.ax.plot(self.Qbkg, self.Ibkg_Q, "g--", label="Bkg")
        self.ui.rawDataPlot.canvas.ax.legend()
        self.ui.rawDataPlot.canvas.draw()

    #---------------------------------------------------------

    def SQ(self):
        """Function to calculate and plot the structure factor S(Q)"""
        
        elementList = Utility.molToelemList("Ar")
        elementParameters = Utility.read_parameters(elementList, "./elementParameters.txt")
        
        self.Q, self.I_Q, self.Qbkg, self.Ibkg_Q = UtilityAnalysis.check_data_length(self.Q, \
            self.I_Q, self.Qbkg, self.Ibkg_Q, self.ui.minQ.value(), self.ui.maxQ.value())
        
        self.fe_Q, self.Ztot = MainFunctions.calc_eeff(elementList, self.Q, elementParameters)
        self.Iincoh_Q = MainFunctions.calc_Iincoh(elementList, self.Q, elementParameters)
        self.J_Q = MainFunctions.calc_JQ(self.Iincoh_Q, self.Ztot, self.fe_Q)
        self.Sinf, Sinf_Q = MainFunctions.calc_Sinf(elementList, self.fe_Q, \
            self.Q, self.Ztot, elementParameters)
        
        Isample_Q = MainFunctions.calc_IsampleQ(self.I_Q, self.ui.sfValue.value(), self.Ibkg_Q)
        alpha = MainFunctions.calc_alpha(self.J_Q[self.Q<=self.ui.QmaxInt.value()], self.Sinf, \
            self.Q[self.Q<=self.ui.QmaxInt.value()], Isample_Q[self.Q<=self.ui.QmaxInt.value()], \
        self.fe_Q[self.Q<=self.ui.QmaxInt.value()], self.Ztot, self.ui.rho0Value.value())
        Icoh_Q = MainFunctions.calc_Icoh(alpha, Isample_Q, self.Iincoh_Q)
        
        S_Q = MainFunctions.calc_SQ(Icoh_Q, self.Ztot, self.fe_Q, self.Sinf, self.Q, \
            self.ui.minQ.value(), self.ui.QmaxInt.value(), self.ui.maxQ.value())
        Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(self.Q, S_Q, self.Sinf, \
            self.ui.smoothFactor.value(), self.ui.minQ.value(), self.ui.QmaxInt.value(), \
            self.ui.maxQ.value())
        self.SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, self.Q, self.Sinf, \
            self.ui.QmaxInt.value(), self.ui.dampFactor.value())
        
        self.ui.factorPlot.canvas.ax.plot(self.Q, self.SsmoothDamp_Q, label="S(Q)")
        self.ui.factorPlot.canvas.draw()

    #---------------------------------------------------------

    def Fr(self):
        """Function to calculte and plot F(r)"""
        
        self.i_Q = MainFunctions.calc_iQ(self.SsmoothDamp_Q, self.Sinf)
        self.r = MainFunctions.calc_r(self.Q)
        self.F_r = MainFunctions.calc_Fr(self.r, self.Q[self.Q<=self.ui.QmaxInt.value()], \
            self.i_Q[self.Q<=self.ui.QmaxInt.value()])
        
        self.ui.distfuncPlot.canvas.ax.plot(self.r, self.F_r, label="F(r)")
        self.ui.distfuncPlot.canvas.draw()
    
    #---------------------------------------------------------

    def Optimization(self):
        """Function to optimize and plot F(r)"""
        
        Fintra_r = np.zeros(self.r.size)
        
        
        F_rIt, deltaF_rIt = Optimization.calc_optimize_Fr(self.ui.iterations.value(), self.F_r, \
            Fintra_r, self.ui.rho0Value.value(), self.i_Q[self.Q<=self.ui.QmaxInt.value()], \
            self.Q[self.Q<=self.ui.QmaxInt.value()], self.Sinf, \
            self.J_Q[self.Q<=self.ui.QmaxInt.value()], self.r, self.ui.rmin.value(), "n")
        
        self.ui.distfuncPlot.canvas.ax.plot(self.r, F_rIt, label=r"F_{opt}(r)")
        self.ui.distfuncPlot.canvas.draw()
    
    #---------------------------------------------------------
    
    # def importXYZFile(self):
        # """Function to import the XYZ file and calculate the intramolecular
           # component"""
        
        # path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load XYZ File", \
            # r".\xyzFiles", "XYZ File(*xyz)")
        
        # numAtoms, element, x, y, z = Utility.read_xyz_file(path)
        
        # r, iintradamp_Q, Fintra_r = Optimization.calc_intraComponent(self.Q, fe_Q, Ztot, \
            # variables.QmaxIntegrate, variables.maxQ, elementList, element, \
        # x, y, z, elementParameters, variables.damping_factor)
    
    #---------------------------------------------------------
    
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
