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
            r"../data/cea_files/Ar", "Data File(*chi *xy)")
            # r"C:\Users\devoto\work\ID27\data\cea_files\Ar", "Data File(*chi *xy)")

        
        self.Q, self.I_Q = Utility.read_file(path)
        
        self.ui.fileSampleName.setPlainText(path)
        
        self.ui.rawDataPlot.canvas.ax.plot(self.Q, self.I_Q, label="Data")
        self.ui.rawDataPlot.canvas.ax.legend()
        self.ui.rawDataPlot.canvas.draw()

    #---------------------------------------------------------
        
    def import_bkg(self):
        """Function to load and plot the bkg file"""
        
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load Bkg File", \
            r"../data/cea_files/Ar", "Data File(*chi *xy)")
        	# r"C:\Users\devoto\work\ID27\data\cea_files\Ar", "Data File(*chi *xy)")
        
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
        
        # density, scaleFactor = Minimization.chi2_minimization(self.ui.sfValue.value(), 
            # self.Q, self.I_Q, self.Ibkg_Q, 
            # self.J_Q, self.fe_Q, self.Iincoh_Q, self.Sinf, self.Ztot,
            # self.ui.rho0Value.value(), Fintra_r, self.r, self.ui.minQ.value(), self.ui.QmaxInt.value(), variables.maxQ, 
            # variables.smoothFactor, variables.dampFactor, variables.iteration, variables.rmin)
        
        scaleStep = 0.05
        densityStep = 0.025
        numSample = 23
        numLoopIteration = 0
        
        # plt.ion()
        # figure, ax = plt.subplots()
        
        scaleFactor = self.ui.sfValue.value()
        density = self.ui.rho0Value.value()
        
        while True:
            self.ui.chi2_plot.canvas.ax.cla()
            self.ui.chi2_plot.canvas.ax.grid(True)
            scaleArray = UtilityAnalysis.make_array_loop(scaleFactor, scaleStep, numSample)
            
            chi2Array = np.zeros(numSample)
            
            self.ui.chi2_plot.canvas.ax.set_xlabel("Scale")
            self.ui.chi2_plot.canvas.ax.relim()
            self.ui.chi2_plot.canvas.ax.autoscale_view()
            for i in range(len(scaleArray)):
                chi2Array[i], SsmoothDamp_Q, F_r, Fopt_r = KaplowMethod.Kaplow_method(self.Q, self.I_Q,
                    self.Ibkg_Q, self.J_Q, self.fe_Q, self.Iincoh_Q, self.Sinf, self.Ztot, scaleArray[i], density, Fintra_r, self.r,
                    self.ui.minQ.value(), self.ui.QmaxInt.value(), self.ui.maxQ.value(), self.ui.smoothFactor.value(),
                    self.ui.dampFactor.value(), self.ui.iterations.value(), self.ui.rmin.value())
                
                self.ui.chi2_plot.canvas.ax.scatter(scaleArray[i], chi2Array[i])
                self.ui.chi2_plot.canvas.draw()
            
            xfit, yfit, scaleFactor = Minimization.chi2_fit(scaleArray, chi2Array)
            self.ui.chi2_plot.canvas.ax.plot(xfit, yfit)
            self.ui.chi2_plot.canvas.draw()
            
            self.ui.chi2_plot.canvas.ax.cla()
            self.ui.chi2_plot.canvas.ax.grid(True)

            density0 = density
            densityArray = UtilityAnalysis.make_array_loop(density, densityStep, numSample)
            chi2Array = np.zeros(numSample)
            
            self.ui.chi2_plot.canvas.ax.set_xlabel("Density")
            self.ui.chi2_plot.canvas.ax.relim()
            self.ui.chi2_plot.canvas.ax.autoscale_view()
            for i in range(len(densityArray)):
                chi2Array[i], SsmoothDamp_Q, F_r, Fopt_r = KaplowMethod.Kaplow_method(self.Q, self.I_Q,
                    self.Ibkg_Q, self.J_Q, self.fe_Q, self.Iincoh_Q, self.Sinf, self.Ztot, scaleFactor, densityArray[i], Fintra_r, self.r,
                    self.ui.minQ.value(), self.ui.QmaxInt.value(), self.ui.maxQ.value(), self.ui.smoothFactor.value(),
                    self.ui.dampFactor.value(), self.ui.iterations.value(), self.ui.rmin.value())
                
                self.ui.chi2_plot.canvas.ax.scatter(densityArray[i], chi2Array[i])
                self.ui.chi2_plot.canvas.draw()
                
            
            xfit, yfit, density = Minimization.chi2_fit(densityArray, chi2Array)
            self.ui.chi2_plot.canvas.ax.plot(xfit, yfit)
            self.ui.chi2_plot.canvas.draw()
            
            if np.abs(density-density0) > density0/25:
                scaleStep = 0.006
                densityStep = density0/10
            elif np.abs(density-density0) > density0/75:
                scaleStep = 0.0006
                densityStep = density0/100
            else:
                scaleStep = 0.00006
                densityStep = density0/1000

            numLoopIteration += 1
            if (np.abs(density-density0) < density0/2500 or numLoopIteration > 30):
                break
        
        S_Q = UtilityAnalysis.S_QCalculation(self.Q, self.I_Q, self.Ibkg_Q, scaleFactor, 
            self.J_Q, self.Sinf, self.fe_Q, self.Ztot, density, self.Iincoh_Q, 
            self.ui.minQ.value(), self.ui.QmaxInt.value(), self.ui.maxQ.value(), self.ui.smoothFactor.value(), self.ui.dampFactor.value())
        
        i_Q = MainFunctions.calc_iQ(S_Q, self.Sinf)
        F_r = MainFunctions.calc_Fr(self.r, self.Q[self.Q<=self.ui.QmaxInt.value()], i_Q[self.Q<=self.ui.QmaxInt.value()])
        
        Fopt_r, deltaFopt_r = Optimization.calc_optimize_Fr(self.ui.iterations.value(), self.F_r, \
            Fintra_r, density, i_Q[self.Q<=self.ui.QmaxInt.value()], \
            self.Q[self.Q<=self.ui.QmaxInt.value()], self.Sinf, \
            self.J_Q[self.Q<=self.ui.QmaxInt.value()], self.r, self.ui.rmin.value(), "n")
            
        Sopt_Q = MainFunctions.calc_SQCorr(Fopt_r, self.r, self.Q, self.Sinf)
            
        # # F_rIt, deltaF_rIt = Optimization.calc_optimize_Fr(self.ui.iterations.value(), self.F_r, \
            # # Fintra_r, self.ui.rho0Value.value(), self.i_Q[self.Q<=self.ui.QmaxInt.value()], \
            # # self.Q[self.Q<=self.ui.QmaxInt.value()], self.Sinf, \
            # # self.J_Q[self.Q<=self.ui.QmaxInt.value()], self.r, self.ui.rmin.value(), "n")
        
        self.ui.distfuncPlot.canvas.ax.plot(self.r, Fopt_r, label=r"F_{opt}(r)")
        self.ui.distfuncPlot.canvas.draw()
        
        self.ui.factorPlot.canvas.ax.plot(self.Q, Sopt_Q, label=r"S_{opt}(Q)")
        self.ui.factorPlot.canvas.draw()
    
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
