# The MIT License (MIT)

# Copyright (c) 2016 Francesco Devoto

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

from __future__ import (absolute_import, division, print_function, unicode_literals)
import six

from PySide import QtGui, QtCore

import sys
import os

import scipy.constants as sc
from scipy import fftpack
from scipy import signal
from scipy.integrate import simps
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use('Qt4Agg')
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np
import time
import math

from modules.MainFunctions import *
from modules.UtilityAnalysis import *
from modules.Utility import *
from modules.Optimization import *
from modules.Minimization import *
from modules.Formalism import *
from modules.IgorFunctions import *

from modules.gui.mplwidget import MplWidget
from modules.gui.gui import *

class LASDiA(QtGui.QMainWindow, Ui_LASDiAGui):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)

        # self.figure = plt.figure()
        # self.canvas = FigureCanvas(self.figure)
        #self.toolbar = NavigationToolbar(self.canvas, self)
        self.ui = Ui_LASDiAGui()
        self.ui.setupUi(self)
        # self.RawData = MplWidget()

        # Set the variable
        self.Q = None
        self.I_Q = None
        self.Qbkg = None
        self.I_Qbkg = None
        self.S_Q = None
        self.elementList = None
        self.molecule = None
        self.min_index = None
        self.max_index = None
        self.Ztot = None
        self.fe_Q = None
        self.r = None
        self.F_r = None
        self.Sinf = None
        self.i_Q = None
        self.rho0 = None
        self.Iincoh_Q = None
        self.J_Q = None
        self.minQ = None
        self.QmaxIntegrate = None
        self.maxQ = None
        self.validation_index = None
        self.integration_index = None
        self.calculation_index = None
        self.damp_factor = None
        self.newQ = None
        self.S_QsmoothedDamp = None

        # Set the buttons
        self.ui.LoadData.clicked.connect(self.load_data)
        self.ui.LoadBkg.clicked.connect(self.load_bkg)
        self.ui.CalcSQ.clicked.connect(self.calcSQ)
        # self.ui.CalcFr.clicked.connect(self.CalcFr)
        # self.ui.Optimization.clicked.connect(self.calcOptimization)
        #self.ui.Minimization.clicked.connect(self.calcMinimization)

        # Set plots

    #---------------------------------------------------------

    def load_data(self):
        '''load and plot the file'''
        # Open data file
        path = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '.')

        self.Q, self.I_Q = read_file(path[0])

    #---------------------------------------------------------

    def load_bkg(self):
        '''load and plot the file'''
        # Open bkg file
        path = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '.')

        # Modify the variables as numpy array
        self.Qbkg, self.I_Qbkg = read_file(path[0])


    #---------------------------------------------------------

    def calcSQ(self):
        self.molecule = str(self.ui.Molecule.currentText())
        # multiplicity = self.ui.Multiplicity.value()
        # self.elementList = {element:multiplicity}
        self.elementList = molToelemList(molecule)

        #numAtoms =  sc.N_A
        N = 1
        if self.ui.sValue.value() == 0.00:
            s = 1
        else:
            s = self.ui.sValue.value()

        if self.ui.rho0Value.value() == 0.00:
            self.rho0 = 25.0584
        else:
            self.rho0 = self.ui.rho0Value.value()

        self.minQ = self.ui.minQ.value()
        self.QmaxIntegrate = self.ui.QmaxIntegrate.value()
        self.maxQ = self.ui.maxQ.value()

        self.min_index, self.max_index = calc_indices(self.Q, self.minQ, self.QmaxIntegrate, self.maxQ)
        self.validation_index, self.integration_index, self.calculation_index = calc_ranges(self.Q, self.minQ, self.QmaxIntegrate, self.maxQ)

        self.fe_Q, self.Ztot = calc_eeff(self.elementList, self.Q)
        self.Iincoh_Q = calc_Iincoh(self.elementList, self.Q)
        self.J_Q = calc_JQ(self.Iincoh_Q, self.Ztot, self.fe_Q)
        self.Sinf = calc_Sinf(self.elementList, self.fe_Q, self.Q, self.Ztot)
        Isample_Q = calc_IsampleQ(self.I_Q, s, self.I_Qbkg)

        alpha = calc_alpha(self.J_Q[self.integration_index], self.Sinf, \
            self.Q[self.integration_index], Isample_Q[self.integration_index], \
            self.fe_Q[self.integration_index], self.Ztot, self.rho0)
        Icoh_Q = calc_Icoh(N, alpha, Isample_Q, self.Iincoh_Q)

        self.S_Q = calc_SQ(N, Icoh_Q, self.Ztot, self.fe_Q, self.Sinf, self.Q, self.max_index, self.integration_index)

        smooth_factor = self.ui.smoothFactor.value()
        damp_factor = self.ui.dampingFactor.value()

        self.newQ, S_Qsmoothed = calc_SQsmoothing(self.Q[self.validation_index], \
            self.S_Q[self.validation_index], self.Sinf, smooth_factor, \
            self.min_index, self.minQ, self.QmaxIntegrate, self.maxQ, 550)
        self.S_QsmoothedDamp = calc_SQdamp(S_Qsmoothed, self.newQ, self.Sinf, \
            self.QmaxIntegrate, self.damp_factor)

    #---------------------------------------------------------

    def CalcFr(self):
        Qi_Q = calc_QiQ(newQ, S_QsmoothedDamp, Sinf)
        i_Q = calc_iQ(S_QsmoothedDamp, Sinf)

        validation_indexSmooth, integration_indexSmooth, calculation_indexSmooth = calc_ranges(newQ, minQ, QmaxIntegrate, maxQ)
        min_indexSmooth, max_indexSmooth = calc_indices(newQ, minQ, QmaxIntegrate, maxQ)

        r = calc_r(newQ)
        F_r = calc_Fr(r, newQ[integration_indexSmooth], Qi_Q[integration_indexSmooth])

    #---------------------------------------------------------

    # def calcOptimization(self):
    #     iteration = self.ui.Iteration.value()
    #     r_cutoff = self.ui.rcutoff.value()
    #
    #     Fintra_r = calc_Fintra()
    #     optF_r = calc_optimize_Fr(iteration, self.F_r, Fintra_r, self.rho0, self.i_Q, self.Q, self.Sinf, self.J_Q, self.r, r_cutoff)
    #
    #     #self.ui.Fr.canvas.ax.clear()
    #     self.ui.Fr.canvas.ax.plot(self.r, optF_r)
    #     self.ui.Fr.canvas.draw()
    #
    #     # SQ_F = calc_SQ_F(optF_r, self.r, self.Q, self.Sinf)
    #
    #     # self.ui.SQ.canvas.ax.plot(self.Q, SQ_F)
    #     # self.ui.SQ.canvas.draw()

    #---------------------------------------------------------

#    def calcMinimization(self):


    #---------------------------------------------------------


def main():
    app = QtGui.QApplication(sys.argv)
    form = LASDiA()
    form.show()
    sys.exit(app.exec_())
    #app.exec_()



if __name__ == '__main__':
    main()
