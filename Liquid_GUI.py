# The MIT License (MIT)

# Copyright (c) 2015 Francesco Devoto

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

# First test of a GUI
# This script load a file, read the data inside and plot them on a graph
# Example with class and functions

# ToDo:
# The data saved in the files are Q and I(Q), modify calsSQ for the 2_theta
# Don't use A, it's not a formal part of SI, move to nm

# The nomenclature and the procedures follow the article: Eggert et al. 2002 PRB, 65, 174105.

from __future__ import (absolute_import, division, print_function, unicode_literals)
import six

import sys
from PyQt4 import QtGui

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import matplotlib.pyplot as plt

import numpy as np
import scipy.constants as sc
from scipy import fftpack
from scipy.integrate import simps

from modules.MainModules import *

import os

class Window(QtGui.QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        # Set the variable from the file
        # transfMomVect: numpy array of transfer momentum Q
        # intensVectData: numpy array of data intensity I(Q)
        # intensVectBkg: numpy array of bkg intensity I(Q)
        self.transfMomVectData = None
        self.intensVectData = None
        self.transfMomVectBkg = None
        self.intensVectBkg = None
        self.S_Q = None
                
        # Set the buttons
        self.button_load = QtGui.QPushButton('Load Data')
        self.button_load.clicked.connect(self.load_data)
        self.button_loadBkg = QtGui.QPushButton('Load Bkg')
        self.button_loadBkg.clicked.connect(self.load_bkg)

        self.button_calcSQ = QtGui.QPushButton('Calc S(Q)')
        self.button_calcSQ.clicked.connect(self.calc_sq)
        self.button_calcgr = QtGui.QPushButton('Calc g(r)')
        self.button_calcgr.clicked.connect(self.calc_gr)

        # Set the layout
        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.button_load)
        layout.addWidget(self.button_loadBkg)
        layout.addWidget(self.button_calcSQ)
        layout.addWidget(self.button_calcgr)
        self.setLayout(layout)

    #---------------------------------------------------------        
 
    def load_data(self):
        '''load and plot the file'''
        # Open data file
        path = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '.')
        
        self.transfMomVectData, self.intensVectData = read_file(path)

        # create an axis
        ax = self.figure.add_subplot(311)

        # discards the old graph
        ax.hold(False)

        # plot data
        ax.plot(self.transfMomVectData, self.intensVectData, label='Data')
        plt.legend()
        plt.xlabel('Q (1/A)')
        plt.ylabel('I(Q)')
        
        # refresh canvas
        self.canvas.draw()
     
    #---------------------------------------------------------
    
    def load_bkg(self):
        '''load and plot the file'''
        # Open bkg file
        path = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '.')
        
        # Modify the variables as numpy array
        self.transfMomVectBkg, self.intensVectBkg = read_file(path)

        # create an axis
        ax = self.figure.add_subplot(311)

        # discards the old graph
        ax.hold(True)

        # plot data
        ax.plot(self.transfMomVectBkg, self.intensVectBkg, 'k--', label='Bkg')
        plt.legend()

        # refresh canvas
        self.canvas.draw()
    
    #---------------------------------------------------------
    
    def calc_sq(self):
        """Function to calculate the S(Q)
        Where S(Q) is: S(Q) = I_coh(Q) / (N * f^2(Q))
        """
        
        # Read the parameters to calculate the atomic form factor for the element
        # The atomic form factor formula and the parameters values are taken from:
        # http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php

        element = "Ar"
        
        file = open("./affParam.txt", "r")
        header1 = file.readline()
        lines = file.readlines()
        file.close()

        for line in lines:
            columns = line.split()
            if columns[0] == element:
                a1 = float(columns[1])
                b1 = float(columns[2])
                a2 = float(columns[3])
                b2 = float(columns[4])
                a3 = float(columns[5])
                b3 = float(columns[6])
                a4 = float(columns[7])
                b4 = float(columns[8])
                c = float(columns[9])
                break
        
        # Calculate the atomic form factor
        # f(Q) = f1(Q) + f2(Q) + f3(Q) + f4(Q) + c
        # fi(Q) = ai * exp(-bi * (Q/4pi)^2)
        
        f1_Q = a1 * np.exp(-b1 * (self.transfMomVectData / (4*np.pi))**2)
        f2_Q = a2 * np.exp(-b2 * (self.transfMomVectData / (4*np.pi))**2)
        f3_Q = a3 * np.exp(-b3 * (self.transfMomVectData / (4*np.pi))**2)
        f4_Q = a4 * np.exp(-b4 * (self.transfMomVectData / (4*np.pi))**2)
        
        f_Q = f1_Q + f2_Q + f3_Q + f4_Q + c
        
        # Set for now the number of atoms as Avogadro constant
        numAtoms = sc.N_A
        
        # Calculate the structure factor S(Q)
        intensVect = self.intensVectData - self.intensVectBkg
        self.S_Q = intensVect / (numAtoms * f_Q**2)
        
        # Create an axis
        ax = self.figure.add_subplot(312)

        # Discards the old graph
        ax.hold(False)

        # Plot data
        ax.plot(self.transfMomVectData, self.S_Q, label='S(Q)')
        plt.legend()
        plt.xlabel('Q (1/A)')
        plt.ylabel('S(Q)')

        # Refresh canvas
        self.canvas.draw()

    #---------------------------------------------------------
        
    def calc_gr(self):

        #---------------------------------------------
        # Test with the integral
        r = np.linspace(0.5, 10, self.transfMomVectData.size)
        rho_0 = 1.4
        i_Q = self.S_Q - 1
#        Q_max = np.amax(self.transfMomVect)

        F_r = (2.0 / np.pi) * simps(self.transfMomVectData * i_Q * np.array(np.sin(np.mat(self.transfMomVectData).T * np.mat(r))).T, self.transfMomVectData)

        g_r = 1 + F_r / (4 * np.pi * r * rho_0)
        
        #---------------------------------------------
        # # Test with the FFT
        
        # r = fftpack.fftfreq(self.transfMomVect.size)
        # F_r = fftpack.fft(self.S_Q)
        # F_r_abs = np.abs(F_r)
        
        # rho = 1.4
        # g_r_abs = 1 + F_r_abs * (4 * np.pi * r * rho)

        #---------------------------------------------
        # # Test with the sum:

        # r = np.linspace(0.5, 10, self.transfMomVectData.size)
        # rho = 1.4
        # i_Q = self.S_Q - 1
        # F_r = np.zeros(self.transfMomVectData.size)
        # for idx_r, ri in enumerate(r):
        #     for idx_Q, Qj in enumerate(self.transfMomVectData):
        #         F_r[idx_r] = (2.0 / np.pi) * ( Qj * i_Q[idx_Q] * np.sin(Qj * ri))

        # g_r = 1 + F_r / (4 * np.pi * r * rho)
        
        #---------------------------------------------

        # create an axis
        ax = self.figure.add_subplot(313)

        # discards the old graph
        ax.hold(False)

        # plot data
        ax.plot(r, g_r, label='g(r)')
        plt.legend()
        plt.xlabel('r (A)')
        plt.ylabel('g(r)')

        # refresh canvas
        self.canvas.draw()

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)

    main = Window()
    main.show()

    sys.exit(app.exec_())
