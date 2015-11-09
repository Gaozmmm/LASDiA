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
        self.button_load.clicked.connect(self.loadData)
        self.button_loadBkg = QtGui.QPushButton('Load Bkg')
        self.button_loadBkg.clicked.connect(self.loadBkg)

        self.button_calcSQ = QtGui.QPushButton('Calc S(Q)')
        self.button_calcSQ.clicked.connect(self.calcSQ)
        self.button_calcgr = QtGui.QPushButton('Calc g(r)')
        self.button_calcgr.clicked.connect(self.calcgr)

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
 
    def loadData(self):
        '''load and plot the file'''
        # Open data file
        path = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '.')
        file = open(path, "r")
        file_name = file.name
        ext = os.path.splitext(file_name)[1]
        
        if ext == ".chi":
            # Read and ignore the first 4 lines
            header1 = file.readline()
            header2 = file.readline()
            header3 = file.readline()
            header4 = file.readline()
        elif ext == ".xy":
            # Read and ignore the first 17 lines
            header1 = file.readline()
            header2 = file.readline()
            header3 = file.readline()
            header4 = file.readline()
            header5 = file.readline()
            header6 = file.readline()
            header7 = file.readline()
            header8 = file.readline()
            header9 = file.readline()
            header10 = file.readline()
            header11 = file.readline()
            header12 = file.readline()
            header13 = file.readline()
            header14 = file.readline()
            header15 = file.readline()
            header16 = file.readline()
            header17 = file.readline()


        # Read all the file in one time and close it
        lines = file.readlines()
        file.close()

        # Set the variables
        # transfMom: transfer momentum Q
        # intens: intensity I(Q)
        transfMomData = []
        intensData = []

        # Scan the rows of the file stored in lines, and put the values into the variables
        for line in lines:
            columns = line.split()
            transfMomData.append(float(columns[0]))
            intensData.append(float(columns[1]))

        # Modify the variables as numpy array
        self.transfMomVectData = np.array(transfMomData)
        self.intensVectData = np.array(intensData)

        # create an axis
        ax = self.figure.add_subplot(311)

        # discards the old graph
        ax.hold(False)

        # plot data
        ax.plot(self.transfMomVectData, self.intensVectData)

        # refresh canvas
        self.canvas.draw()
     
    #---------------------------------------------------------
    
    def loadBkg(self):
        '''load and plot the file'''
        # Open bkg file
        path = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '.')
        file = open(path, "r")
        file_name = file.name
        ext = os.path.splitext(file_name)[1]
        
        if ext == ".chi":
            # Read and ignore the first 4 lines
            header1 = file.readline()
            header2 = file.readline()
            header3 = file.readline()
            header4 = file.readline()
        elif ext == ".xy":
            # Read and ignore the first 17 lines
            header1 = file.readline()
            header2 = file.readline()
            header3 = file.readline()
            header4 = file.readline()
            header5 = file.readline()
            header6 = file.readline()
            header7 = file.readline()
            header8 = file.readline()
            header9 = file.readline()
            header10 = file.readline()
            header11 = file.readline()
            header12 = file.readline()
            header13 = file.readline()
            header14 = file.readline()
            header15 = file.readline()
            header16 = file.readline()
            header17 = file.readline()
        

        # Read all the file in one time and close it
        lines = file.readlines()
        file.close()

        # Set the variables
        # transfMom: transfer momentum Q
        # intens: intensity I(Q)
        transfMomBkg = []
        intensBkg = []

        # Scan the rows of the file stored in lines, and put the values into the variables
        for line in lines:
            columns = line.split()
            transfMomBkg.append(float(columns[0]))
            intensBkg.append(float(columns[1]))

        # Modify the variables as numpy array
        self.transfMomVectBkg = np.array(transfMomBkg)
        self.intensVectBkg = np.array(intensBkg)

        # create an axis
        ax = self.figure.add_subplot(311)

        # discards the old graph
        ax.hold(True)

        # plot data
        ax.plot(self.transfMomVectBkg, self.intensVectBkg, 'k--')

        # refresh canvas
        self.canvas.draw()
     
    #---------------------------------------------------------
        
    def calcSQ(self):
        # Function to calculate the S(Q)
        # Where S(Q) is: S(Q) = I_coh(Q) / (N * f^2(Q))

        # Set the parameters to calculate the atomic form factor for Argon
        # The atomic form factor formula and the parameters values are taken from:
        # http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
        a1 = 7.4845
        b1 = 0.9072
        a2 = 6.7723
        b2 = 14.8407
        a3 = 0.6539
        b3 = 43.8983
        a4 = 1.6442
        b4 = 33.3929
        c = 1.4445

        # a1 = 4.43918
        # b1 =1.64167
        # a2 = 3.20345
        # b2 = 3.43757
        # a3 = 1.19453
        # b3 = 0.2149
        # a4 = 0.41653
        # b4 = 6.65365
        # c = 0.746297

        
        # Calculate the atomic form factor
        # f(Q) = f1(Q) + f2(Q) + f3(Q) + f4(Q) + c
        # fi(Q) = ai * exp(-bi * (Q/4pi)^2)
        
        f1_Q = a1 * np.exp(-b1 * ( self.transfMomVectData / (4*np.pi)**2 ) )
        f2_Q = a2 * np.exp(-b2 * ( self.transfMomVectData / (4*np.pi)**2 ) )
        f3_Q = a3 * np.exp(-b3 * ( self.transfMomVectData / (4*np.pi)**2 ) )
        f4_Q = a4 * np.exp(-b4 * ( self.transfMomVectData / (4*np.pi)**2 ) )
        
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
        ax.plot(self.transfMomVectData, self.S_Q)

        # Refresh canvas
        self.canvas.draw()

    #---------------------------------------------------------
        
    def calcgr(self):

        #---------------------------------------------
        # Test with the integral
        r = np.linspace(0.5, 10, self.transfMomVectData.size)
        pi2 = 2.0 / np.pi
        rho_0 = 1.4
        i_Q = self.S_Q - 1
#        Q_max = np.amax(self.transfMomVect)

        g_r = pi2 * simps(self.transfMomVectData * i_Q * np.array(np.sin(np.mat(self.transfMomVectData).T * np.mat(r))).T, self.transfMomVectData)

        #g_r = 1 + f_r / (4 * np.pi * r * rho_0)
        
        #---------------------------------------------
        # # Test with the FFT
        
        # r = fftpack.fftfreq(self.transfMomVect.size)
        # f_r = fftpack.fft(self.S_Q)
        # f_r_abs = np.abs(f_r)
        
        # rho = 1.4
        # g_r_abs = 1 + f_r_abs * (4 * np.pi * r * rho)

        #---------------------------------------------
        # # Test with the sum:

        # r = np.linspace(0.5, 10, self.transfMomVectData.size)
        # pi2 = 2.0 / np.pi
        # rho = 1.4
        # i_Q = self.S_Q - 1
        # f_r = np.zeros(self.transfMomVectData.size)
        # for idx_r, ri in enumerate(r):
        #     for idx_Q, Qj in enumerate(self.transfMomVectData):
        #         f_r[idx_r] = pi2 * ( Qj * i_Q[idx_Q] * np.sin(Qj * ri))

        # g_r = 1 + f_r / (4 * np.pi * r * rho)
        
        #---------------------------------------------

        # create an axis
        ax = self.figure.add_subplot(313)

        # discards the old graph
        ax.hold(False)

        # plot data
        ax.plot(r, g_r)

        # refresh canvas
        self.canvas.draw()

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)

    main = Window()
    main.show()

    sys.exit(app.exec_())
