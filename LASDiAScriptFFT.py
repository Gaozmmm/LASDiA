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

"""LASDiA main script file.
This script is mainly used for testing the software, but it can be used to run LASDiA
in text mode.

The nomenclature and the procedures follow the article: Eggert et al. 2002 PRB, 65, 174105.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by
an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""


from __future__ import (absolute_import, division, print_function, unicode_literals)
import six

import sys
import matplotlib.pyplot as plt
import numpy as np
import time
# from itertools import product
# from timeit import default_timer as timer
from scipy.integrate import simps
from scipy import fftpack
import math

from modules import Formalism
from modules import Geometry
from modules import IgorFunctions
from modules import KaplowMethod
from modules import MainFunctions
from modules import Minimization
from modules import Optimization
from modules import Utility
from modules import UtilityAnalysis

from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget


if __name__ == "__main__":
    
    #---------------------------Files reading----------------------------------
    
    variables = Utility.read_inputFile("./inputFile.txt")
    
    elementList = Utility.molToelemList(variables.molecule)
    elementParameters = Utility.read_parameters(elementList, variables.element_params_path)
    
    path = Utility.path_xyz_file(variables.molecule)
    numAtoms, element, x, y, z = Utility.read_xyz_file(path)
    
    Q, I_Q = Utility.read_file(variables.data_file)
    Qbkg, Ibkg_Q  = Utility.read_file(variables.bkg_file)
    
    
    #--------------------Preliminary calculation-------------------------------

    Q, I_Q, Qbkg, Ibkg_Q = UtilityAnalysis.check_data_length(Q, I_Q, Qbkg, Ibkg_Q,
        variables.minQ, variables.maxQ)
        
    
    fe_Q, Ztot = MainFunctions.calc_eeff(elementList, Q, elementParameters)
    Iincoh_Q = MainFunctions.calc_Iincoh(elementList, Q, elementParameters)
    J_Q = MainFunctions.calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf, Sinf_Q = MainFunctions.calc_Sinf(elementList, fe_Q, Q, Ztot, elementParameters)
    
    dampingFunction = UtilityAnalysis.calc_dampingFunction(Q, variables.dampingFactor,
        variables.QmaxIntegrate, variables.typeFunction)

    #-------------------Intra-molecular components-----------------------------

    iintra_Q = Optimization.calc_iintra(Q, fe_Q, Ztot, variables.QmaxIntegrate, 
        variables.maxQ, elementList, element, x, y, z, elementParameters)
    iintradamp_Q = UtilityAnalysis.calc_iintradamp(iintra_Q, Q, variables.QmaxIntegrate, 
        dampingFunction)
    r = MainFunctions.calc_r(Q)
    Fintra_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], 
        iintradamp_Q[Q<=variables.QmaxIntegrate])

    # ------------------------Starting minimization----------------------------

    scaleFactor = variables.scaleFactor
    density = variables.density
    
    Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scaleFactor, Ibkg_Q)
    alpha = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, 
        Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate],
        fe_Q[Q<=variables.QmaxIntegrate], Ztot, density)
    Icoh_Q = MainFunctions.calc_Icoh(alpha, Isample_Q, Iincoh_Q)

    S_Q = MainFunctions.calc_SQ(Icoh_Q, Ztot, fe_Q, Sinf, Q, variables.minQ, 
        variables.QmaxIntegrate, variables.maxQ)
    Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, 
        variables.smoothingFactor, 
        variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Sinf,
        dampingFunction)
    
    # Utility.plot_data(Q, SsmoothDamp_Q, "S_Q", "Q", "S(Q)", "S(Q)", "y")
    
    i_Q = MainFunctions.calc_iQ(SsmoothDamp_Q, Sinf)
    F_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], 
        i_Q[Q<=variables.QmaxIntegrate])
    
    Num_SS =550
    
    Q, i_Q = UtilityAnalysis.rebinning(Q, i_Q, 0.0, 
        variables.maxQ, Num_SS)
    
    # Utility.plot_data(Q, SsmoothDamp_Q, "S_Q", "Q", "S(Q)", "S(Q)", "y")
    # Utility.plot_data(r, F_r, "F_r", "r", "F(r)", "F(r)", "n")
    
    # ----------------------------FFT calculation------------------------------
    
    pMax, elem = UtilityAnalysis.find_nearest(Q, variables.QmaxIntegrate)
    NumPoints = 2*2*2**math.ceil(math.log(5*(pMax+1))/math.log(2))
    DelR = 2*np.pi/(np.mean(np.diff(Q))*NumPoints)
    Qi_Q = np.resize(Q*i_Q, NumPoints)
    Q = np.linspace(np.amin(Q), np.amax(Q), len(Qi_Q), endpoint=True)
    Qi_Q[pMax+1:] = 0.0
    DeltaQ = np.diff(Q)
    meanDeltaQ = np.mean(DeltaQ)
    r2 = fftpack.fftfreq(Q.size, meanDeltaQ)
    F_r2 = fftpack.fft(Qi_Q)
    F_r2 = F_r2[np.where(r2>=0)]
    F_r2 = -np.imag(F_r2)*meanDeltaQ*2/np.pi
    r2 = np.arange(0, 0+DelR*len(F_r2), DelR)
    
    # Utility.plot_data(r, F_r, "F_r", "r", "F(r)", "F(r)", "y")
    # Utility.plot_data(r2, F_r2, "F_r2", "r", "F(r)", "F2(r)", "y")
    
    # ----------------------------IFFT calculation-----------------------------
    
    # print(len(F_r2))
    
    NumPoints = 2**math.ceil(math.log(len(F_r2)-1)/math.log(2))
    for n in range(len(F_r2), NumPoints):
        F_r2.append(0)
    
    # print(len(F_r2))
    
    # Utility.plot_data(r2, F_r2, "F_r2", "r", "F(r)", "F2(r)", "y")
    
    DelG = np.zeros(len(F_r2))
    print(len(DelG))
    # Utility.plot_data(r2, DelG, "F_r2", "r", "F(r)", "F2(r)", "y")
    DelG[r2<variables.rmin] = F_r2[r2<variables.rmin]-4*np.pi*r2[r2<variables.rmin]*density
    # plt.plot(DelG)
    # Utility.plot_data(r2, DelG, "F_r2", "r", "F(r)", "F2(r)", "y")
    # print(len(DelG))
    
    F_r2 = DelG
    
    Q2 = np.linspace(0.0, variables.maxQ, Num_SS, endpoint=True)
    DelQ = 2*np.pi/(np.mean(np.diff(r2))*NumPoints)
    print(DelQ)
    
    Deltar = np.diff(r2)
    meanDeltar = np.mean(Deltar)
    Q2 = fftpack.fftfreq(r2.size, meanDeltar)
    QiQ = fftpack.fft(F_r2)
    QiQ = QiQ[np.where(Q2>=0)]
    print(len(QiQ))
    QiQ = -np.imag(QiQ)*meanDeltar
    # Q2 = np.arange(0, 0+DelQ*Num_SS, DelQ)
    
    # Utility.plot_data(Q2, QiQ, "QiQ", "r", "F(r)", "F2(r)", "y")
    
    plt.show()