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
    rintra, Fintra_r = UtilityAnalysis.calc_FFT_QiQ(Q, iintradamp_Q, variables.QmaxIntegrate)
    
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

    Q, SsmoothDamp_Q = UtilityAnalysis.rebinning(Q, SsmoothDamp_Q, 0.0, 
        variables.maxQ, variables.NumPoints)
    
    i_Q = MainFunctions.calc_iQ(SsmoothDamp_Q, Sinf)
    Qi_Q = Q*i_Q
    r, GR = UtilityAnalysis.calc_FFT_QiQ(Q, Qi_Q, variables.QmaxIntegrate)
    _, Fintra_r = UtilityAnalysis.rebinning(rintra, Fintra_r, np.amin(Fintra_r), 
        np.amax(Fintra_r), len(GR))
    
    # ------------------------------Qi(Q) with IFFT----------------------------
    
    QiQ1 = np.zeros(len(SsmoothDamp_Q))
    idx, _ = UtilityAnalysis.find_nearest(Qi_Q, variables.QmaxIntegrate)
    QiQ1[Q<variables.QmaxIntegrate] = Qi_Q[Q<variables.QmaxIntegrate]
    QiQ1[0] = 0.0
    
    GR1 = GR
    DelG = np.zeros(len(GR))
    Rnn = variables.rmin
    DelG[r<Rnn] = GR1[r<Rnn]-(Fintra_r[r<Rnn]-4*np.pi*r[r<Rnn]*density)
    
    # Utility.plot_data(r, DelG, "GR", "r", "G(r)", "G(r)", "y")
    
    for i in range(variables.iterations):
        Q1, QiQCorr = UtilityAnalysis.calc_IFFT_Fr(r, DelG)
        mask = np.where((Q1>0.0) & (Q1<variables.QmaxIntegrate))
        QiQ1[mask] = QiQ1[mask] - (QiQ1[mask] / 
            (Q1[mask] *(Sinf + J_Q[:len(Q1[mask])])) + 1) * QiQCorr[mask]
        r, GR1 = UtilityAnalysis.calc_FFT_QiQ(Q1, QiQ1, variables.QmaxIntegrate)
        
        DelG = np.zeros(len(GR1))
        DelG[r<Rnn] = GR1[r<Rnn]-(Fintra_r[r<Rnn]-4*np.pi*r[r<Rnn]*density)
        
        Utility.plot_data(r, GR1, "GR", "r", "G(r)", "G(r)", "y")
        plt.show()
        
        Rnn = 0.99*r[np.where(GR1==np.amin(GR1[r>0.95*Rnn]))[0][0]]
        print(Rnn)

    SQCorr = np.zeros(len(QiQ1))
    SQCorr[1:] = QiQ1[1:]/Q1[1:]+Sinf
    
    DTemp = DelG
    DTemp = DTemp**2
    print(np.mean(DTemp))
    
    # Utility.plot_data(r, DTemp, "GR", "r", "G(r)", "G(r)", "y")
    
    plt.show()