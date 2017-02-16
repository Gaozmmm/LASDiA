# The MIT License (MIT)

# Copyright (c) 2015-2016 European Synchrotron Radiation Facility

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""LASDiA main script file.
This script is mainly used for testing the software, but it can be used to run
LASDiA in text mode.

The nomenclature and the procedures follow the article:
    Eggert et al. 2002 PRB, 65, 174105.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by
an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""


from __future__ import (absolute_import, division, print_function, unicode_literals)

# import sys
import matplotlib.pyplot as plt
# import numpy as np
# import time
# from itertools import product
# from timeit import default_timer as timer
# from scipy.integrate import simps
# from scipy import fftpack
# from scipy import interpolate
# from scipy import signal
# import math

# from modules import Formalism
# from modules import Geometry
from modules import IgorFunctions
# from modules import KaplowMethod
from modules import MainFunctions
# from modules import Minimization
from modules import Optimization
from modules import Utility
from modules import UtilityAnalysis

# from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget


if __name__ == "__main__":

    # ---------------------------Files reading---------------------------------

    variables = Utility.read_inputFile("./inputFile.txt")

    elementList = Utility.molToelemList(variables.molecule)
    elementParameters = Utility.read_parameters(elementList, variables.element_params_path)

    path = Utility.path_xyz_file(variables.molecule)
    numAtoms, element, x, y, z = Utility.read_xyz_file(path)

    Q, I_Q = Utility.read_file(variables.data_file)
    Qbkg, Ibkg_Q = Utility.read_file(variables.bkg_file)

    # -------------------Preliminary calculation-------------------------------

    fe_Q, Ztot = MainFunctions.calc_eeff(elementList, Q, elementParameters)
    Iincoh_Q = MainFunctions.calc_Iincoh(elementList, Q, elementParameters)
    J_Q = IgorFunctions.calc_JQ(Iincoh_Q, fe_Q)
    Sinf = MainFunctions.calc_Sinf(elementList, fe_Q, Q, Ztot, elementParameters)

    dampingFunction = UtilityAnalysis.calc_dampingFunction(Q, variables.dampingFactor,
        variables.QmaxIntegrate, variables.typeFunction)

    # ------------------Intra-molecular components-----------------------------

    iintra_Q = Optimization.calc_iintra(Q, fe_Q, Ztot, variables.QmaxIntegrate,
        variables.maxQ, elementList, element, x, y, z, elementParameters)
    iintradamp_Q = UtilityAnalysis.calc_iintradamp(iintra_Q, Q, variables.QmaxIntegrate, 
        dampingFunction)
    rintra, Fintra_r = IgorFunctions.calc_FFT_QiQ(Q, iintradamp_Q, variables.QmaxIntegrate)

    _, dampingFunction = UtilityAnalysis.rebinning(Q, dampingFunction, 0.0,
        variables.maxQ, variables.NumPoints)

    # ---------------------Geometrical correction------------------------------

    absCorrFactor = IgorFunctions.absorption(Q)
    # I_Q = I_Q /(absCorrFactor)
    # Ibkg_Q  = Ibkg_Q / (absCorrFactor)

    # ------------------------Starting minimization----------------------------

    scaleFactor = variables.scaleFactor
    density = variables.density
    smoothingFactor = 0.2

    Subt = IgorFunctions.calc_Subt(I_Q, Ibkg_Q, scaleFactor, absCorrFactor)
    alpha = IgorFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, 
        Q[Q<=variables.QmaxIntegrate], Subt[Q<=variables.QmaxIntegrate],
        fe_Q[Q<=variables.QmaxIntegrate], Ztot, density)
    
    S_Q = IgorFunctions.calc_SQ(Q, Subt, alpha, fe_Q, J_Q, Ztot, Sinf, 
        variables.QmaxIntegrate)
    
    Utility.plot_data(Q, S_Q, "S_Q", "Q", "S_Q", "S_Q", "y")
    
    Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, 
        smoothingFactor, 
        variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    
    Utility.plot_data(Q, Ssmooth_Q, "S_Q", "Q", "S_Q", "S_Q", "y")
    
    # Q, Ssmooth_Q = UtilityAnalysis.rebinning(Q, Ssmooth_Q, 0.0,
        # variables.maxQ, variables.NumPoints)
    
    # # Q, Ssmooth_Q = Utility.read_file("../data/cea_files/ar/HT2_034SofQ_SS.txt")
    # # Utility.plot_data(Q, Ssmooth_Q, "I_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"igor", "y")
    
    # SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Sinf,
        # dampingFunction)
        
    # chi2 = IgorFunctions.FitRemoveGofRPeaks(Q, SsmoothDamp_Q, Sinf, 
        # variables.QmaxIntegrate, rintra, Fintra_r, variables.iterations, 
        # variables.rmin, density, J_Q, Ztot)

    # print(chi2)
    
    plt.show()