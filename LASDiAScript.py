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
# from scipy.integrate import simps

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


if __name__ == '__main__':
    
    variables = Utility.read_inputFile("./inputFile.txt")
    
    elementList = Utility.molToelemList(variables.molecule)
    elementParameters = Utility.read_parameters(elementList, variables.element_params_path)
    
    path = Utility.path_xyz_file(variables.molecule)
    numAtoms, element, x, y, z = Utility.read_xyz_file(path)
    
    Q, I_Q = Utility.read_file(variables.data_file)
    Qbkg, Ibkg_Q  = Utility.read_file(variables.bkg_file)
    
    # I_Q = UtilityAnalysis.remove_peaks(Q, I_Q)
    
    lorch = UtilityAnalysis.calc_dampingFunction(Q, variables.dampFactor, variables.QmaxIntegrate, "Lorch Function")
    lorch2 = UtilityAnalysis.calc_dampingFunction(Q, variables.dampFactor, variables.maxQ, "Lorch Function")

    expon = UtilityAnalysis.calc_dampingFunction(Q, variables.dampFactor, variables.QmaxIntegrate, "Exponential")
    
    # S_Qdamp = UtilityAnalysis.calc_SQdamp(I_Q, 1, dampFunc)

    # S_Qdamp2 = UtilityAnalysis.calc_SQdamp2(I_Q, Q, 1, variables.QmaxIntegrate, variables.dampFactor)


    # Q, I_Q, Qbkg, Ibkg_Q = UtilityAnalysis.check_data_length(Q, I_Q, Qbkg, Ibkg_Q, \
    #     variables.minQ, variables.maxQ)
    
    # fe_Q, Ztot = MainFunctions.calc_eeff(elementList, Q, elementParameters)
    # Iincoh_Q = MainFunctions.calc_Iincoh(elementList, Q, elementParameters)
    # J_Q = MainFunctions.calc_JQ(Iincoh_Q, Ztot, fe_Q)
    # Sinf, Sinf_Q = MainFunctions.calc_Sinf(elementList, fe_Q, Q, Ztot, elementParameters)
    
    # r, iintradamp_Q, Fintra_r = Optimization.calc_intraComponent(Q, fe_Q, Ztot, \
    #     variables.QmaxIntegrate, variables.maxQ, elementList, element, \
    #     x, y, z, elementParameters, variables.dampFactor)
    
    # scaleFactor = variables.sfValue
    # density = variables.rho0Value
    
    # density, scaleFactor = Minimization.chi2_minimization(scaleFactor, Q, I_Q, Ibkg_Q, 
    #     J_Q, fe_Q, Iincoh_Q, Sinf, Ztot,
    #     density, Fintra_r, r, variables.minQ, variables.QmaxIntegrate, variables.maxQ, 
    #     variables.smoothFactor, variables.dampFactor, variables.iterations, variables.rmin)
    
    # print("Final values ", density, scaleFactor)
    
    # S_Q = UtilityAnalysis.S_QCalculation(Q, I_Q, Ibkg_Q, scaleFactor, J_Q, Sinf, fe_Q, Ztot, density, Iincoh_Q, 
    #     variables.minQ, variables.QmaxIntegrate, variables.maxQ, variables.smoothFactor, variables.dampFactor)
    
    # i_Q = MainFunctions.calc_iQ(S_Q, Sinf)
    # F_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], i_Q[Q<=variables.QmaxIntegrate])
    
    # Fopt_r, deltaFopt_r = Optimization.calc_optimize_Fr(variables.iterations, F_r, 
    #         Fintra_r, density, i_Q[Q<=variables.QmaxIntegrate],
    #         Q[Q<=variables.QmaxIntegrate], Sinf,
    #         J_Q[Q<=variables.QmaxIntegrate], r, variables.rmin, "n")
        
    # Sopt_Q = MainFunctions.calc_SQCorr(Fopt_r, r, Q, Sinf)
    
    # Utility.plot_data(Q, S_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S(Q)$", "y")
    # Utility.plot_data(r, F_r, "F_r", r"$r(nm)$", r"$F(r)$", r"$F(r)$", "y")
    # Utility.plot_data(r, Fopt_r, "F_r", r"$r(nm)$", r"$F(r)$", r"$F_{opt}(r)$", "y")
    # Utility.plot_data(Q, Sopt_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S_{opt}(Q)$", "y")
    # Utility.plot_data(Q, I_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$I(Q)$", "y")
    # Utility.plot_data(Q, S_Qdamp, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$I(Q)damp$", "y")
    # Utility.plot_data(Q, S_Qdamp2, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$I(Q)damp2$", "y")

    Utility.plot_data(Q, lorch, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$lorch$", "y")
    Utility.plot_data(Q, lorch2, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$lorch2$", "y")

    # Utility.plot_data(Q, expon, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$exp$", "y")
   

    plt.show()
