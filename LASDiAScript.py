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
    
    Q, I_Q, Qbkg, Ibkg_Q = UtilityAnalysis.check_data_length(Q, I_Q, Qbkg, Ibkg_Q, \
        variables.minQ, variables.maxQ)
    
    fe_Q, Ztot = MainFunctions.calc_eeff(elementList, Q, elementParameters)
    Iincoh_Q = MainFunctions.calc_Iincoh(elementList, Q, elementParameters)
    J_Q = MainFunctions.calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf, Sinf_Q = MainFunctions.calc_Sinf(elementList, fe_Q, Q, Ztot, elementParameters)
    
    r, iintradamp_Q, Fintra_r = Optimization.calc_intraComponent(Q, fe_Q, Ztot, \
        variables.QmaxIntegrate, variables.maxQ, elementList, element, \
        x, y, z, elementParameters, variables.dampFactor)
    
    scaleFactor = variables.sfValue
    density = variables.rho0Value
    
    scaleStep = 0.05
    densityStep = 0.025
    numSample = 23
    numLoopIteration = 0

    
    plt.ion()
    figure, ax = plt.subplots()

    while True:
        plt.cla()
        ax.grid()
        scaleArray = UtilityAnalysis.make_array_loop(scaleFactor, scaleStep, numSample)

        chi2Array = np.zeros(numSample)

        plt.xlabel("Scale")
        for i in range(len(scaleArray)):
            chi2Array[i], SsmoothDamp_Q, F_r, Fopt_r = KaplowMethod.Kaplow_method(variables, Q, I_Q, \
                Ibkg_Q, J_Q, fe_Q, Iincoh_Q, Sinf, Ztot, scaleArray[i], density, Fintra_r, r)
            
            # plt.figure("chi2")
            plt.scatter(scaleArray[i], chi2Array[i])
            figure.canvas.draw()
        
        xfit, yfit, scaleFactor = Minimization.chi2_fit(scaleArray, chi2Array)
        plt.plot(xfit, yfit)
        figure.canvas.draw()
        print(scaleFactor)
        
        # time.sleep(1.0)
        plt.cla()
        ax.grid()

        density0 = density
        densityArray = UtilityAnalysis.make_array_loop(density, densityStep, numSample)
        chi2Array = np.zeros(numSample)
        
        plt.xlabel("Density")
        for i in range(len(densityArray)):
            chi2Array[i], SsmoothDamp_Q, F_r, Fopt_r = KaplowMethod.Kaplow_method(variables, Q, I_Q, \
                Ibkg_Q, J_Q, fe_Q, Iincoh_Q, Sinf, Ztot, scaleFactor, densityArray[i], Fintra_r, r)
            
            plt.scatter(densityArray[i], chi2Array[i])
            figure.canvas.draw()
            
        
        xfit, yfit, density = Minimization.chi2_fit(densityArray, chi2Array)
        plt.plot(xfit, yfit)
        figure.canvas.draw()
        print(density0, density)
        # time.sleep(1.0)

        if np.abs(density-density0) > density0/25:
            scaleStep = 0.006
            densityStep = density0/10
            print(1, np.abs(density-density0), density0/25)
        elif np.abs(density-density0) > density0/75:
            scaleStep = 0.0006
            densityStep = density0/100
            print(2, np.abs(density-density0), density0/75)
        else:
            scaleStep = 0.00006
            densityStep = density0/1000
            print(3, np.abs(density-density0))

        numLoopIteration += 1
        print(numLoopIteration)
        # if (numLoopIteration == 2):
        if (np.abs(density-density0) < density0/2500 or numLoopIteration > 30):
            print(4, np.abs(density-density0), density0/2500)
            break

    plt.ioff()
        
    # # Utility.plot_data(Q, SsmoothDamp_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S(Q)$", "y")
    # # Utility.plot_data(r, F_r, "F_r", r"$r(nm)$", r"$F(r)$", r"$F(r)$", "y")
    # # Utility.plot_data(r, Fopt_r, "F_r", r"$r(nm)$", r"$F(r)$", r"$F_{opt}(r)$", "y")
    
    plt.show()
