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

    Q, I_Q, Qbkg, Ibkg_Q = UtilityAnalysis.check_data_length(Q, I_Q, Qbkg, Ibkg_Q, \
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
    
    scaleStep = 0.05
    scaleStepEnd = 0.00006
    densityStep = 0.025
    numSample = 23
    loopIteration = 0
    StopLoop = 1
    
    # plt.ion()
    # figure, ax = plt.subplots()

    while True: # Loop for the step changing
        NoPeak = 0

        # --------------------Scale minimization---------------------------
        Flag = 0
        while True: # Loop for the range shifting
            # ax.cla()
            # ax.grid(True)
            scaleArray = UtilityAnalysis.make_array_loop(scaleFactor, scaleStep, numSample)

            chi2Array = np.zeros(numSample)

            # plt.xlabel("Scale")
            # ax.relim()
            # ax.autoscale_view()
            for i in range(len(scaleArray)):
                
                # ------------------Kaplow method for scale--------------------

                Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scaleArray[i], Ibkg_Q)
                alpha = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, \
                    Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], \
                    fe_Q[Q<=variables.QmaxIntegrate], Ztot, density)
                Icoh_Q = MainFunctions.calc_Icoh(alpha, Isample_Q, Iincoh_Q)

                S_Q = MainFunctions.calc_SQ(Icoh_Q, Ztot, fe_Q, Sinf, Q, variables.minQ, \
                    variables.QmaxIntegrate, variables.maxQ)
                Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, 
                    variables.smoothingFactor, \
                    variables.minQ, variables.QmaxIntegrate, variables.maxQ)
                SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Sinf, \
                    dampingFunction)

                # Utility.plot_data(Q, SsmoothDamp_Q, "I_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$I(Q)$", "y")

                i_Q = MainFunctions.calc_iQ(SsmoothDamp_Q, Sinf)
                F_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], \
                    i_Q[Q<=variables.QmaxIntegrate])

                Fopt_r, deltaFopt_r = Optimization.calc_optimize_Fr(variables.iterations, F_r, \
                    Fintra_r, density, i_Q[Q<=variables.QmaxIntegrate], Q[Q<=variables.QmaxIntegrate], \
                    Sinf, J_Q[Q<=variables.QmaxIntegrate], r, variables.rmin, "n")

                chi2Array[i] = simps(deltaFopt_r[r < variables.rmin]**2, r[r < variables.rmin])
            
                # plt.scatter(scaleArray[i], chi2Array[i])
                # figure.canvas.draw()

                # ------------------End Kaplow method for scale----------------
            
            if np.amax(chi2Array) > 10**8:
                chi2Array = chi2Array[0:np.argmax(chi2Array)]

            scaleFactor = scaleArray[np.argmin(chi2Array)] - scaleStep*1.1
            
            nearIdx, nearEl = UtilityAnalysis.find_nearest(scaleArray, scaleFactor)
            
            # print("---------")
            # print("chi2Array ", chi2Array)
            # print("scaleArray ", scaleArray)
            # print("chi2 min ", np.amin(chi2Array))
            # print("scale min ", scaleArray[np.argmin(chi2Array)])
            # print("scaleFactor ", scaleFactor)

            if nearIdx == 0:
                scaleFactor-=scaleStep*10
                scaleStep*=10
                NoPeak+=1
            if nearIdx >= numSample-2:
                scaleFactor+=scaleStep*10
                scaleStep*=10
                NoPeak+=1

            scaleStep /= 10
            Flag += 1
            loopIteration += 1
            print(loopIteration, scaleFactor)
            if (10*scaleStep<=scaleStepEnd)*(NoPeak>=5)*((Flag!=1)+(scaleFactor+scaleStep*1.1<0)):
                break
        if (StopLoop == 1): 
            break
    
    
    # plt.ioff()
    # plt.show()
