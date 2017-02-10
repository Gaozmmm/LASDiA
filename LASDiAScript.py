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
    
    # Q, I_Q, Qbkg, Ibkg_Q = UtilityAnalysis.check_data_length(Q, I_Q, Qbkg, Ibkg_Q, \
        # variables.minQ, variables.maxQ)
    
    fe_Q, Ztot = MainFunctions.calc_eeff(elementList, Q, elementParameters)
    Iincoh_Q = MainFunctions.calc_Iincoh(elementList, Q, elementParameters)
    J_Q = IgorFunctions.calc_JQ(Iincoh_Q, fe_Q)
    Sinf = MainFunctions.calc_Sinf(elementList, fe_Q, Q, Ztot, elementParameters)
    
    dampingFunction = UtilityAnalysis.calc_dampingFunction(Q, variables.dampingFactor,
        variables.QmaxIntegrate, variables.typeFunction)

    #-------------------Intra-molecular components-----------------------------

    iintra_Q = Optimization.calc_iintra(Q, fe_Q, Ztot, variables.QmaxIntegrate, 
        variables.maxQ, elementList, element, x, y, z, elementParameters)
    iintradamp_Q = UtilityAnalysis.calc_iintradamp(iintra_Q, Q, variables.QmaxIntegrate, 
        dampingFunction)
    rintra, Fintra_r = IgorFunctions.calc_FFT_QiQ(Q, iintradamp_Q, variables.QmaxIntegrate)
    
    _, dampingFunction = UtilityAnalysis.rebinning(Q, dampingFunction, 0.0, 
        variables.maxQ, variables.NumPoints)
    
    Qorg = Q
    
    # ---------------------Geometrical correction------------------------------
    
    absCorrFactor = IgorFunctions.absorption(Q)
    # I_Q = I_Q /(absCorrFactor)
    # Ibkg_Q  = Ibkg_Q / (absCorrFactor)
    
    # ------------------------Starting minimization----------------------------

    scaleFactor = variables.scaleFactor
    density = variables.density
    
    scaleStep = 0.05
    scaleStepEnd = 0.00006
    densityStep = density/50
    densityStepEnd = density/250
    numSample = 23
    loopIteration = 0
    NoPeak = 0
    scaleFactor = scaleFactor-scaleStep*11
    
    # ----------------------First scale minimization---------------------------
    Flag = 0
    while True: # Loop for the range shifting
        
        # print("------")
        # print(Flag, scaleFactor, scaleStep, scaleStepEnd)
        # scaleArray = np.zeros(numSample)
        scaleArray = UtilityAnalysis.make_array_loop(scaleFactor, scaleStep, numSample)
        # print(scaleArray)
        chi2Array = np.zeros(numSample)
        
        for i in range(len(scaleArray)):
            
            # ------------------Kaplow method for scale--------------------
            Q = Qorg
            Subt = IgorFunctions.calc_Subt(I_Q, Ibkg_Q, scaleArray[i], absCorrFactor)
            alpha = IgorFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, 
                Q[Q<=variables.QmaxIntegrate], Subt[Q<=variables.QmaxIntegrate],
                fe_Q[Q<=variables.QmaxIntegrate], Ztot, density)
            
            S_Q = IgorFunctions.calc_SQ(Q, Subt, alpha, fe_Q, J_Q, Ztot, Sinf, 
                variables.QmaxIntegrate)
                
            Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, 
                variables.smoothingFactor, 
                variables.minQ, variables.QmaxIntegrate, variables.maxQ)
            Q, Ssmooth_Q = UtilityAnalysis.rebinning(Q, Ssmooth_Q, 0.0, 
                variables.maxQ, variables.NumPoints)
            SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Sinf,
                dampingFunction)
            
            Q, SsmoothDamp_Q = UtilityAnalysis.rebinning(Q, SsmoothDamp_Q, 0.0, 
                variables.maxQ, variables.NumPoints)
            
            chi2Array[i] = IgorFunctions.FitRemoveGofRPeaks(Q, SsmoothDamp_Q, Sinf, 
                variables.QmaxIntegrate, rintra, Fintra_r, variables.iterations, 
                variables.rmin, density, J_Q, Ztot)
            # print(chi2Array[i])
        
        # --------------------Range shifting selection --------------------
        
        # Utility.write_file("./chi2Array"+str(Flag)+".txt", scaleArray, chi2Array)
        
        if np.amax(chi2Array) > 10**8:
            scaleFactor = scaleArray[np.argmin(chi2Array[0:np.argmax(chi2Array)])] # - scaleStep*1.1
        else:
            scaleFactor = scaleArray[np.argmin(chi2Array)] #- scaleStep*1.1
        
        print(Flag, scaleFactor, scaleStep, scaleStepEnd)
        
        nearIdx, nearEl = UtilityAnalysis.find_nearest(scaleArray, scaleFactor)
        
        if nearIdx == 0:
            scaleFactor -= scaleStep*10
            scaleStep *= 10
            NoPeak += 1
        if nearIdx >= numSample-2:
            scaleFactor += scaleStep*10
            scaleStep *= 10
            NoPeak += 1

        scaleStep /= 10
        Flag += 1

        
        # print("------")
        #if (10*scaleStep<=scaleStepEnd) and (NoPeak>=5) and ((Flag!=1) or (scaleFactor+scaleStep*1.1<0)):
        if (Flag==2):
            break
        
    # ------------------------chi2 curve fit for scale-------------------------
    
    print(scaleFactor)
    
    # print("finale scale array ", scaleArray)
    
    # if scaleFactor < 0:
        # print("Scale factor < 0")
        # # break
    # else:
        # left = scaleArray[0]
        # right = scaleArray[-1]
        # coeffs = np.polyfit(scaleArray, chi2Array, 3)
        # if (4*coeffs[1]**2 - 12*coeffs[2]*coeffs[0] < 0):
            # scaleFactor = (left+right)/2
        # else:
            # x1=(-2*coeffs[1]+(4*coeffs[1]**2-12*coeffs[2]*coeffs[0])**0.5)/(6*coeffs[0])
            # x2=(-2*coeffs[1]-(4*coeffs[1]**2-12*coeffs[2]*coeffs[0])**0.5)/(6*coeffs[0])
            # if (2*coeffs[1]+6*coeffs[0]*x1 > 2*coeffs[1]+6*coeffs[0]*x2):
                # scaleFactor = left + np.diff(scaleArray)[0]*x1
            # else:
                # scaleFactor = left + np.diff(scaleArray)[0]*x2
    
    # plt.ion()
    # figure, ax = plt.subplots()
    # ax.cla()
    # ax.grid(True)
    # plt.xlabel("Scale")
    # plt.plot(scaleArray, chi2Array, "o")
    # pol_fit = np.poly1d(np.polyfit(scaleArray, chi2Array, 3))
    # x_fit = np.linspace(scaleArray[0], scaleArray[-1], 1000)
    # y_fit = pol_fit(x_fit)
    # plt.plot(x_fit, y_fit)
    # figure.canvas.draw()

    # print("finale scale factor", scaleFactor)
    
    # # ----------------------First density minimization-------------------------
    
    # NoPeak = 0
    # while True: # Loop for the range shifting
        
        # density0 = density
        # densityArray = UtilityAnalysis.make_array_loop(density, densityStep, numSample)
        # chi2Array = np.zeros(numSample)
        
        # for i in range(len(densityArray)):
            
            # # ----------------Kaplow method for density--------------------

            # Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scaleFactor, Ibkg_Q)
            # alpha = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, \
                # Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], \
                # fe_Q[Q<=variables.QmaxIntegrate], Ztot, densityArray[i])
            # Icoh_Q = MainFunctions.calc_Icoh(alpha, Isample_Q, Iincoh_Q)

            # S_Q = MainFunctions.calc_SQ(Icoh_Q, Ztot, fe_Q, Sinf, Q, variables.minQ, \
                # variables.QmaxIntegrate, variables.maxQ)
            # Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, 
                # variables.smoothingFactor, \
                # variables.minQ, variables.QmaxIntegrate, variables.maxQ)
            # SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Sinf, \
                # dampingFunction)

            # # Utility.plot_data(Q, SsmoothDamp_Q, "I_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$I(Q)$", "y")

            # i_Q = MainFunctions.calc_iQ(SsmoothDamp_Q, Sinf)
            # F_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], \
                # i_Q[Q<=variables.QmaxIntegrate])

            # Fopt_r, deltaFopt_r = Optimization.calc_optimize_Fr(variables.iterations, F_r, \
                # Fintra_r, densityArray[i], i_Q[Q<=variables.QmaxIntegrate], Q[Q<=variables.QmaxIntegrate], \
                # Sinf, J_Q[Q<=variables.QmaxIntegrate], r, variables.rmin, "n")

            # chi2Array[i] = simps(deltaFopt_r[r < variables.rmin]**2, r[r < variables.rmin])
        
        # # --------------------Range shifting selection---------------------
        
        # density = densityArray[np.argmin(chi2Array)] - scaleStep*1.1
        
        # nearIdx, nearEl = UtilityAnalysis.find_nearest(densityArray, density)
        
        # if nearIdx == 0:
            # density-=scaleStep*10
            # scaleStep*=10
            # NoPeak+=1
        # if nearIdx >= numSample-2:
            # density+=scaleStep*10
            # scaleStep*=10
            # NoPeak+=1

        # scaleStep /= 10

        # # if (10*densityStep<=densityStepEnd)*(NoPeak>=5):
        # if (1==1):
            # break
        
    # # -------------------chi2 curve fit for density--------------------
    
    # left = densityArray[0]
    # right = densityArray[-1]
    # coeffs = np.polyfit(densityArray, chi2Array, 3)
    # if (4*coeffs[1]**2 - 12*coeffs[2]*coeffs[0] < 0):
            # density = (left+right)/2
    # else:
        # x1=(-2*coeffs[1]+(4*coeffs[1]**2-12*coeffs[2]*coeffs[0])**0.5)/(6*coeffs[0])
        # x2=(-2*coeffs[1]-(4*coeffs[1]**2-12*coeffs[2]*coeffs[0])**0.5)/(6*coeffs[0])
        # if (2*coeffs[1]+6*coeffs[0]*x1 > 2*coeffs[1]+6*coeffs[0]*x2):
            # density = left + np.diff(densityArray)[0]*x1
        # else:
            # density = left + np.diff(densityArray)[0]*x2

    # ax.cla()
    # ax.grid(True)
    # plt.xlabel("Density")
    # # ax.relim()
    # # ax.autoscale_view()
    # plt.plot(densityArray, chi2Array, "o")
    # pol_fit = np.poly1d(np.polyfit(densityArray, chi2Array, 3))
    # x_fit = np.linspace(densityArray[0], densityArray[-1], 1000)
    # y_fit = pol_fit(x_fit)
    # plt.plot(x_fit, y_fit)
    # figure.canvas.draw()

    
    
    
    # while True: # Loop for the step changing
        # NoPeak = 0

        # # --------------------Scale minimization---------------------------
        # Flag = 0
        # while True: # Loop for the range shifting
            
            # scaleArray = UtilityAnalysis.make_array_loop(scaleFactor, scaleStep, numSample)
            # chi2Array = np.zeros(numSample)

            
            # for i in range(len(scaleArray)):
                
                # # ------------------Kaplow method for scale--------------------

                # Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scaleArray[i], Ibkg_Q)
                # alpha = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, \
                    # Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], \
                    # fe_Q[Q<=variables.QmaxIntegrate], Ztot, density)
                # Icoh_Q = MainFunctions.calc_Icoh(alpha, Isample_Q, Iincoh_Q)

                # S_Q = MainFunctions.calc_SQ(Icoh_Q, Ztot, fe_Q, Sinf, Q, variables.minQ, \
                    # variables.QmaxIntegrate, variables.maxQ)
                # Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, 
                    # variables.smoothingFactor, \
                    # variables.minQ, variables.QmaxIntegrate, variables.maxQ)
                # SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Sinf, \
                    # dampingFunction)

                # # Utility.plot_data(Q, SsmoothDamp_Q, "I_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$I(Q)$", "y")

                # i_Q = MainFunctions.calc_iQ(SsmoothDamp_Q, Sinf)
                # F_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], \
                    # i_Q[Q<=variables.QmaxIntegrate])

                # Fopt_r, deltaFopt_r = Optimization.calc_optimize_Fr(variables.iterations, F_r, \
                    # Fintra_r, density, i_Q[Q<=variables.QmaxIntegrate], Q[Q<=variables.QmaxIntegrate], \
                    # Sinf, J_Q[Q<=variables.QmaxIntegrate], r, variables.rmin, "n")

                # chi2Array[i] = simps(deltaFopt_r[r < variables.rmin]**2, r[r < variables.rmin])
            
            # # --------------------Range shifting selection --------------------
            
            # if np.amax(chi2Array) > 10**8:
                # scaleFactor = scaleArray[np.argmin(chi2Array[0:np.argmax(chi2Array)])] - scaleStep*1.1
            # else:
                # scaleFactor = scaleArray[np.argmin(chi2Array)] - scaleStep*1.1
            
            # nearIdx, nearEl = UtilityAnalysis.find_nearest(scaleArray, scaleFactor)
            
            # if nearIdx == 0:
                # scaleFactor-=scaleStep*10
                # scaleStep*=10
                # NoPeak+=1
            # if nearIdx >= numSample-2:
                # scaleFactor+=scaleStep*10
                # scaleStep*=10
                # NoPeak+=1

            # scaleStep /= 10
            # Flag += 1

            # if (10*scaleStep<=scaleStepEnd)*(NoPeak>=5)*((Flag!=1)+(scaleFactor+scaleStep*1.1<0)):
                # break
            
            # # ---------------------chi2 curve fit for scale--------------------
            
            # if scaleFactor < 0:
                # print("Scale factor < 0")
                # break
            # else:
                # left = scaleArray[0]
                # right = scaleArray[-1]
                # coeffs = np.polyfit(scaleArray, chi2Array, 3)
                # if (4*coeffs[1]**2 - 12*coeffs[2]*coeffs[0] < 0):
                    # scaleFactor = (left+right)/2
                # else:
                    # x1=(-2*coeffs[1]+(4*coeffs[1]**2-12*coeffs[2]*coeffs[0])**0.5)/(6*coeffs[0])
                    # x2=(-2*coeffs[1]-(4*coeffs[1]**2-12*coeffs[2]*coeffs[0])**0.5)/(6*coeffs[0])
                    # if (2*coeffs[1]+6*coeffs[0]*x1 > 2*coeffs[1]+6*coeffs[0]*x2):
                        # scaleFactor = left + np.diff(scaleArray)[0]*x1
                    # else:
                        # scaleFactor = left + np.diff(scaleArray)[0]*x2
            
        # # ----------------------Density minimization---------------------------
        
        # NoPeak = 0
        # while True: # Loop for the range shifting
            
            # density0 = density
            # densityArray = UtilityAnalysis.make_array_loop(density, densityStep, numSample)
            # chi2Array = np.zeros(numSample)
            
            # for i in range(len(densityArray)):
                
                # # ----------------Kaplow method for density--------------------

                # Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scaleFactor, Ibkg_Q)
                # alpha = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, \
                    # Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], \
                    # fe_Q[Q<=variables.QmaxIntegrate], Ztot, densityArray[i])
                # Icoh_Q = MainFunctions.calc_Icoh(alpha, Isample_Q, Iincoh_Q)

                # S_Q = MainFunctions.calc_SQ(Icoh_Q, Ztot, fe_Q, Sinf, Q, variables.minQ, \
                    # variables.QmaxIntegrate, variables.maxQ)
                # Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, 
                    # variables.smoothingFactor, \
                    # variables.minQ, variables.QmaxIntegrate, variables.maxQ)
                # SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Sinf, \
                    # dampingFunction)

                # # Utility.plot_data(Q, SsmoothDamp_Q, "I_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$I(Q)$", "y")

                # i_Q = MainFunctions.calc_iQ(SsmoothDamp_Q, Sinf)
                # F_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], \
                    # i_Q[Q<=variables.QmaxIntegrate])

                # Fopt_r, deltaFopt_r = Optimization.calc_optimize_Fr(variables.iterations, F_r, \
                    # Fintra_r, densityArray[i], i_Q[Q<=variables.QmaxIntegrate], Q[Q<=variables.QmaxIntegrate], \
                    # Sinf, J_Q[Q<=variables.QmaxIntegrate], r, variables.rmin, "n")

                # chi2Array[i] = simps(deltaFopt_r[r < variables.rmin]**2, r[r < variables.rmin])
            
            # # --------------------Range shifting selection---------------------
            
            # density = densityArray[np.argmin(chi2Array)] - scaleStep*1.1
            
            # nearIdx, nearEl = UtilityAnalysis.find_nearest(densityArray, density)
            
            # if nearIdx == 0:
                # density-=scaleStep*10
                # scaleStep*=10
                # NoPeak+=1
            # if nearIdx >= numSample-2:
                # density+=scaleStep*10
                # scaleStep*=10
                # NoPeak+=1

            # scaleStep /= 10

            # if (10*densityStep<=densityStepEnd)*(NoPeak>=5):
                # break
            
            # # -------------------chi2 curve fit for density--------------------
            
            # left = densityArray[0]
            # right = densityArray[-1]
            # coeffs = np.polyfit(densityArray, chi2Array, 3)
            # if (4*coeffs[1]**2 - 12*coeffs[2]*coeffs[0] < 0):
                    # density = (left+right)/2
            # else:
                # x1=(-2*coeffs[1]+(4*coeffs[1]**2-12*coeffs[2]*coeffs[0])**0.5)/(6*coeffs[0])
                # x2=(-2*coeffs[1]-(4*coeffs[1]**2-12*coeffs[2]*coeffs[0])**0.5)/(6*coeffs[0])
                # if (2*coeffs[1]+6*coeffs[0]*x1 > 2*coeffs[1]+6*coeffs[0]*x2):
                    # density = left + np.diff(densityArray)[0]*x1
                # else:
                    # density = left + np.diff(densityArray)[0]*x2

        # # -----------------------Step changing selection-----------------------
        
        # if np.abs(density-density0) > density0/25:
            # scaleStep = 0.006
            # densityStep = density0/10
        # elif np.abs(density-density0) > density0/75:
            # scaleStep = 0.0006
            # densityStep = density0/100
        # else:
            # scaleStep = 0.00006
            # densityStep = density0/1000

        # numLoopIteration += 1
        # if (np.abs(density-density0) < density0/2500 or numLoopIteration > 30):
            # break
    
    # print("scaleFactor final ", scaleFactor)
    # print("density final ", density)
    
    # plt.ioff()
    plt.show()