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

"""Set of modules used in LASDiA to calculate the chi2 minimization.

The nomenclature and the procedure follow the article:
Eggert et al. 2002 PRB, 65, 174105.

For the functions arguments and the returns I followed this convetion for the notes:
arguments: description - type
returns: description - type.

For the variables name I used this convention:
if the variable symbolizes a function, its argument is preceded by an underscore: f(x) -> f_x
otherwise it is just the name.
"""


import numpy as np
import matplotlib.pyplot as plt

from modules import IgorFunctions
from modules import KaplowMethod
from modules import MainFunctions
from modules import Optimization
from modules import UtilityAnalysis
import time


def chi2Fit(variableArray, chi2Array):
    """Function to fit the chi2 function to find the minimal value.
    
    Parameters
    ----------
    variableArray : numpy array
                    array with variable values
    chi2Array     : numpy array
                    array with chi2 values
    
    Returns
    -------
    xFit          : numpy array
                    array with the x values to plot the fit results
    yFit          : numpy array
                    array with the fit y values
    absXMin       : float
                    fit curve minimum value abscissa
    absYMin       : float
                    fit curve minimum value ordinate
    """
    
    polFit = np.poly1d(np.polyfit(variableArray, chi2Array, 3))
    
    polRoots = polFit.deriv().r
    realPolRoots = polRoots[polRoots.imag==0].real
    test = polFit.deriv(2)(realPolRoots)
    
    xMin = realPolRoots[test>0]
    yMin = polFit(xMin)
    
    absXMin = xMin[np.argmin(yMin)]
    absYMin = np.amin(yMin)
    
    # absXMin = variableArray[0]+xMin[np.argmin(yMin)]*np.diff(variableArray)[0]
    # absYMin = polFit(absXMin)
    
    xFit = np.linspace(variableArray[0], variableArray[-1], 1000)
    yFit = polFit(xFit)
    
    return (xFit, yFit, absXMin, absYMin)


def OptimizeScale(Q, I_Q, Ibkg_Q, J_Q, Iincoh_Q, fe_Q, maxQ, minQ, QmaxIntegrate,
    Ztot, density, scaleFactor, Sinf, smoothingFactor, rmin, dampingFunction,
    Fintra_r, iterations, scaleStep):
    """Function for the scale factor optimization.
    """
    
    Flag = 0
    NoPeak = 0
    scaleStepEnd = 0.00006
    numSample = 23
    
    # plt.ion()
    # figure, ax = plt.subplots()

    # Loop for the range shifting
    # print("cond", (10*scaleStep>scaleStepEnd) * (NoPeak<5) * ((Flag==1) + (scaleFactor+scaleStep*1.1>=0)))
    # while ((10*scaleStep>scaleStepEnd) * (NoPeak<5) * ((Flag==1) + (scaleFactor+scaleStep*1.1>=0))):
    while 1:
        scaleArray = UtilityAnalysis.makeArrayLoop(scaleFactor, scaleStep)
        chi2Array = np.zeros(numSample)
        
        # ax.cla()
        # ax.grid(True)
        # ax.relim()
        # ax.autoscale_view()
        
        for i in range(numSample):

            # ------------------Kaplow method for scale--------------------

            Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scaleArray[i], Ibkg_Q)
            alpha = MainFunctions.calc_alpha(J_Q[Q<=QmaxIntegrate], Sinf, 
                Q[Q<=QmaxIntegrate], Isample_Q[Q<=QmaxIntegrate], 
                fe_Q[Q<=QmaxIntegrate], Ztot, density)
            Icoh_Q = MainFunctions.calc_Icoh(alpha, Isample_Q, Iincoh_Q)

            S_Q = MainFunctions.calc_SQ(Icoh_Q, Ztot, fe_Q, Sinf, Q, minQ, 
                QmaxIntegrate, maxQ)
            
            Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, 
                smoothingFactor, minQ, QmaxIntegrate, maxQ)
            
            SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Sinf,
                dampingFunction)
            
            i_Q = MainFunctions.calc_iQ(SsmoothDamp_Q, Sinf)
            
            Qi_Q = Q*i_Q
            r, F_r = MainFunctions.calc_Fr(Q[Q<=QmaxIntegrate], 
                Qi_Q[Q<=QmaxIntegrate])

            Fopt_r, deltaFopt_r = Optimization.calc_optimize_Fr(iterations, F_r,
                Fintra_r, density, i_Q[Q<=QmaxIntegrate], Q[Q<=QmaxIntegrate],
                Sinf, J_Q[Q<=QmaxIntegrate], r, rmin, "n")
            
            deltaFopt_r[np.where(r>=rmin)] = 0.0
            chi2Array[i] = np.mean(deltaFopt_r**2)

        # --------------------Range shifting selection --------------------
        
        if np.amax(chi2Array) > 10**8:
            scaleFactorIdx = np.argmin(chi2Array[0:np.argmax(chi2Array)])
            scaleFactor = scaleArray[scaleFactorIdx]
        else:
            scaleFactorIdx = np.argmin(chi2Array)
            scaleFactor = scaleArray[scaleFactorIdx]
        
        # nearIdx, nearEl = UtilityAnalysis.find_nearest(scaleArray, scaleFactor)
        
        print(scaleFactor)
        
        if scaleFactorIdx == 0:
            print("out1")
            scaleFactor -= scaleStep*10
            scaleStep *= 10
            NoPeak += 1
        if scaleFactorIdx >= numSample-2:
            print("out2")
            scaleFactor += scaleStep*10
            scaleStep *= 10
            NoPeak += 1
        
        plt.scatter(scaleArray, chi2Array)
        plt.grid(True)
        plt.show()
        
        scaleStep /= 10
        Flag += 1
        print(Flag, scaleFactor, scaleStep)
        if (scaleFactorIdx>=6 and scaleFactorIdx<=16):
            break

    # ------------------------chi2 curve fit for scale-------------------------
    # plt.ioff()

    xFit, yFit, scaleFactor, chi2Min = chi2Fit(scaleArray, chi2Array)
    # plt.plot(xFit, yFit)
    # plt.plot(scaleFactor, chi2Min, "r")
    # plt.grid(True)
    
    print("final scale factor", scaleFactor)
    # plt.show()
    
    return scaleFactor


def OptimizeDensity(Q, I_Q, Ibkg_Q, J_Q, Iincoh_Q, fe_Q, maxQ, minQ, QmaxIntegrate,
    Ztot, density, scaleFactor, Sinf, smoothingFactor, rmin, dampingFunction,
    Fintra_r, iterations, densityStep, densityStepEnd):
    """Function for the density optimization.
    """

    Flag = 0
    NoPeak = 0
    numSample = 23
    
    # Loop for the range shifting
    # while ((10*densityStep>densityStepEnd) and (NoPeak<5)):
    while 1:
        densityArray = UtilityAnalysis.makeArrayLoop(density, densityStep)
        chi2Array = np.zeros(numSample)

        for i in range(numSample):

            # ------------------Kaplow method for scale--------------------

            Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scaleFactor, Ibkg_Q)
            alpha = MainFunctions.calc_alpha(J_Q[Q<=QmaxIntegrate], Sinf, 
                Q[Q<=QmaxIntegrate], Isample_Q[Q<=QmaxIntegrate], 
                fe_Q[Q<=QmaxIntegrate], Ztot, densityArray[i])
            Icoh_Q = MainFunctions.calc_Icoh(alpha, Isample_Q, Iincoh_Q)

            S_Q = MainFunctions.calc_SQ(Icoh_Q, Ztot, fe_Q, Sinf, Q, minQ, 
                QmaxIntegrate, maxQ)

            Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, 
                smoothingFactor, minQ, QmaxIntegrate, maxQ)

            SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Sinf,
                dampingFunction)

            i_Q = MainFunctions.calc_iQ(SsmoothDamp_Q, Sinf)
            
            Qi_Q = Q*i_Q
            r, F_r = MainFunctions.calc_Fr(Q[Q<=QmaxIntegrate], 
                Qi_Q[Q<=QmaxIntegrate])

            Fopt_r, deltaFopt_r = Optimization.calc_optimize_Fr(iterations, F_r,
                Fintra_r, densityArray[i], i_Q[Q<=QmaxIntegrate], Q[Q<=QmaxIntegrate],
                Sinf, J_Q[Q<=QmaxIntegrate], r, rmin, "n")
            
            deltaFopt_r[np.where(r>=rmin)] = 0.0
            chi2Array[i] = np.mean(deltaFopt_r**2)
            
        # --------------------Range shifting selection --------------------
        
        densityIdx = np.argmin(chi2Array)
        density = densityArray[densityIdx]
        
        # nearIdx, nearEl = UtilityAnalysis.find_nearest(densityArray, density)
        
        if densityIdx == 0:
            print("out3")
            density -= densityStep*10
            densityStep *= 10
            NoPeak += 1
        if densityIdx >= numSample-2:
            print("out4")
            density += densityStep*10
            densityStep *= 10
            NoPeak += 1

        densityStep /= 10
        
        Flag += 1
        print(Flag, density, densityStep)
        plt.scatter(densityArray, chi2Array)
        plt.grid(True)
        plt.show()
        if (densityIdx>=6 and densityIdx<=16):
            break
        
        
    # ------------------------chi2 curve fit for scale-------------------------
    
    xFit, yFit, density, chi2Min = chi2Fit(densityArray, chi2Array)

    print("final density", density)
    return density


def chi2_minimization(scaleFactor, Q, I_Q, Ibkg_Q, J_Q, fe_Q, Iincoh_Q, Sinf, Ztot,
    density, Fintra_r, r, minQ, QmaxIntegrate, maxQ, smoothFactor, dampFactor, iteration, rmin):
    """Function to calculate the whole loop for the chi2 minimization.
    """
    
    scaleStep = 0.05
    densityStep = 0.025
    numSample = 23
    numLoopIteration = 0
    
    plt.ion()
    figure, ax = plt.subplots()

    while True:
        NoPeak = 0
        while True:
            ax.cla()
            ax.grid(True)
            scaleArray = UtilityAnalysis.make_array_loop(scaleFactor, scaleStep, numSample)

            chi2Array = np.zeros(numSample)

            plt.xlabel("Scale")
            ax.relim()
            ax.autoscale_view()
            for i in range(len(scaleArray)):
                chi2Array[i], SsmoothDamp_Q, F_r, Fopt_r = KaplowMethod.Kaplow_method(Q, I_Q,
                    Ibkg_Q, J_Q, fe_Q, Iincoh_Q, Sinf, Ztot, scaleArray[i], density, Fintra_r, r,
                    minQ, QmaxIntegrate, maxQ, smoothFactor, dampFactor, iteration, rmin)
                
                plt.scatter(scaleArray[i], chi2Array[i])
                figure.canvas.draw()
            
            if np.amax(chi2Array) > 10**8:
                chi2Array = chi2Array[0:np.amax(chi2Array)]
            
            scaleFactor = scaleArray[np.argmin(chi2Array)] - scaleStep*1.1
            
        ax.cla()
        ax.grid(True)

        density0 = density
        densityArray = UtilityAnalysis.make_array_loop(density, densityStep, numSample)
        chi2Array = np.zeros(numSample)
        
        plt.xlabel("Density")
        ax.relim()
        ax.autoscale_view()
        for i in range(len(densityArray)):
            chi2Array[i], SsmoothDamp_Q, F_r, Fopt_r = KaplowMethod.Kaplow_method(Q, I_Q,
                Ibkg_Q, J_Q, fe_Q, Iincoh_Q, Sinf, Ztot, scaleFactor, densityArray[i], Fintra_r, r,
                minQ, QmaxIntegrate, maxQ, smoothFactor, dampFactor, iteration, rmin)
            
            plt.scatter(densityArray[i], chi2Array[i])
            figure.canvas.draw()
            
        
        xfit, yfit, density = chi2_fit(densityArray, chi2Array)
        plt.plot(xfit, yfit)
        figure.canvas.draw()
        
        if np.abs(density-density0) > density0/25:
            scaleStep = 0.006
            densityStep = density0/10
        elif np.abs(density-density0) > density0/75:
            scaleStep = 0.0006
            densityStep = density0/100
        else:
            scaleStep = 0.00006
            densityStep = density0/1000

        numLoopIteration += 1
        if (np.abs(density-density0) < density0/2500 or numLoopIteration > 30):
            break
    
    plt.ioff()
    
    return (density, scaleFactor)