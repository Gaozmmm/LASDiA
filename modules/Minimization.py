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

from modules import Geometry
from modules import IgorFunctions
from modules import KaplowMethod
from modules import MainFunctions
from modules import Optimization
from modules import UtilityAnalysis
import time
from scipy.integrate import simps


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
    
    print("xMin", xMin)
    print("yMin", yMin)
    
    absXMin = xMin[np.argmin(yMin)]
    absYMin = np.amin(yMin)
    
    # absXMin = variableArray[0]+xMin[np.argmin(yMin)]*np.diff(variableArray)[0]
    # absYMin = polFit(absXMin)
    
    xFit = np.linspace(variableArray[0], variableArray[-1], 1000)
    yFit = polFit(xFit)
    
    return (xFit, yFit, absXMin, absYMin)


def OptimizeScale(Q, I_Q, Ibkg_Q, J_Q, Iincoh_Q, fe_Q, minQ, QmaxIntegrate, maxQ,
    Ztot, density, scaleFactor, Sinf, smoothingFactor, rmin, dampingFunction,
    Fintra_r, iterations, scaleStep, sth, s0th, thickness_sampling, phi_matrix,
    MCC_flag):
    """Function for the scale factor optimization.

    Q                  : numpy array
                         momentum transfer (nm^-1)
    I_Q                : numpy array
                         measured scattering intensity 
    Ibkg_Q             : numpy array
                         background scattering intensity
    J_Q                : numpy array
                         J(Q)
    Iincoh_Q           : numpy array
                         incoherent scattering intensity
    fe_Q               : numpy array
                         effective electric form factor
    minQ               : float
                         minimum Q value
    QmaxIntegrate      : float
                         maximum Q value for the intagrations
    maxQ               : float
                         maximum Q value
    Ztot               : int
                         total Z number
    density            : float
                         average atomic density
    scaleFactor        : float
                         scale factor
    Sinf               : float
                         value of S(Q) for Q->inf
    smoothingFactor    : float
                         smoothing factor
    rmin               : float
                         r cut-off value (nm)
    dampingFunction    : numpy array
                         damping function
    Fintra_r           : numpy array
                         intramolecular contribution of F(r)
    iterations         : int
                         number of iterations
    scaleStep
    sth
    s0th
    thickness_sampling : float
    phi_matrix
    MCC_flag
    """
    
    numSample = 23
    
    if MCC_flag.lower() == "y":
        T_MCC_sth, T_MCC_corr_factor_bkg = Geometry.MCCCorrection(sth, s0th,
            thickness_sampling, phi_matrix)
        
        I_Q = I_Q /T_MCC_sth
        Ibkg_Q  = Ibkg_Q * T_MCC_corr_factor_bkg / (T_MCC_sth)
    
    flag=0
    # Loop for the range shifting
    while 1:
        scaleArray = UtilityAnalysis.makeArrayLoop(scaleFactor, scaleStep)
        chi2Array = np.zeros(numSample)
        # flag+=1
        # print("iter flag ", flag)

        for i in range(numSample):
            # print("sample ", i)
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
            
            deltaFopt_r[np.where(r>=rmin)] = 0.0 # Igor version
            chi2Array[i] = np.mean(deltaFopt_r**2) # Igor version
            # Fth_r = Fintra_r - 4*np.pi*r*density
            # plt.plot(r, Fth_r)
            # plt.show()
            # print(Fth_r)
            # print(Fintra_r)
            # chi2Array[i] = simps(deltaFopt_r[np.where((r>0)&(r<rmin))]**2)

        # --------------------Range shifting selection --------------------
        
        plt.scatter(scaleArray, chi2Array)
        plt.grid(True)
        plt.show()

        maxLimit = 10**5
        if np.amax(chi2Array) >= maxLimit:
            scaleArray = scaleArray[0:np.argmax(chi2Array>=maxLimit)]
            chi2Array = chi2Array[0:np.argmax(chi2Array>=maxLimit)]
        
        if np.argmin(chi2Array) < 6:
            scaleFactor -= scaleStep*10
        if np.argmin(chi2Array) > 16:
            scaleFactor += scaleStep*10

        if (np.argmin(chi2Array)>=6 and np.argmin(chi2Array)<=16):
            break
        else:
            continue

    # ------------------------chi2 curve fit for scale-------------------------
    # plt.ioff()
    
    xFit, yFit, scaleFactor, chi2 = IgorFunctions.chi2Fit(scaleFactor, scaleArray, chi2Array)
    
    print("final scale factor", scaleFactor)

    plt.scatter(scaleArray, chi2Array)
    plt.plot(xFit, yFit)
    plt.grid(True)
    plt.show()
    
    return scaleFactor


def OptimizeDensity(Q, I_Q, Ibkg_Q, J_Q, Iincoh_Q, fe_Q, minQ, QmaxIntegrate, maxQ,
    Ztot, density, scaleFactor, Sinf, smoothingFactor, rmin, dampingFunction,
    Fintra_r, iterations, densityStep, sth, s0th,
    thickness_sampling, phi_matrix, MCC_flag):
    """Function for the density optimization.

    Q                  : numpy array
                         momentum transfer (nm^-1)
    I_Q                : numpy array
                         measured scattering intensity 
    Ibkg_Q             : numpy array
                         background scattering intensity
    J_Q                : numpy array
                         J(Q)
    Iincoh_Q           : numpy array
                         incoherent scattering intensity
    fe_Q               : numpy array
                         effective electric form factor
    minQ               : float
                         minimum Q value
    QmaxIntegrate      : float
                         maximum Q value for the intagrations
    maxQ               : float
                         maximum Q value
    Ztot               : int
                         total Z number
    density            : float
                         average atomic density
    scaleFactor        : float
                         scale factor
    Sinf               : float
                         value of S(Q) for Q->inf
    smoothingFactor    : float
                         smoothing factor
    rmin               : float
                         r cut-off value (nm)
    dampingFunction    : numpy array
                         damping function
    Fintra_r           : numpy array
                         intramolecular contribution of F(r)
    iterations         : int
                         number of iterations
    scaleStep
    sth
    s0th
    thickness_sampling : float
    phi_matrix
    MCC_flag
    """

    numSample = 23

    if MCC_flag.lower() == "y":
        T_MCC_sth, T_MCC_corr_factor_bkg = Geometry.MCCCorrection(sth, s0th,
            thickness_sampling, phi_matrix)
        
        I_Q = I_Q /T_MCC_sth
        Ibkg_Q  = Ibkg_Q * T_MCC_corr_factor_bkg / (T_MCC_sth)
    
    flag =0
    # Loop for the range shifting
    while 1:
        flag+=1
        # print("iter flag ", flag)

        densityArray = UtilityAnalysis.makeArrayLoop(density, densityStep)
        chi2Array = np.zeros(numSample)
        
        for i in range(numSample):
            # print("sample ", i, densityArray[i])

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

            # r, F_r = MainFunctions.calc_Fr(Q[Q<=QmaxIntegrate], 
            #     Qi_Q[Q<=QmaxIntegrate])

            
            Fopt_r, deltaFopt_r = Optimization.calc_optimize_Fr(iterations, F_r,
                Fintra_r, densityArray[i], i_Q[Q<=QmaxIntegrate], Q[Q<=QmaxIntegrate],
                Sinf, J_Q[Q<=QmaxIntegrate], r, rmin, "n")
            
            deltaFopt_r[np.where(r>=rmin)] = 0.0
            chi2Array[i] = np.mean(deltaFopt_r**2)
            
        # --------------------Range shifting selection --------------------
        
        # densityIdx = np.argmin(chi2Array)
        # density = densityArray[densityIdx] - densityStep*1.1
        
        # nearIdx, nearEl = UtilityAnalysis.find_nearest(densityArray, density)
        
        plt.scatter(densityArray, chi2Array)
        plt.grid(True)
        plt.show()

        maxLimit = 10**5
        if np.amax(chi2Array) >= maxLimit:
            densityArray = densityArray[0:np.argmax(chi2Array>=maxLimit)]
            chi2Array = chi2Array[0:np.argmax(chi2Array>=maxLimit)]

        # print(np.argmin(chi2Array))

        if np.argmin(chi2Array) < 6:
            # print("bug")
            density -= densityStep*10
        if np.argmin(chi2Array) > 16:
            density += densityStep*10

        

        if (np.argmin(chi2Array)>=6 and np.argmin(chi2Array)<=16):
            break
        else:
            continue
            
        
    # ------------------------chi2 curve fit for scale-------------------------
    
    #xFit, yFit, density, chi2Min = chi2Fit(densityArray, chi2Array)
    xFit, yFit, density, chi2 = IgorFunctions.chi2Fit(density, densityArray, chi2Array)
    
    print("final density", density)

    plt.scatter(densityArray, chi2Array)
    plt.plot(xFit, yFit)
    plt.grid(True)
    plt.show()

    return density


def OptimizeThickness(Q, I_Q, Ibkg_Q, J_Q, Iincoh_Q, fe_Q, maxQ, minQ, QmaxIntegrate,
    Ztot, density, scaleFactor, sth, s0th, Sinf, smoothingFactor, rmin, dampingFunction,
    Fintra_r, iterations, sthStep, ws1, ws2, r1, r2, d, phi_matrix_flag):
    """Function for the thickness optimization.
    """
    
    Flag = 0
    NoPeak = 0
    sth = max(sth-sthStep*11, 0.0)
    sthStepEnd = 0.0006
    numSample = 23
    
    # Loop for the range shifting
    while 1:
        sthArray = IgorFunctions.makeArrayLoop(sth, sthStep)
        chi2Array = np.zeros(numSample)
        
        for i in range(numSample):

            # ------------------Kaplow method for scale--------------------
            
            T_MCC_sth, T_MCC_corr_factor_bkg = Geometry.MCCCorrection(sthArray[i],
                s0th, thickness_sampling, phi_matrix)
    
            I_Q = I_Q /T_MCC_sth
            Ibkg_Q  = Ibkg_Q * T_MCC_corr_factor_bkg / (T_MCC_sth)
            
            Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scaleFactor, Ibkg_Q)
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
            sthIdx = np.argmin(chi2Array[0:np.argmax(chi2Array)])
        else:
            sthIdx = np.argmin(chi2Array)
        
        sth = sthArray[sthIdx]-sthStep*1.1
        nearIdx, nearEl = UtilityAnalysis.find_nearest(sthArray, sth)
        
        if nearIdx == 0:
            print("out1")
            sth -= sthStep*10
            sthStep *= 10
            NoPeak += 1
        if nearIdx >= numSample-2:
            print("out2")
            sth += sthStep*10
            sthStep *= 10
            NoPeak += 1
        
        sthStep /= 10
        Flag += 1
        
        if((10*sthStep>sthStepEnd) and (NoPeak<5)):
            continue
        else:
            break

    # ------------------------chi2 curve fit for scale-------------------------
    
    xFit, yFit, sth, chi2 = IgorFunctions.chi2Fit(sth, sthArray, chi2Array)
    
    print("final sample thickness", sth)
    
    return sth


def OptimizeThicknessRef(Q, I_Q, Ibkg_Q, J_Q, Iincoh_Q, fe_Q, maxQ, minQ, QmaxIntegrate,
    Ztot, density, scaleFactor, sth, s0th, Sinf, smoothingFactor, rmin, dampingFunction,
    Fintra_r, iterations, s0thStep, ws1, ws2, r1, r2, d, phi_matrix_flag):
    """Function for the thickness optimization.
    """
    
    Flag = 0
    NoPeak = 0
    s0th = max(s0th-s0thStep*11, 0.0)
    s0thStepEnd = 0.0006
    numSample = 23
    
    # Loop for the range shifting
    while 1:
        s0thArray = IgorFunctions.makeArrayLoop(s0th, s0thStep)
        chi2Array = np.zeros(numSample)
        
        for i in range(numSample):

            # ------------------Kaplow method for scale--------------------
            
            T_MCC_sth, T_MCC_corr_factor_bkg = Geometry.MCCCorrection(sth,
                s0thArray[i], thickness_sampling, phi_matrix)
    
            I_Q = I_Q /T_MCC_sth
            Ibkg_Q  = Ibkg_Q * T_MCC_corr_factor_bkg / (T_MCC_sth)
            
            Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scaleFactor, Ibkg_Q)
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
            s0thIdx = np.argmin(chi2Array[0:np.argmax(chi2Array)])
        else:
            s0thIdx = np.argmin(chi2Array)
        
        s0th = s0thArray[s0thIdx]-s0thStep*1.1
        nearIdx, nearEl = UtilityAnalysis.find_nearest(s0thArray, s0th)
        
        if nearIdx == 0:
            print("out1")
            s0th -= s0thStep*10
            s0thStep *= 10
            NoPeak += 1
        if nearIdx >= numSample-2:
            print("out2")
            s0th += s0thStep*10
            s0thStep *= 10
            NoPeak += 1
        
        s0thStep /= 10
        Flag += 1
        
        if((10*s0thStep>s0thStepEnd) and (NoPeak<5)):
            continue
        else:
            break

    # ------------------------chi2 curve fit for scale-------------------------
    
    xFit, yFit, sth, chi2 = IgorFunctions.chi2Fit(s0th, s0thArray, chi2Array)
    
    print("final sample thickness ref", s0th)
    
    return s0th


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