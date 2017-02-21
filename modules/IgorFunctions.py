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

"""Module containing the functions used in Igor Pro code to compare their results 
with those of LASDiA.

The nomenclature and the procedure follow the article:
Eggert et al. 2002 PRB, 65, 174105.

For the functions arguments and the returns I followed this convetion for the
notes:
arguments: description - type
returns: description - type.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by
an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""


import numpy as np
#from scipy.interpolate import UnivariateSpline
from scipy.integrate import simps
from scipy import fftpack
import math
import matplotlib.pyplot as plt

# from modules import Formalism
# from modules import Geometry
# from modules import IgorFunctions
# from modules import KaplowMethod
from modules import MainFunctions
from modules import Minimization
# from modules import Optimization
from modules import Utility
from modules import UtilityAnalysis


def calc_JQ(Iincoh_Q, fe_Q):
    """Function to calculate J(Q) (eq. 35)

    arguments:
    Iincoh_Q: incoherent scattering intensity - array
    fe_Q: effective electric form factor - array

    returns:
    J_Q: J(Q) - array
    """

    J_Q = Iincoh_Q/fe_Q**2

    return J_Q


def calc_Subt(I_Q, Ibkg_Q, scaleFactor):
    """
    """
    
    # Subt = (I_Q - scaleFactor*Ibkg_Q)/absCorrFactor
    Subt = (I_Q - scaleFactor*Ibkg_Q)
    
    return Subt


def calc_alpha(J_Q, Sinf, Q, Isample_Q, fe_Q, Ztot, rho0):
    """Function to calculate the normalization factor alpha (eq. 34)

    arguments:
    J_Q: J(Q) - array
    Sinf: Sinf - number
    Q: momentum transfer - array
    Isample_Q: sample scattering intensity - array
    fe_Q: effective electric form factor - array
    Ztot: total Z number - number
    rho0: average atomic density - number

    returns:
    alpha: normalization factor - number
    """

    Integral1 = simps(Q**2*(J_Q + Sinf*Ztot**2), Q)
    Integral2 = simps(Isample_Q*Q**2/fe_Q**2, Q)
    alpha = ((-2*np.pi**2*Ztot**2*rho0) + Integral1) / Integral2

    return alpha


def calc_SQ(Q, Subt, alpha, fe_Q, J_Q, Ztot, Sinf, QmaxIntegrate):
    """Function to calculate the structure factor S(Q) (eq. 18) with Igor formula

    arguments:
    Ztot: total Z number - number
    fe_Q: effective electric form factor - array
    Sinf: Sinf - number
    Q: momentum transfer - array
    min_index: array index of element with Q<minQ - array
    max_index: array index of element with Q>QmaxIntegrate & Q<=maxQ - array
    calculation_index: array index of element in the calculation range Q>minQ & Q<=QmaxIntegrate - array

    returns:
    S_Q: structure factor - array
    """
    
    S_Q = (alpha*Subt/fe_Q**2-J_Q)/Ztot**2
    S_Q[Q>QmaxIntegrate] = Sinf
    
    return S_Q


def absorption(Q):
    """Function to calculate the absorption correction with Igor formula.
    """
    
    xD=0.0
    Lambda=0.03778
    eD=0.17
    BRad=1.2
    eDAC=1.25
    Gasket=0.00525
    InnerAperature=0.31
    OuterAperature=1.639
    DACAperature=1.639
    muB=0.0
    
    SeatAbs_c = np.array([0,0.997566,0.0,0.0,0.00107,0.000126,0.00044,
        4.6E-05,2E-05,7.1E-05,4.9E-05,6.6E-05,0.000335,3.2E-06,
        2.9E-05,3.4E-05,1.5E-05,0.0,1.9E-06,1.7E-05,9.4E-05,5E-06,1.2E-05,0])

    SeatAbs_x0 = np.array([2.3518,1.0,1.0,10.32,10.32,10.32,10.32,
        10.32,10.32,10.32,10.32,10.32,10.32,10.32,10.32,10.32,10.32,
        10.5744,16.1851,18.0876,13.477,15.1361,15.9401,1.0])

    SeatAbs_y0 = np.array([0.16609,0.016509,0.0360423,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,9.7903,
        -4.87898,0.0,3.76538,-5.24716,0.0360423])
    
    SeatAbs_A1 = np.array([0.40924,0.0,0.0,3.14791,11.8267,14.8745,
        14.8745,50.7179,53.929,63.2412,82.8997,92.8759,102.418,103.483,
        115.173,126.886,132.891,118.909,74.5406,53.752,105.304,96.6244,
        89.9649,0.0])

    SeatAbs_A2 = np.array([5.8379,0.0,0.0,1.68147,6.01518,7.54942,
        7.54942,32.0789,34.5927,36.0429,42.4638,43.5787,51.3385,
        63.0198,75.1925,68.467,79.5591,91.9556,24.7572,44.489,
        69.6294,58.604,62.8023,0.0])

    SeatAbs_t1 = np.array([16.73,1.0,1.0,2.55004,2.67334,2.69543,
        2.69543,2.73906,2.72493,2.80362,2.9916,3.11676,3.0976,2.94736,
        2.86549,3.16199,3.08399,2.83517,6.64092,4.76416,3.87274,10.048,
        4.58574,1.0])

    SeatAbs_t2 = np.array([4.6014,1.0,1.0,8.26453,8.92972,9.0562,9.0562,
        9.00571,9.07263,9.55472,9.97426,10.4699,10.3523,9.72451,9.59644,
        10.3455,10.0967,9.55927,4.76267,20.1567,13.2047,3.39669,18.6351,
        1.0])

    # muB=0.0230077
    if muB==0.0:
        for i in range(len(SeatAbs_y0)):
            muB += SeatAbs_c[i] * muFunc(SeatAbs_y0[i],SeatAbs_x0[i],SeatAbs_A1[i],
                SeatAbs_t1[i],SeatAbs_A2[i],SeatAbs_t2[i])
    
    
    yA=InnerAperature/2
    xA=eD-xD
    Theta1=np.arctan(yA/xA)
    Q1=4*np.pi/Lambda*np.sin(Theta1/2)
    yB=OuterAperature/2
    xB=(BRad**2-yB**2)**0.5-xD
    Theta2=np.arctan(yB/xB)
    Q2=4*np.pi/Lambda*np.sin(Theta2/2)
    yG=DACAperature/2
    xG=(BRad**2-yG**2)**0.5-xD
    ThetaDAC=np.arctan(yG/xG)
    QDAC=4*np.pi/Lambda*np.sin(ThetaDAC/2)
    yH=Gasket/2
    xH=-xD
    if xH==0.0:
        ThetaGasket=np.pi/2
    elif xH<0.0:
        ThetaGasket=np.arctan(yH/xH)+np.pi
    else:
        ThetaGasket=np.arctan(yH/xH)
    QGasket=4*np.pi/Lambda*np.sin(ThetaGasket/2)
    m=(yB-yA)/(xB-xA)
    rhoD=3.5155
    muD=muFunc(SeatAbs_y0[0],SeatAbs_x0[0],SeatAbs_A1[0],SeatAbs_t1[0],SeatAbs_A2[0],SeatAbs_t2[0])
    rhoB=2.34
    alphaDAC=0.8
    
    Theta=2*np.arcsin(Lambda*Q/(4*np.pi))
    dOI=xA/np.cos(Theta)
    dOJ=((yA-m*xA)/(np.tan(Theta)-m))/np.cos(Theta)
    dOF=(BRad**2-xD**2*np.sin(Theta))**0.5-xD*np.cos(Theta)
    dOG=(eDAC-xD)/np.cos(Theta)
    
    AbsRefCalc = np.zeros(len(Q))
    
    if (Theta1<=Theta2):
        AbsRefCalc[Q<=Q1] = np.exp(-muD*rhoD*dOI[Q<=Q1])
        AbsRefCalc[(Q>Q1) & (Q<=Q2)] += np.exp(-(muD*rhoD*dOI[(Q>Q1) & (Q<=Q2)]
            +muB*rhoB*(dOJ[(Q>Q1) & (Q<=Q2)]-dOI[(Q>Q1) & (Q<=Q2)])))
        AbsRefCalc[(Q>Q2) & (Q<=QDAC)] += np.exp(-(muD*rhoD*dOI[(Q>Q2) & (Q<=QDAC)]
            +muB*rhoB*(dOF[(Q>Q2) & (Q<=QDAC)]-dOI[(Q>Q2) & (Q<=QDAC)])))
        AbsRefCalc[Q>QDAC] += np.exp(-(muD*rhoD*dOI[Q>QDAC]+muB*rhoB*(dOF[Q>QDAC]
            -dOI[Q>QDAC])+alphaDAC*(dOG[Q>QDAC]-dOF[Q>QDAC])))
        AbsRefCalc[Q<QGasket]*=1
    else:
        AbsRefCalc[Q<=Q2] = np.exp(-muD*rhoD*dOI[Q<=Q2])
        AbsRefCalc[(Q>Q2) & (Q<=Q1)] += np.exp(-(muD*rhoD*dOI[(Q>Q2) & (Q<=Q1)]
            +muB*rhoB*(dOF[(Q>Q2) & (Q<=Q1)]-dOJ[(Q>Q2) & (Q<=Q1)])))
        AbsRefCalc[(Q>Q1) & (Q<=QDAC)] += np.exp(-(muD*rhoD*dOI[(Q>Q1) & (Q<=QDAC)]
            +muB*rhoB*(dOF[(Q>Q1) & (Q<=QDAC)]-dOI[(Q>Q1) & (Q<=QDAC)])))
        AbsRefCalc[Q>QDAC] += np.exp(-(muD*rhoD*dOI[Q>QDAC]
            +muB*rhoB*(dOF[Q>QDAC]-dOI[Q>QDAC])+alphaDAC*(dOG[Q>QDAC]-dOF[Q>QDAC])))
        AbsRefCalc[Q<QGasket]*=1
    
    return AbsRefCalc


def muFunc(y0,x0,A1,t1,A2,t2):
    """absorption coefficient function"""
    
    x = 33.1952 #x-ray energy (kev)
    mu = y0 + A1*np.exp(-(x-x0)/t1) + A2*np.exp(-(x-x0)/t2) 
    return mu


def calc_FFT_QiQ(Q, Qi_Q, QmaxIntegrate):
    """Function to calculate the FFT following the Igor Pro procedure.
    I do not agree with this procedure, but I follow it to compare my results
    with Igor Pro's ones.
    
    Parameters
    ----------
    Q             : numpy array
                    momentum transfer (nm^-1)
    i_Q           : numpy array
                    i(Q)
    
    QmaxIntegrate : float
                    maximum Q value for the integrations
    
    Returns
    -------
    r             : numpy array
                    atomic distance (nm)
    F_r           : numpy array
                    F(r)
    """

    pMax, elem = UtilityAnalysis.find_nearest(Q, QmaxIntegrate)
    NumPoints = 2*2*2**math.ceil(math.log(5*(pMax+1))/math.log(2))
    DelR = 2*np.pi/(np.mean(np.diff(Q))*NumPoints)
    # Qi_Q = Utility.resize_zero(Q[Q<=QmaxIntegrate]*i_Q[Q<=QmaxIntegrate], NumPoints)
    Qi_Q = Utility.resize_zero(Qi_Q[Q<=QmaxIntegrate], NumPoints)
    Qi_Q[pMax+1:] = 0.0
    Q = np.arange(np.amin(Q), np.amin(Q)+np.mean(np.diff(Q))*NumPoints, np.mean(np.diff(Q)))
    r = MainFunctions.calc_r(Q)
    F_r = fftpack.fft(Qi_Q)
    F_r = F_r[np.where(r>=0.0)]
    F_r = -np.imag(F_r)*np.mean(np.diff(Q))*2/np.pi
    r = np.arange(0.0, 0.0+DelR*len(F_r), DelR)
    
    return (r, F_r)


def calc_IFFT_Fr(r, F_r):
    """Function to calculate the FFT following the IGOR procedure.
    I do not agree with this procedure, but I follow it to compare my results
    with Igor Pro's ones.
    
    Parameters
    ----------
    r      : numpy array
             atomic distance (nm)
    F_r    : numpy array
             F(r)
    
    Returns
    -------
    Q      : numpy array
             momentum transfer (nm^-1)
    Qi_Q   : numpy array
             Qi(Q)
    """
    # Utility.write_file("./F_r.txt", r, F_r)
    NumPoints = 2**math.ceil(math.log(len(F_r)-1)/math.log(2))
    F_r = Utility.resize_zero(F_r, NumPoints)
    Q = np.linspace(0.0, 109, 550)
    DelQ = 2*np.pi/(np.mean(np.diff(r))*NumPoints)
    meanDeltar = np.mean(np.diff(r))
    Q1 = fftpack.fftfreq(r.size, meanDeltar)
    Qi_Q = fftpack.fft(F_r)
    Qi_Q = Qi_Q[np.where(Q1>=0.0)]
    Qi_Q = -np.imag(Qi_Q)*meanDeltar
    # print(meanDeltar)
    Q1 = np.arange(0.0, 0.0+DelQ*len(Qi_Q), DelQ)
    # print(len(Q1), len(Qi_Q))
    # Utility.write_file("./Qi_Q.txt", Q1, Qi_Q)
    
    idxArray = np.zeros(550, dtype=np.int)
    for i in range(len(Q)):
        idxArray[i], _ = UtilityAnalysis.find_nearest(Q1, Q[i])
    Qi_Q = Qi_Q[idxArray]
    
    return (Q, Qi_Q)


def FitRemoveGofRPeaks(Q, SsmoothDamp_Q, Sinf, QmaxIntegrate, Fintra_r,
    iterations, rmin, density, J_Q, Ztot):
    """Equivalent of Igor function
    """
    
    Qi_Q = Q*MainFunctions.calc_iQ(SsmoothDamp_Q, Sinf)
    r, GR = calc_FFT_QiQ(Q, Qi_Q, QmaxIntegrate)
    
    # _, Fintra_r = UtilityAnalysis.rebinning(rintra, Fintra_r, np.amin(Fintra_r), 
        # np.amax(Fintra_r), len(GR))

    QiQ1 = np.zeros(len(SsmoothDamp_Q))
    idx, _ = UtilityAnalysis.find_nearest(Qi_Q, QmaxIntegrate)
    QiQ1[Q<QmaxIntegrate] = Qi_Q[Q<QmaxIntegrate]
    QiQ1[0] = 0.0
    
    GR1 = GR
    DelG = np.zeros(len(GR))
    Rnn = rmin
    DelG[r<Rnn] = GR1[r<Rnn]-(Fintra_r[r<Rnn]-4*np.pi*r[r<Rnn]*density)
    
    GRIdealSmallR = Fintra_r-4*np.pi*r*density
    
    for i in range(iterations):
        Q1, QiQCorr = calc_IFFT_Fr(r, DelG)
        mask = np.where((Q1>0.0) & (Q1<QmaxIntegrate))
        QiQ1[mask] = QiQ1[mask] - (QiQ1[mask] / 
            (Q1[mask] *(Sinf + J_Q[:len(Q1[mask])]/Ztot**2)) + 1) * QiQCorr[mask]
        
        r, GR1 = calc_FFT_QiQ(Q1, QiQ1, QmaxIntegrate)
        
        DelG = np.zeros(len(GR1))
        DelG[r<Rnn] = GR1[r<Rnn]-GRIdealSmallR[r<Rnn]
        
        _, rmin = UtilityAnalysis.find_nearest(r, 0.95*Rnn)
        Rnn = 0.99*r[np.where(GR1==np.amin(GR1[r>=rmin]))[0][0]]

    DTemp = DelG
    DTemp = DTemp**2
    chi2 = np.mean(DTemp)

    return chi2


def OptimizeScaleGofRCorr(Q, I_Q, Ibkg_Q, absCorrFactor, J_Q, fe_Q, maxQ, minQ,
    QmaxIntegrate, Ztot, density, scaleFactor, Sinf, smoothingFactor, rmin, NumPoints,
    dampingFunction, Fintra_r, iterations, scaleStep):
    """Function for the scale factor optimization.
    """
    
    Flag = 0
    NoPeak = 0
    # scaleStep = 0.05
    scaleStepEnd = 0.00006
    numSample = 23
    scaleFactor = scaleFactor-scaleStep*11
    Qorg = Q
    
    # Loop for the range shifting
    while ((10*scaleStep>scaleStepEnd) and (NoPeak<5) and ((Flag==1) or (scaleFactor+scaleStep*1.1>=0))):
        # print("scale loop")
        scaleArray = make_array_loop(scaleFactor, scaleStep, numSample)
        chi2Array = np.zeros(numSample)
        
        for i in range(len(scaleArray)):
            
            # ------------------Kaplow method for scale--------------------
            Q = Qorg
            Subt = calc_Subt(I_Q, Ibkg_Q, scaleArray[i])
            alpha = calc_alpha(J_Q[Q<=QmaxIntegrate], Sinf, 
                Q[Q<=QmaxIntegrate], Subt[Q<=QmaxIntegrate],
                fe_Q[Q<=QmaxIntegrate], Ztot, density)
            
            S_Q = calc_SQ(Q, Subt, alpha, fe_Q, J_Q, Ztot, Sinf, 
                QmaxIntegrate)
                
            Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, 
                smoothingFactor, 
                minQ, QmaxIntegrate, maxQ)
            Q, Ssmooth_Q = UtilityAnalysis.rebinning(Q, Ssmooth_Q, 0.0, 
                maxQ, NumPoints)
            SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Sinf,
                dampingFunction)
            
            chi2Array[i] = FitRemoveGofRPeaks(Q, SsmoothDamp_Q, Sinf, 
                QmaxIntegrate, Fintra_r, iterations, 
                rmin, density, J_Q, Ztot)
        
        # --------------------Range shifting selection --------------------
        
        if np.amax(chi2Array) > 10**8:
            scaleFactor = scaleArray[np.argmin(chi2Array[0:np.argmax(chi2Array)])] - scaleStep*1.1
        else:
            scaleFactor = scaleArray[np.argmin(chi2Array)] - scaleStep*1.1
        
        nearIdx, nearEl = UtilityAnalysis.find_nearest(scaleArray, scaleFactor)

        if nearIdx == 0:
            print("out1")
            scaleFactor -= scaleStep*10
            scaleStep *= 10
            NoPeak += 1
        if nearIdx >= numSample-2:
            print("out2")
            scaleFactor += scaleStep*10
            scaleStep *= 10
            NoPeak += 1

        scaleStep /= 10
        Flag += 1
        print(Flag)

    # ------------------------chi2 curve fit for scale-------------------------

    xFit, yFit, scaleFactor, chi2Min = chi2Fit(scaleFactor, scaleArray, chi2Array)
    # plt.plot(xFit, yFit)
    # plt.grid(True)
    
    print("final scale factor", scaleFactor)
    # plt.show()
    
    return scaleFactor


def OptimizeDensityGofRCorr(Q, I_Q, Ibkg_Q, absCorrFactor, J_Q, fe_Q, maxQ, minQ,
    QmaxIntegrate, Ztot, density, scaleFactor, Sinf, smoothingFactor, rmin, NumPoints,
    dampingFunction, Fintra_r, iterations, densityStep, densityStepEnd):
    """Function for the density optimization.
    """
    
    Flag = 0
    NoPeak = 0
    # densityStep = density/50
    # densityStepEnd = density/250
    numSample = 23
    density = density-densityStep*11
    Qorg = Q
    
    while ((10*densityStep>densityStepEnd) and (NoPeak<5)): # Loop for the range shifting
        densityArray = make_array_loop(density, densityStep, numSample)
        chi2Array = np.zeros(numSample)
        
        for i in range(len(densityArray)):
            
            # ------------------Kaplow method for scale--------------------
            Q = Qorg
            Subt = calc_Subt(I_Q, Ibkg_Q, scaleFactor)
            alpha = calc_alpha(J_Q[Q<=QmaxIntegrate], Sinf, 
                Q[Q<=QmaxIntegrate], Subt[Q<=QmaxIntegrate],
                fe_Q[Q<=QmaxIntegrate], Ztot, densityArray[i])
            
            S_Q = calc_SQ(Q, Subt, alpha, fe_Q, J_Q, Ztot, Sinf, 
                QmaxIntegrate)
                
            Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, 
                smoothingFactor, 
                minQ, QmaxIntegrate, maxQ)
            Q, Ssmooth_Q = UtilityAnalysis.rebinning(Q, Ssmooth_Q, 0.0, 
                maxQ, NumPoints)
            SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Sinf,
                dampingFunction)
            
            chi2Array[i] = FitRemoveGofRPeaks(Q, SsmoothDamp_Q, Sinf, 
                QmaxIntegrate, Fintra_r, iterations, 
                rmin, densityArray[i], J_Q, Ztot)
        
        # --------------------Range shifting selection --------------------

        density = densityArray[np.argmin(chi2Array)] - densityStep*1.1
        
        nearIdx, nearEl = UtilityAnalysis.find_nearest(densityArray, density)
        
        if nearIdx == 0:
            print("out3")
            density -= densityStep*10
            densityStep *= 10
            NoPeak += 1
        if nearIdx >= numSample-2:
            print("out4")
            density += densityStep*10
            densityStep *= 10
            NoPeak += 1

        densityStep /= 10
        
        Flag += 1
        print(Flag)
        plt.scatter(densityArray, chi2Array)
        plt.grid(True)
        # plt.show
        # if ((10*densityStep<=densityStepEnd) and (NoPeak>=5)):
            # break
        
    # ------------------------chi2 curve fit for scale-------------------------
    
    xFit, yFit, density, chi2Min = chi2Fit(density, densityArray, chi2Array)

    print("final density", density)
    return density


def chi2Fit(value, valueArray, chi2Array):
    """Function to fit the chi2 plot.
    """
    
    if value < 0:
        print("Scale factor < 0")
        # break
    else:
        left = valueArray[0]
        right = valueArray[-1]
        coeffs = np.polyfit(valueArray, chi2Array, 3)
        if (4*coeffs[1]**2 - 12*coeffs[2]*coeffs[0] < 0):
            value = (left+right)/2
        else:
            x1=(-2*coeffs[1]+(4*coeffs[1]**2-12*coeffs[2]*coeffs[0])**0.5)/(6*coeffs[0])
            x2=(-2*coeffs[1]-(4*coeffs[1]**2-12*coeffs[2]*coeffs[0])**0.5)/(6*coeffs[0])
            if (2*coeffs[1]+6*coeffs[0]*x1 > 2*coeffs[1]+6*coeffs[0]*x2):
                value = x1
            else:
                value = x2
    
    xFit = np.linspace(valueArray[0], valueArray[-1], 1000)
    p = np.poly1d(coeffs)
    yFit = p(xFit)
    
    valueY = p(value)
    
    return (xFit, yFit, value, valueY)


def make_array_loop(varValue, step, numSample):
    """Function to create an array given its middle value and the percentage
    of the extreme.
    
    Parameters
    ----------
    varValue  : float
                variable's value to generate the array
    step      : float
                array step
    numSample : int
                number of sample
    
    Returns
    -------
    varArray  : numpy array
                variable final array
    """

    lowExtreme = varValue
    highExtreme = varValue+step*22
    
    varArray = np.linspace(lowExtreme, highExtreme, numSample)
    
    return varArray