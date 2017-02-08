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
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

# from modules import Formalism
# from modules import Geometry
# from modules import IgorFunctions
# from modules import KaplowMethod
from modules import MainFunctions
# from modules import Minimization
# from modules import Optimization
from modules import Utility
from modules import UtilityAnalysis


def calc_JQIgor(Iincoh_Q, fe_Q):
    """Function to calculate J(Q) (eq. 35)

    arguments:
    Iincoh_Q: incoherent scattering intensity - array
    Ztot: total Z number - number
    fe_Q: effective electric form factor - array

    returns:
    J_Q: J(Q) - array
    """

    J_Q = Iincoh_Q/(fe_Q**2) # Igor formula

    return J_Q


def calc_alphaIgor(J_Q, Sinf, Q, Isample_Q, fe_Q, Ztot, rho0, index):
    """Function to calculate the normalization factor alpha (eq. 34)

    arguments:
    J_Q: J(Q) - array
    Sinf: Sinf - number
    Q: momentum transfer - array
    Isample_Q: sample scattering intensity - array
    fe_Q: effective electric form factor - array
    Ztot: total Z number - number
    rho0: average atomic density - number
    index: array index of element in the calculation range - array

    returns:
    alpha: normalization factor - number
    """

    Integral1 = simps(Q[index]**2*(J_Q[index] + Sinf*Ztot**2), Q[index])
    Integral2 = simps(Isample_Q[index]*Q[index]**2/fe_Q[index]**2,Q[index])
    alpha = ((-2*np.pi**2*Ztot**2*rho0) + Integral1) / Integral2

    return alpha


def calc_SQIgor(Isample_Q, J_Q, Ztot, fe_Q, Sinf, Q, alpha, min_index, max_index, calculation_index): # Igor formula
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

    S_Q = (alpha * Isample_Q[calculation_index]/fe_Q[calculation_index]**2 - J_Q[calculation_index]) / Ztot**2

    S_Qmin = np.zeros(Q[min_index].size)
    S_Q = np.concatenate([S_Qmin, S_Q])

    S_Qmax = np.zeros(Q[max_index].size)
    S_Qmax.fill(Sinf)
    S_Q = np.concatenate([S_Q, S_Qmax])

    return S_Q


def absorptionIgor(Q):
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
            +muB*rhoB*(dOF[Q>QDAC]-dOI[Q>QDAC])+AlphaDAC*(dOG[Q>QDAC]-dOF[Q>QDAC])))
        AbsRefCalc[Q<QGasket]*=1
    
    return AbsRefCalc


def muFunc(y0,x0,A1,t1,A2,t2):
    """absorption coefficient function"""
    
    x = 33.1952 #x-ray energy (kev)
    mu = y0 + A1*np.exp(-(x-x0)/t1) + A2*np.exp(-(x-x0)/t2) 
    return mu


def FitRemoveGofRPeaks(Q, SsmoothDamp_Q, Sinf, QmaxIntegrate, rintra, Fintra_r,
    iterations, rmin, density, J_Q):
    """Equivalent of Igor function
    """
    
    Qi_Q = Q*MainFunctions.calc_iQ(SsmoothDamp_Q, Sinf)
    # Qi_Q = Q*i_Q
    r, GR = UtilityAnalysis.calc_FFT_QiQ(Q, Qi_Q, QmaxIntegrate)
    Utility.write_file("./GR.txt", r, GR)
    Utility.plot_data(r, GR, "GR", "r", "GR", "GR", "y")
    plt.show()
    
    _, Fintra_r = UtilityAnalysis.rebinning(rintra, Fintra_r, np.amin(Fintra_r), 
        np.amax(Fintra_r), len(GR))
    
    QiQ1 = np.zeros(len(SsmoothDamp_Q))
    idx, _ = UtilityAnalysis.find_nearest(Qi_Q, QmaxIntegrate)
    QiQ1[Q<QmaxIntegrate] = Qi_Q[Q<QmaxIntegrate]
    QiQ1[0] = 0.0
    
    GR1 = GR
    DelG = np.zeros(len(GR))
    Rnn = rmin
    DelG[r<Rnn] = GR1[r<Rnn]-(Fintra_r[r<Rnn]-4*np.pi*r[r<Rnn]*density)
    
    for i in range(iterations):
        Q1, QiQCorr = UtilityAnalysis.calc_IFFT_Fr(r, DelG)
        mask = np.where((Q1>0.0) & (Q1<QmaxIntegrate))
        QiQ1[mask] = QiQ1[mask] - (QiQ1[mask] / 
            (Q1[mask] *(Sinf + J_Q[:len(Q1[mask])])) + 1) * QiQCorr[mask]
        r, GR1 = UtilityAnalysis.calc_FFT_QiQ(Q1, QiQ1, QmaxIntegrate)
        
        DelG = np.zeros(len(GR1))
        DelG[r<Rnn] = GR1[r<Rnn]-(Fintra_r[r<Rnn]-4*np.pi*r[r<Rnn]*density)
        
        Rnn = 0.99*r[np.where(GR1==np.amin(GR1[r>0.95*Rnn]))[0][0]]

    # SQCorr = np.zeros(len(QiQ1))
    # SQCorr[1:] = QiQ1[1:]/Q1[1:]+Sinf
    
    DTemp = DelG
    DTemp = DTemp**2
    chi2 = np.mean(DTemp)
    
    return chi2