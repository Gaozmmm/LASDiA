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

"""Module containing the main functions used in LASDiA.

The nomenclature and the procedure follow the article:
Eggert et al. 2002 PRB, 65, 174105.

For the functions arguments and the returns I followed this convetion for the notes:
arguments: description - type
returns: description - type.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""

import matplotlib.pyplot as plt

import sys
import os

import numpy as np
import scipy.constants as sc
from scipy import fftpack
from scipy.integrate import simps
from scipy.interpolate import UnivariateSpline

from modules.Utility import *
from modules.UtilityAnalysis import *


def calc_aff(element, Q, elementParameters):
    """Function to calculate the Atomic Form Factor.
    The atomic form factor is calculated with the formula from:
    http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php

    and the parameters from the article:
    Hajdu Acta Cryst. (1972). A28, 250

    Parameters
    ----------
    element           : string
                        chemical element
    Q                 : numpy array
                        momentum transfer (nm^-1)
    elementParameters : dictionary("element": parameters)
                        chemical elements of the sample with their parameters
                        element    : string
                                     chemical element
                        parameters : list
                                     list of the parameters
                                     (Z, a1, b1, a2, b2, a3, b3, a4, b4, c, M, K, L)
    
    Returns
    -------
    f_Q              : numpy array
                       atomic form factor
    """
    
    for elementDict, paramList in elementParameters.items():
        if elementDict == element:
            a1 = paramList[1]
            b1 = paramList[2]
            a2 = paramList[3]
            b2 = paramList[4]
            a3 = paramList[5]
            b3 = paramList[6]
            a4 = paramList[7]
            b4 = paramList[8]
            c = paramList[9]
            
    f1_Q = a1 * np.exp(-b1 * (Q/(4*10*np.pi))**2)
    f2_Q = a2 * np.exp(-b2 * (Q/(4*10*np.pi))**2)
    f3_Q = a3 * np.exp(-b3 * (Q/(4*10*np.pi))**2)
    f4_Q = a4 * np.exp(-b4 * (Q/(4*10*np.pi))**2)

    f_Q = f1_Q + f2_Q + f3_Q + f4_Q + c

    return f_Q


def calc_eeff(elementList, Q, elementParameters):
    """Function to calculate the effective electron Form Factor, fe (eq. 10).

    Parameters
    ----------
    elementList       : dictionary("element": multiplicity)
                        chemical elements of the sample with their multiplicity
                        element      : string
                                       chemical element
                        multiplicity : int
                                       chemical element multiplicity
    Q                 : numpy array
                        momentum transfer (nm^-1)
    elementParameters : dictionary("element": parameters)
                        chemical elements of the sample with their parameters
                        element    : string
                                     chemical element
                        parameters : list
                                     list of the parameters
                                     (Z, a1, b1, a2, b2, a3, b3, a4, b4, c, M, K, L)
    
    Returns
    -------
    fe_Q              : numpy array
                        effective electric form factor
    Ztot              : int
                        total Z number
    """

    fe_Q = np.zeros(Q.size)
    Ztot = 0
    
    for element, multiplicity in elementList.items():
        for elementDict, paramList in elementParameters.items():
            if elementDict == element:
                Ztot += multiplicity * paramList[0]
        fe_Q += multiplicity * calc_aff(element, Q, elementParameters)

    fe_Q /= Ztot

    return (fe_Q, Ztot)


def calc_aff2(element, Q, aff_path):
    """Function to calculate the Atomic Form Factor.
    The atomic form factor is calculated with the formula from:
    http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php

    and the parameters from the article:
    Hajdu Acta Cryst. (1972). A28, 250

    Parameters
    ----------
    element  : string
               chemical element
    Q        : numpy array
               momentum transfer (nm^-1)
    aff_path : string
               path of atomic scattering form factor parameters
    
    Returns
    -------
    f_Q      : numpy array
               atomic form factor
    """

    # open, read and close the parameters file
    # the first line contain the header and it is useless for the calculation
    # file = open("./affParamCEA.txt", "r")
    file = open(aff_path, "r")
    header1 = file.readline()
    lines = file.readlines()
    file.close()

    # scan the lines and when it find the right element, save the parameters in variables
    for line in lines:
        columns = line.split()
        if columns[0] == element:
            a1 = float(columns[1])
            b1 = float(columns[2])
            a2 = float(columns[3])
            b2 = float(columns[4])
            a3 = float(columns[5])
            b3 = float(columns[6])
            a4 = float(columns[7])
            b4 = float(columns[8])
            c = float(columns[9])
            break

    # Calculate the atomic form factor as:
    # f(Q) = f1(Q) + f2(Q) + f3(Q) + f4(Q) + c
    # fi(Q) = ai * exp(-bi * (Q/4pi)^2)

    f1_Q = a1 * np.exp(-b1 * (Q/(4*10*np.pi))**2)
    f2_Q = a2 * np.exp(-b2 * (Q/(4*10*np.pi))**2)
    f3_Q = a3 * np.exp(-b3 * (Q/(4*10*np.pi))**2)
    f4_Q = a4 * np.exp(-b4 * (Q/(4*10*np.pi))**2)

    f_Q = f1_Q + f2_Q + f3_Q + f4_Q + c

    return f_Q


def calc_eeff2(elementList, Q, incoh_path, aff_path):
    """Function to calculate the effective electron Form Factor, fe (eq. 10).

    Parameters
    ----------
    elementList  : dictionary("element": multiplicity)
                   chemical elements of the sample with their multiplicity
                   element      : string
                                  chemical element
                   multiplicity : int
                                  chemical element multiplicity
    Q            : numpy array
                   momentum transfer (nm^-1)
    incoh_path   : string
                   path of incoherent scattered intensities parameters
    aff_path     : string
                   path of atomic scattering form factor parameters
    
    Returns
    -------
    fe_Q         : numpy array
                   effective electric form factor
    Ztot         : int
                   total Z number
    """

    fe_Q = 0
    Ztot = 0

    # file = open("./incohParamCEA.txt", "r")
    file = open(incoh_path, "r")
    header1 = file.readline()
    lines = file.readlines()
    file.close()

    # scan the lines and when it find the right element take the Z number
    for element, multiplicity in elementList.items():
        # print (element, multiplicity)
        for line in lines:
            columns = line.split()
            if columns[0] == element:
                Ztot += multiplicity * float(columns[1])
                break

    # print (Ztot)

    for element, multiplicity in elementList.items():
        fe_Q += multiplicity * calc_aff2(element, Q, aff_path)

    fe_Q /= Ztot

    return (fe_Q, Ztot)


def calc_Kp2(fe_Q, element, Q, aff_path):
    """Function to calculate the average of effective atomic number Kp (eq. 11, 14).

    Parameters
    ----------
    fe_Q    : numpy array
              effective electric form factor
    element : string
              chemical element of the sample
    Q       : numpy array
              momentum transfer (nm^-1)
    
    Returns
    -------
    Kp      : float
              average of effective atomic number
    """

    # effective atomic number
    Kp_Q = calc_aff2(element, Q, aff_path)/fe_Q

    # average effective atomic number
    Kp = np.mean(Kp_Q)

    return Kp


def calc_Iincoh(elementList, Q, elementParameters):
    """Function to calculate the incoherent scattering intensity Iincoh(Q).
    The incoherent scattering intensity is calculated with the formula from the article:
    Hajdu Acta Cryst. (1972). A28, 250

    Parameters
    ----------
    elementList       : dictionary("element": multiplicity)
                        chemical elements of the sample with their multiplicity
                        element      : string
                                       chemical element
                        multiplicity : int
                                       chemical element multiplicity
    Q                 : numpy array
                        momentum transfer (nm^-1)
    
    elementParameters : dictionary("element": parameters)
                        chemical elements of the sample with their parameters
                        element    : string
                                     chemical element
                        parameters : list
                                     list of the parameters
                                     (Z, a1, b1, a2, b2, a3, b3, a4, b4, c, M, K, L)
    
    Returns
    -------
    Iincoh_Q         : numpy array
                       incoherent scattering intensity
    """
    
    Iincoh_Q = np.zeros(Q.size)
    
    # scan the lines and when it find the right element take the Z number
    for element, multiplicity in elementList.items():
        aff = calc_aff(element, Q, elementParameters)
        
        for elementDict, paramList in elementParameters.items():
            if elementDict == element:
                Z = paramList[0]
                M = paramList[10]
                K = paramList[11]
                L = paramList[12]
                break
        
        Iincoh_Q += multiplicity * ((Z - aff**2/Z ) * (1 - M * (np.exp(-K*Q/(4*10*np.pi)) - np.exp(-L*Q/(4*10*np.pi)))))

    return Iincoh_Q


def calc_JQ(Iincoh_Q, Ztot, fe_Q):
    """Function to calculate J(Q) (eq. 35).
    
    Parameters
    ----------
    Iincoh_Q : numpy array
               incoherent scattering intensity
    Ztot     : int
               total Z number
    fe_Q     : numpy array
               effective electric form factor
    
    Returns
    -------
    J_Q      : numpy array
               J(Q)
    """
    
    J_Q = Iincoh_Q/(Ztot**2 * fe_Q**2)
    
    return J_Q


def calc_Kp(fe_Q, element, Q, elementParameters):
    """Function to calculate the average of effective atomic number Kp (eq. 11, 14).

    Parameters
    ----------
    fe_Q              : numpy array
                        effective electric form factor
    element           : string
                        chemical element of the sample
    Q                 : numpy array
                        momentum transfer (nm^-1)
    elementParameters : dictionary("element": parameters)
                        chemical elements of the sample with their parameters
                        element    : string
                                     chemical element
                        parameters : list
                                     list of the parameters
                                     (Z, a1, b1, a2, b2, a3, b3, a4, b4, c, M, K, L)
    
    Returns
    -------
    Kp                : float
                        average of effective atomic number
    """

    # effective atomic number
    Kp_Q = calc_aff(element, Q, elementParameters)/fe_Q

    # average effective atomic number
    Kp = np.mean(Kp_Q)

    return Kp


def calc_Sinf(elementList, fe_Q, Q, Ztot, elementParameters):
    """Function to calculate Sinf (eq. 19).

    Parameters
    ----------
    elementList       : dictionary("element": multiplicity)
                        chemical elements of the sample with their multiplicity
                        element      : string
                                       chemical element
                        multiplicity : int
                                       chemical element multiplicity
    fe_Q              : numpy array
                        effective electric form factor
    Q                 : numpy array
                        momentum transfer (nm^-1)
    Ztot              : int
                        total Z number
    elementParameters : dictionary("element": parameters)
                        chemical elements of the sample with their parameters
                        element    : string
                                     chemical element
                        parameters : list
                                     list of the parameters
                                     (Z, a1, b1, a2, b2, a3, b3, a4, b4, c, M, K, L)
    
    Returns
    -------
    Sinf              : float
                        Sinf
    """
    
    sum_Kp2 = 0
    
    for element, multiplicity in elementList.items():
        sum_Kp2 += multiplicity * calc_Kp(fe_Q, element, Q, elementParameters)**2

    Sinf = sum_Kp2 / Ztot**2

    return Sinf


def calc_IsampleQ(I_Q, s, Ibkg_Q ):
    """Function to calculate the sample scattering intensity Isample(Q) (eq. 28).

    Parameters
    ----------
    I_Q       : numpy array
                measured scattering intensity
    s         : float
                scale factor
    Ibkg_Q     : numpy array
                background scattering intensity
    
    Returns
    -------
    Isample_Q : numpy array
                sample scattering intensity
    """

    Isample_Q = I_Q - s * Ibkg_Q 

    return Isample_Q


def calc_alpha(J_Q, Sinf, Q, Isample_Q, fe_Q, Ztot, rho0):
    """Function to calculate the normalization factor alpha (eq. 34).

    Parameters
    ----------
    J_Q       : numpy array
                J(Q)
    Sinf      : float
                Sinf
    Q         : numpy array
                momentum transfer (nm^-1)
    Isample_Q : numpy array
                sample scattering intensity
    fe_Q      : numpy array
                effective electric form factor
    Ztot      : int 
                total Z number
    rho0      : float
                average atomic density
    
    Returns
    -------
    alpha     : float
                normalization factor
    """

    Integral1 = simps((J_Q + Sinf) * Q**2, Q)
    Integral2 = simps((Isample_Q/fe_Q**2) * Q**2,Q)
    alpha = Ztot**2 * (((-2*np.pi**2*rho0) + Integral1) / Integral2)
    
    return alpha


def calc_Icoh(N, alpha, Isample_Q, Iincoh_Q):
    """Function to calcultate the cohrent scattering intensity Icoh(Q) (eq. 10).

    Parameters
    ----------
    N         : int
                number of atoms
    alpha     : float
                normalization factor
    Isample_Q : numpy array
                sample scattering intensity
    Iincoh_Q  : numpy array
                incoherent scattering intensity
    
    Returns
    -------
    Icoh_Q    : numpy array
                cohrent scattering intensity
    """

    Icoh_Q = N * ((alpha * Isample_Q) - Iincoh_Q)

    return Icoh_Q


def calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, minQ, QmaxIntegrate, maxQ):
    """Function to calculate the structure factor S(Q) (eq. 18) with Igor range.
    
    Parameters
    ----------
    N             : int
                    number of atoms
    Icoh_Q        : numpy array
                    cohrent scattering intensity
    Ztot          : int
                    total Z number
    fe_Q          : numpy array 
                    effective electric form factor
    Sinf          : float
                    Sinf
    Q             : numpy array 
                    momentum transfer (nm^-1)
    minQ          : float
                    minimum Q value
    QmaxIntegrate : float
                    maximum Q value for the intagrations
    maxQ          : float
                    maximum Q value
    
    Returns
    -------
    S_Q           : numpy array
                    structure factor
    """
    
    S_Q = np.zeros(Q.size)
    S_Q[(Q>minQ) & (Q<=QmaxIntegrate)] = Icoh_Q[(Q>minQ) & (Q<=QmaxIntegrate)] / (N * Ztot**2 * fe_Q[(Q>minQ) & (Q<=QmaxIntegrate)]**2)
    S_Q[(Q>QmaxIntegrate) & (Q<=maxQ)] = Sinf
    
    return S_Q


def calc_iQ(S_Q, Sinf):
    """Function to calculate i(Q) (eq. 21).

    Parameters
    ----------
    S_Q  : numpy array
           structure factor
    Sinf : float
           Sinf
    
    Returns
    -------
    i_Q  : numpy array
           i(Q)
    """

    i_Q = S_Q - Sinf

    return i_Q


def calc_QiQ(Q, S_Q, Sinf):
    """Function to calculate Qi(Q) (eq. 7, 20).

    Parameters
    ----------
    Q    : numpy array
           momentum transfer (nm^-1)
    S_Q  : numpy array
           structure factor
    Sinf : float
           Sinf
    
    Returns
    -------
    Qi_Q : numpy array
           Qi(Q)
    """

    Qi_Q = Q*(S_Q - Sinf)

    return Qi_Q


def calc_r(Q):
    """Function to calculate the r value and range used into F(r) calculation.
    
    Parameters
    ----------
    Q : numpy array
        momentum transfer (nm^-1)
    
    Returns
    -------
    r : numpy array
        atomic distance (nm)
    """

    DeltaQ = np.diff(Q)
    meanDeltaQ = np.mean(DeltaQ)
    r = fftpack.fftfreq(Q.size, meanDeltaQ)
    mask = np.where(r>=0)

    return r[mask]


def calc_Fr(r, Q, i_Q):
    """Function to calculate F(r) (eq. 20).
    
    Parameters
    ----------
    r   : numpy array
          atomic distance (nm)
    Q   : numpy array
          momentum transfer (nm^-1)
    i_Q : numpy array 
          i(Q)
    
    Returns
    -------
    F_r  : numpy array
           F(r)
    """

    DeltaQ = np.diff(Q)
    meanDeltaQ = np.mean(DeltaQ)
    rQ = np.outer(r,Q)
    sinrQ = np.sin(rQ)
    F_r = (2.0 / np.pi) * np.sum(Q*i_Q * sinrQ, axis=1) * meanDeltaQ

    return F_r


def calc_Fr2(Q, Qi_Q, QmaxIntegrate):
    """Function to calculate F(r) (eq. 20) with the FFT.

    Parameters
    ----------
    Q             : numpy array
                    momentum transfer (nm^-1)
    Qi_Q          : numpy array 
                    Qi(Q)
    QmaxIntegrate : float
                    maximum Q value for the intagration
    
    Returns
    -------
    F_r  : numpy array
           F(r)
    """

    DeltaQ = np.diff(Q)

    # mask = np.where(Q<=QmaxIntegrate)
    # print(len(mask[0]))
    # print(len(Qi_Q))
    # zeroPad = np.zeros(2**275)
    # Qi_Q = np.concatenate([Qi_Q, zeroPad])

    # F_r2 = np.fft.fft(Qi_Q)
    F_r2 = fftpack.fft(Qi_Q, 2**12)
    F_r2 = np.imag(F_r2)
    F_r2 *= DeltaQ[0] *2/np.pi

    F_r3 = fftpack.dst(Qi_Q, type=2, n=2**11)
    F_r3 *= DeltaQ[0] * 2/np.pi

    # for i in range(len(F_r2)):
        # print(F_r2[i], "----", F_r2imag[i], "----", F_r3[i])

    return (F_r2, F_r3)


def calc_NewDimFFT(newQ, maxQ, QmaxIntegrate, Qi_Q):
    """Function to redefine the array dimension to use the FFT.
    """

    idx, elem = find_nearest(newQ, QmaxIntegrate)
    newDim = 2*2*2**math.ceil(math.log(5*(idx+1))/math.log(2))
    Qi_Q2 = np.resize(Qi_Q, newDim)
    Qi_Q2[idx:] = 0.0

    print(len(Qi_Q2))
    print(newDim)

    newQ2 = np.linspace(np.amin(newQ), maxQ, newDim, endpoint=True)
    print(len(newQ2))
    # DeltaQ = np.diff(newQ)
    # deltaR = 2*np.pi/(DeltaQ[0]*newDim)

    return (newQ2, Qi_Q2)


def calc_SQCorr(F_r, r, Q, Sinf):
    """Function to calculate S(Q) Corr from F(r) optimal.
    
    Parameters
    ----------
    F_r  : numpy array
           F(r)
    r    : numpy array
           atomic distance (nm)
    Q    : numpy array
           momentum transfer (nm^-1)
    Sinf : float
           Sinf
    
    Returns
    -------
    S_Q  : numpy array
           corrected structure factor
    """

    # Qi_Q =  simps(F_r * (np.array(np.sin(np.mat(r).T *  np.mat(Q)))).T, r)
    Deltar = np.diff(r)
    meanDeltar = np.mean(Deltar)
    Qr = np.outer(Q,r)
    sinQr = np.sin(Qr)
    Qi_Q =  np.sum(sinQr * F_r, axis=1) * meanDeltar
    S_Q = Qi_Q/Q +Sinf

    return S_Q
