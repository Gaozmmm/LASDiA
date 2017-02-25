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

"""Module containing the main functions used in LASDiA.

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
from scipy import fftpack
from scipy.integrate import simps
import math

from modules import Utility
from modules import UtilityAnalysis


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


def calc_Iincoh(elementList, Q, elementParameters):
    """Function to calculate the incoherent scattering intensity Iincoh(Q).
    The incoherent scattering intensity is calculated with the formula from the article:
    Hajdu Acta Cryst. (1972). A28, 250.
    For some atoms (as Bi) the incoherent scattering intensity is calculated with the
    formula from the article:
    Palinkas Acta Cryst. (1973). A29, 10.
    To be consistent with the other formula, the Palinkas' parameters a, b, and c are
    renamed M, K, and L, respectly.
    
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
    
    for element, multiplicity in elementList.items():
        aff = calc_aff(element, Q, elementParameters)
        
        for elementDict, paramList in elementParameters.items():
            if elementDict == element:
                Z = paramList[0]
                M = paramList[10]
                K = paramList[11]
                L = paramList[12]
                break
        
        if Z >= 37:
            Iincoh_Q += multiplicity * (Z*(1 - M/(1+K*Q)**L))
        else:
            Iincoh_Q += multiplicity * ((Z - aff**2/Z) * (1 - M * (np.exp(-K*Q/(4*10*np.pi)) \
                - np.exp(-L*Q/(4*10*np.pi)))))
    
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
                        value of S(Q) for Q->inf
    """
    
    sum_Kp2 = 0

    for element, multiplicity in elementList.items():
        sum_Kp2 += multiplicity * calc_Kp(fe_Q, element, Q, elementParameters)**2

    Sinf = sum_Kp2 / Ztot**2

    return Sinf


def calc_IsampleQ(I_Q, scale_factor, Ibkg_Q):
    """Function to calculate the sample scattering intensity Isample(Q) (eq. 28).
    
    Parameters
    ----------
    I_Q          : numpy array
                   measured scattering intensity
    scale_factor : float
                   scale factor
    Ibkg_Q       : numpy array
                   background scattering intensity
    
    Returns
    -------
    Isample_Q    : numpy array
                   sample scattering intensity
    """
    
    Isample_Q = I_Q - scale_factor * Ibkg_Q
    
    return Isample_Q


def calc_alpha(J_Q, Sinf, Q, Isample_Q, fe_Q, Ztot, rho0):
    """Function to calculate the normalization factor alpha (eq. 34).
    
    Parameters
    ----------
    J_Q       : numpy array
                J(Q)
    Sinf      : float
                value of S(Q) for Q->inf
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
    Integral2 = simps((Isample_Q/fe_Q**2) * Q**2, Q)
    alpha = Ztot**2 * (((-2*np.pi**2*rho0) + Integral1) / Integral2)
    
    return alpha


def calc_Icoh(alpha, Isample_Q, Iincoh_Q):
    """Function to calcultate the cohrent scattering intensity Icoh(Q) (eq. 27).
    
    Parameters
    ----------
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
    
    Icoh_Q = (alpha * Isample_Q) - Iincoh_Q
    
    return Icoh_Q


def calc_SQ(Icoh_Q, Ztot, fe_Q, Sinf, Q, minQ, QmaxIntegrate, maxQ):
    """Function to calculate the structure factor S(Q) (eq. 18) with Igor range.
    
    Parameters
    ----------
    Icoh_Q        : numpy array
                    cohrent scattering intensity
    Ztot          : int
                    total Z number
    fe_Q          : numpy array
                    effective electric form factor
    Sinf          : float
                    value of S(Q) for Q->inf
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
    S_Q[(Q > minQ) & (Q <= QmaxIntegrate)] = Icoh_Q[(Q > minQ) & (Q <= QmaxIntegrate)] \
        / (Ztot**2 * fe_Q[(Q > minQ) & (Q <= QmaxIntegrate)]**2)
    S_Q[Q > QmaxIntegrate] = Sinf
    
    return S_Q


def calc_iQ(S_Q, Sinf):
    """Function to calculate i(Q) (eq. 21).
    
    Parameters
    ----------
    S_Q  : numpy array
           structure factor
    Sinf : float
           value of S(Q) for Q->inf
    
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
    i_Q  : numpy array
           i(Q)
    
    Returns
    -------
    Qi_Q : numpy array
           Qi(Q)
    """
    
    i_Q = S_Q - Sinf
    Qi_Q = Q*i_Q
    
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


def calc_Fr(Q, Qi_Q):
    """Function to calculate F(r) (eq. 20) with the FFT.

    Parameters
    ----------
    Q    : numpy array
           momentum transfer (nm^-1)
    Qi_Q : numpy array
           Qi(Q)
    
    Returns
    -------
    r    : numpy array
           atomic distance (nm)
    F_r  : numpy array
           F(r) distribution function
    """

    meanDeltaQ = np.mean(np.diff(Q))
    r = fftpack.fftfreq(Q.size, meanDeltaQ)
    mask = np.where(r>=0)
    
    rQ = np.outer(r[mask], Q)
    sinrQ = np.sin(rQ)
    
    F_r = (2.0 / np.pi) * simps(Qi_Q * sinrQ, Q)
    
    return (r[mask], F_r)
    
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
           value of S(Q) for Q->inf

    Returns
    -------
    S_Q  : numpy array
           corrected structure factor
    """

    # Qi_Q =  simps(F_r * (np.array(np.sin(np.mat(r).T *  np.mat(Q)))).T, r)
    Deltar = np.diff(r)
    meanDeltar = np.mean(Deltar)
    Qr = np.outer(Q, r)
    sinQr = np.sin(Qr)
    Qi_Q = np.sum(sinQr * F_r, axis=1) * meanDeltar
    S_Q = Qi_Q / Q + Sinf
    S_Q[0] = 0.0

    return S_Q
