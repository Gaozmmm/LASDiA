# The MIT License (MIT)

# Copyright (c) 2016 Francesco Devoto

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


def calc_aff(element, Q):
    """Function to calculate the Atomic Form Factor.
    The atomic form factor is calculated with the formula from:
    http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
    
    and the parameters from the article:
    Hajdu Acta Cryst. (1972). A28, 250
    
    arguments:
    element: chemical element - string
    Q: momentum transfer - array
    
    returns:
    f_Q: atomic form factor - array
    """
    
    # open, read and close the parameters file
    # the first line contain the header and it is useless for the calculation
    file = open("./affParamCEA.txt", "r")
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
    
    
def calc_eeff(elementList, Q):
    """Function to calculate the effective electron Form Factor, fe (eq. 10)
    
    arguments:
    elementList: contains the chemical elements of the sample with their multiplicity - dictionary("element": multiplicity)
        element - string; multiplicity - number
    Q: momentum transfer - array
    
    functions used into the calculation:
    calc_aff: function to calculate the atomic form factor - function
    
    returns:
    fe_Q: effective electric form factor - array
    Ztot: total Z number - number
    """
    
    fe_Q = 0
    Ztot = 0

    file = open("./incohParamCEA.txt", "r")
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
        fe_Q += multiplicity * calc_aff(element, Q)
    
    fe_Q /= Ztot

    return (fe_Q, Ztot)
    
    
def calc_Iincoh(elementList, Q):
    """Function to calculate the incoherent scattering intensity Iincoh(Q)
    The incoherent scattering intensity is calculated with the formula from the article:
    Hajdu Acta Cryst. (1972). A28, 250
    
    arguments:
    elementList: contains the chemical elements of the sample with their multiplicity - dictionary("element": multiplicity)
        element - string; multiplicity - number
    Q: momentum transfer - array
    
    functions used into the calculation:
    calc_aff: function to calculate the atomic form factor - function
    
    returns:
    Iincoh_Q: incoherent scattering intensity - array
    """
    
    file = open("./incohParamCEA.txt", "r")
    header1 = file.readline()
    lines = file.readlines()
    file.close()
    
    Iincoh_Q = 0
    
    # scan the lines and when it find the right element take the Z number
    for element, multiplicity in elementList.items():
        aff = calc_aff(element, Q)
        for line in lines:
            columns = line.split()
            if columns[0] == element:
                Z = float(columns[1])
                M = float(columns[2])
                K = float(columns[3])
                L = float(columns[4])
                break

        Iincoh_Q += multiplicity * ((Z - aff**2/Z ) * (1 - M * (np.exp(-K*Q/(4*10*np.pi)) - np.exp(-L*Q/(4*10*np.pi)))))
        
    return Iincoh_Q

    
def calc_JQ(Iincoh_Q, Ztot, fe_Q):
    """Function to calculate J(Q) (eq. 35)
    
    arguments:
    Iincoh_Q: incoherent scattering intensity - array
    Ztot: total Z number - number
    fe_Q: effective electric form factor - array
    
    returns:
    J_Q: J(Q) - array
    """
    
    J_Q = Iincoh_Q/(Ztot**2 * fe_Q**2)
    
    return J_Q
    
    
def calc_Kp(fe_Q, element, Q):
    """Function to calculate the average of effective atomic number Kp (eq. 11, 14)

    arguments:
    fe_Q: effective electric form factor - array
    element: chemical element of the sample - string
    Q: momentum transfer - array
    
    functions used into the calculation:
    calc_aff: function to calculate the atomic form factor - function
    
    returns:
    Kp: average of effective atomic number - number
    """

    # effective atomic number
    Kp_Q = calc_aff(element,Q)/fe_Q
    
    # average effective atomic number
    Kp = np.mean(Kp_Q)

    return Kp


def calc_Sinf(elementList, fe_Q, Q, Ztot):
    """Function to calculate Sinf (eq. 19)
    
    arguments:
    elementList: contains the chemical elements of the sample with their multiplicity - dictionary("element": multiplicity)
        element - string; multiplicity - number
    fe_Q: effective electric form factor - array
    Q: momentum transfer - array
    Ztot: total Z number - number
    
    functions used into the calculation:
    calc_Kp: function to calculate the average effective atomic number - function
    calc_aff: function to calculate the atomic form factor - function
    
    returns:
    Sinf: Sinf - number
    """
    
    sum_Kp2 = 0

    for element, multiplicity in elementList.items():
        sum_Kp2 += multiplicity * calc_Kp(fe_Q, element, Q)**2
    
    Sinf = sum_Kp2 / Ztot**2

    return Sinf
    
    
def calc_IsampleQ(I_Q, s, I_Qbkg):
    """Function to calculate the sample scattering intensity Isample(Q) (eq. 28)
    
    arguments:
    I_Q: measured scattering intensity - array
    s: scale factor - number
    I_Qbkg: background scattering intensity - array
    
    returns:
    Isample_Q: sample scattering intensity - array
    """
    
    Isample_Q = I_Q - s * I_Qbkg
    
    return Isample_Q
    
    
def calc_alpha(J_Q, Sinf, Q, Isample_Q, fe_Q, Ztot, rho0, index):
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
    
    Integral1 = simps((J_Q[index] + Sinf) * Q[index]**2, Q[index])
    Integral2 = simps((Isample_Q[index]/fe_Q[index]**2) * Q[index]**2,Q[index])
    alpha = Ztot**2 * (((-2*np.pi**2*rho0) + Integral1) / Integral2)
    
    # DeltaQ = np.diff(Q)
    # meanDeltaQ = np.mean(DeltaQ)
    # Int1 = np.sum((J_Q[index] + Sinf) * Q[index]**2) * meanDeltaQ
    # Int2 = np.sum( (Isample_Q[index]/fe_Q[index]**2) * Q[index]**2  ) * meanDeltaQ
    # alpha = Ztot**2 * (((-2*np.pi**2*rho0) + Int1) / Int2)

    return alpha
    
    
def calc_Icoh(N, alpha, Isample_Q, Iincoh_Q):
    """Function to calcultate the cohrent scattering intensity Icoh(Q) (eq. 10)
    
    arguments:
    N: number of atoms - number
    alpha: normalization factor - number
    Isample_Q: sample scattering intensity - array
    Iincoh_Q: incoherent scattering intensity - array
    
    returns:
    Icoh_Q: cohrent scattering intensity - array
    """
    
    Icoh_Q = N * ((alpha * Isample_Q) - Iincoh_Q)
    
    return Icoh_Q
    
    
def calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, max_index, integration_index):
    """Function to calculate the structure factor S(Q) (eq. 18) with Igor range.
    This function doesn't set the value 0 for Q<minQ!!!

    arguments:
    N: number of atoms - number
    Icoh_Q: cohrent scattering intensity - array
    Ztot: total Z number - number
    fe_Q: effective electric form factor - array
    Sinf: Sinf - number
    Q: momentum transfer - array
    max_index: array index of element with Q>QmaxIntegrate & Q<=maxQ - array
    integration_index: array index of element in the integration range Q<=QmaxIntegrate - array
    
    returns:
    S_Q: structure factor - array
    """
    
    S_Q = Icoh_Q[integration_index] / (N * Ztot**2 * fe_Q[integration_index]**2)
    
    S_Qmax = np.zeros(Q[max_index].size)
    S_Qmax.fill(Sinf)
    S_Q = np.concatenate([S_Q, S_Qmax])
    
    return S_Q
    
    
def calc_SQdamp(S_Q, Q, Sinf, QmaxIntegrate, damping_factor):
    """
    """
    
    exponent_factor = damping_factor / QmaxIntegrate**2
    damp_Q = np.exp(-exponent_factor * Q**2)
    
    S_Qdamp = (damp_Q * (S_Q - Sinf)) + Sinf
    
    return S_Qdamp
    
    
def calc_iQ(S_Q, Sinf):
    """Function to calculate i(Q) (eq. 21)
    
    arguments:
    S_Q: structure factor - array
    Sinf: Sinf - number
    
    returns:
    i_Q: i(Q) - array
    """
    
    i_Q = S_Q - Sinf
    
    return i_Q
    
    
def calc_QiQ(Q, S_Q, Sinf):
    """Function to calculate Qi(Q) (eq. 7, 20)
    
    arguments:
    Q: momentum transfer - array
    S_Q: structure factor - array
    Sinf: Sinf - number
    
    returns:
    Qi_Q: Qi(Q) - array
    """
    
    Qi_Q = Q*(S_Q - Sinf)
    
    return Qi_Q
    

def calc_Fr(r, Q, i_Q, Qi_Q):
    """Function to calculate F(r) (eq. 20)
    
    arguments:
    r: radius - array
    Q: momentum transfer - array
    Qi_Q: Qi(Q) - array
    
    returns:
    F_r: F(r) - array
    """
    
    DeltaQ = np.diff(Q)
    meanDeltaQ = np.mean(DeltaQ)
    rQ = np.outer(r,Q)
    sinrQ = np.sin(rQ)
    F_r = (2.0 / np.pi) * np.sum(Q*i_Q * sinrQ, axis=1) * meanDeltaQ
    
    F_r2 = np.fft.fft(Qi_Q)
    F_r3 = np.fft.rfft(Qi_Q)
    
    F_r2 = np.imag(F_r2) * DeltaQ[0] *2/np.pi
    
    return (F_r, F_r2, F_r3)
    
    
# def calc_gr(r, F_r, rho0):
    # """Function to calculate g(r)
    
    # """
    
    # g_r = 1 + F_r / (4 * np.pi * r * rho0)
    
    # return g_r
    
    
def calc_SQ_F(F_r, r, Q, Sinf):
    """Function to calculate S(Q) from F(r)

    """

    # Qi_Q =  simps(F_r * (np.array(np.sin(np.mat(r).T *  np.mat(Q)))).T, r)
    Deltar = np.diff(r)
    meanDeltar = np.mean(Deltar)
    Qr = np.outer(Q,r)
    sinQr = np.sin(Qr)
    Qi_Q =  np.sum(sinQr * F_r, axis=1) * meanDeltar
    S_Q = Qi_Q/Q +Sinf

    return S_Q