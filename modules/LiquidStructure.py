# The MIT License (MIT)

# Copyright (c) 2015 Francesco Devoto

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

"""Set of modules used in the Liquid program

The nomenclature and the procedures follow the article:
Eggert et al. 2002 PRB, 65, 174105

For the functions arguments and the returns I followed this convetion for the notes:
arguments: description - type
returns: description - type
"""

import matplotlib.pyplot as plt

import sys
import os

import numpy as np
import scipy.constants as sc
from scipy import fftpack
from scipy.integrate import simps

# from modules.MainModules import *

def calc_aff(element, Q):
    """Function to calculate the Atomic Form Factor
    The atomic form factor is calculated with the formula from:
    http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
    
    and the parameters from the article:
    Hajdu Acta Cryst. (1972). A28, 250
    
    arguments:
    element: chemical element - string
    Q: transferse moment - array
    
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
    
    f1_Q = a1 * np.exp(-b1 * (Q/(4*np.pi))**2)
    f2_Q = a2 * np.exp(-b2 * (Q/(4*np.pi))**2)
    f3_Q = a3 * np.exp(-b3 * (Q/(4*np.pi))**2)
    f4_Q = a4 * np.exp(-b4 * (Q/(4*np.pi))**2)
    
    f_Q = f1_Q + f2_Q + f3_Q + f4_Q + c
    
    return f_Q
      
    
def calc_eeff(elementList, Q, calc_aff):
    """Function to calculate the effective electron Form Factor, fe
    
    arguments:
    elementList: contains the chemical elements of the sample with their multiplicity - dictionary("element": multiplicity)
        element - string; multiplicity - number
    Q: transferse moment - array
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
    
    
def calc_Iincoh(elementList, Q, calc_aff):
    """ Function to calculate the incoherent scattering intensity Iincoh(Q)
    The incoherent scattering intensity is calculated with the formula and parameters from the article:
    Hajdu Acta Cryst. (1972). A28, 250
    
    arguments:
    elementList: contains the chemical elements of the sample with their multiplicity - dictionary("element": multiplicity)
        element - string; multiplicity - number
    Q: transferse moment - array
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

        Iincoh_Q += multiplicity * ((Z - aff**2/Z ) * (1 - M * (np.exp(-K*Q/(4*np.pi)) - np.exp(-L*Q/(4*np.pi)))))
        
    return Iincoh_Q

    
def calc_JQ(Iincoh_Q, Ztot, fe_Q):
    """Function to calculate J(Q)
    
    arguments:
    Iincoh_Q: incoherent scattering intensity - array
    Ztot: total Z number - number
    fe_Q: effective electric form factor - array
    
    returns:
    J_Q: J(Q) - array
    """
    
    J_Q = Iincoh_Q/(Ztot**2 * fe_Q**2)
    #J_Q = Iincoh_Q/(fe_Q**2)
    
    return J_Q
    
    
def calc_Kp(fe_Q, element, Q, calc_aff):
    """Function to calculate the effective atomic number. Kp_Q, and its average, Kp

    arguments:
    fe_Q: effective electric form factor - array
    element: chemical element of the sample
    Q: transferse moment - array
    calc_aff: function to calculate the atomic form factor - function
    
    returns:
    Kp:  - number
    """

    # effective atomic number
    Kp_Q = calc_aff(element,Q)/fe_Q
    
    # average effective atomic number
    Kp = np.mean(Kp_Q)

    return Kp


def calc_Sinf(elementList, fe_Q, Q, Ztot, calc_Kp, calc_aff):
    """Function to calculate Sinf = sum(Kp^2)/Ztot^2
    
    arguments:
    elementList: contains the chemical elements of the sample with their multiplicity - dictionary("element": multiplicity)
        element - string; multiplicity - number
    fe_Q: effective electric form factor - array
    Q: transferse moment - array
    Ztot: total Z number - number
    calc_Kp: function to calculate the average effective atomic number - function
    calc_aff: function to calculate the atomic form factor - function
    
    returns:
    Sinf: Sinf - number
    """
    
    sum_Kp2 = 0

    for element, multiplicity in elementList.items():
        sum_Kp2 += multiplicity * calc_Kp(fe_Q, element, Q, calc_aff)**2
    
    Sinf = sum_Kp2 / Ztot**2

    return Sinf
