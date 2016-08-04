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

"""Set of modules used in LASDiA to calculate the relations among the different set of Structure Factors.

The formalisms used are Ashcroft-Langreth and Faber-Ziman.

The nomenclature and the procedure follow the book:
Waseda - The Structure of Non-Crystalline Materials.

For the functions arguments and the returns I followed this convetion for the notes:
argument: description - type
return: description - type.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by an underscore: f(x) -> f_x
otherwise it is just the name.
"""


import sys
import os
import math
import numpy as np
import scipy.constants as sc
from scipy import fftpack
from scipy.integrate import simps
from scipy.stats import chisquare
import matplotlib.pyplot as plt

from modules import MainFunctions


def calc_alphaFZ(Q, Isample_Q, Iincoh_Q, rho0, eeff_squared_mean, eeff_mean_squared):
    """Function to calcultate alpha for the FZ formalism.
    
    """
    
    Integral1 = simps((Iincoh_Q + (eeff_squared_mean/eeff_mean_squared)) * Q**2, Q)
    Integral2 = simps((Isample_Q/eeff_mean_squared) * Q**2,Q)
    alpha = ((-2*np.pi**2*rho0) + Integral1) / Integral2
    
    return alpha


def calc_eeff_squared_mean(numAtoms, elementList, Q, elementParameters):
    """Function to calculate the squared mean effective electron Form Factor, <f^2> (eq. FZ-4).
        
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

    eeff_squared_mean = np.zeros(Q.size)
    
    for element, multiplicity in elementList.items():
        eeff_squared_mean += multiplicity * MainFunctions.calc_aff(element, Q, elementParameters)**2
    
    eeff_squared_mean /= numAtoms
    
    return eeff_squared_mean


def calc_eeff_mean_squared(numAtoms, elementList, Q, elementParameters):
    """Function to calculate the mean squared effective electron Form Factor, <f^2> (eq. FZ-5).
        
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

    eeff_mean_squared = np.zeros(Q.size)
    
    for element, multiplicity in elementList.items():
        eeff_mean_squared += multiplicity * MainFunctions.calc_aff(element, Q, elementParameters)
    
    eeff_mean_squared /= numAtoms
    
    eeff_mean_squared = eeff_mean_squared**2
    
    return eeff_mean_squared


def calc_S_QFZ(Q, Icoh_Q, eeff_squared_mean, eeff_mean_squared, minQ, QmaxIntegrate, maxQ):
    """Function to calculate S(Q) with Faber-Ziman formalism.
    
    Parameters
    ----------
    Icoh_Q            : numpy array (nm)
                        coherent scattering intensity
    Q                 : numpy array
                        momentum transfer (nm)
    Sinf              : float
                        Sinf
    min_index         : numpy array
                        array index of element with Q<=minQ
    max_index         : numpy array
                        array index of element with Q>QmaxIntegrate & Q<=maxQ
    calculation_index : numpy array
                        range where S(Q) is calculated
    
    
    Returns
    -------
    S_Q               : numpy array
                        structure factor in FZ formalism
    """
    
    S_Q = np.zeros(Q.size)
    S_Q[(Q>minQ) & (Q<=QmaxIntegrate)] = (Icoh_Q[(Q>minQ) & (Q<=QmaxIntegrate)] - \
        (eeff_squared_mean[(Q>minQ) & (Q<=QmaxIntegrate)] - eeff_mean_squared[(Q>minQ) & (Q<=QmaxIntegrate)])) / (eeff_mean_squared[(Q>minQ) & (Q<=QmaxIntegrate)])
    S_Q[(Q>QmaxIntegrate) & (Q<=maxQ)] = 1
    
    return S_Q


def calc_S_QAL(Q, Icoh_Q, Sinf, eeff_squared_mean, minQ, QmaxIntegrate, maxQ):
    """Function to calculate S(Q) with Ashcroft-Langreth formalism.
    
    Parameters
    ----------
    Icoh_Q            : numpy array (nm)
                        coherent scattering intensity
    Q                 : numpy array
                        momentum transfer (nm)
    Sinf              : float
                        Sinf
    min_index         : numpy array
                        array index of element with Q<=minQ
    max_index         : numpy array
                        array index of element with Q>QmaxIntegrate & Q<=maxQ
    calculation_index : numpy array
                        range where S(Q) is calculated
    
    
    Returns
    -------
    S_Q               : numpy array
                        structure factor in FZ formalism
    """
    
    S_Q = np.zeros(Q.size)
    S_Q[(Q>minQ) & (Q<=QmaxIntegrate)] = Icoh_Q[(Q>minQ) & (Q<=QmaxIntegrate)] / eeff_squared_mean[(Q>minQ) & (Q<=QmaxIntegrate)]
    S_Q[(Q>QmaxIntegrate) & (Q<=maxQ)] = Sinf
    
    return S_Q


# Functions to jump from a formalism to another, need to be fixed!
def ALtoFZ(N1, N2, S11, S12, S22):
    """Function to transform the AL to FZ (eq 1.2.17)

    arguments:
    N1:
    N2:
    S11:
    S12:
    S22:

    returns:
    a11:
    a12:
    a22:
    """

    c1 = N1/(N1+N2)
    c2 = N2/(N1+N2)

    a11_Q = (S11 - c2)/c1
    a12_Q = S12/math.sqrt(c1*c2) + 1
    a22_Q = (S22 - c1)/c2

    return (a11_Q, a12_Q, a22_Q)


def FZtoAL(N1, N2, a11, a12, a22):
    """Function to transform the AL to FZ (eq 1.2.17)

    arguments:
    N1:
    N2:
    a11:
    a12:
    a22:

    returns:
    S11:
    S12:
    S22:
    """

    c1 = N1/(N1+N2)
    c2 = N2/(N1+N2)

    S11_Q = 1 + c1*(a11_Q - 1)
    S12_Q = math.sqrt(c1*c2) * (a12_Q - 1)
    S11_Q = 1 + c2*(a22_Q - 1)

    return (S11_Q, S12_Q, S22_Q)