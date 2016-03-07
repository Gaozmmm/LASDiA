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

"""Set of modules used in LASDiA to calculate all the used functions.

The nomenclature and the procedure follow the article:
Eggert et al. 2002 PRB, 65, 174105.

For the functions arguments and the returns I followed this convetion for the notes:
arguments: description - type
returns: description - type.

For the variables name I used this convention:
if the variable symbolizes a function, its argument is preceded by an underscore: f(x) -> f_x
otherwise it is just the name.
"""

import matplotlib.pyplot as plt

import sys
import os

import numpy as np
import scipy.constants as sc
from scipy import fftpack
from scipy.integrate import simps
from scipy.interpolate import UnivariateSpline


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
    
    
def calc_SQIgor2(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, max_index, integration_index): # Igor formula
    """Function to calculate the structure factor S(Q) (eq. 18) with Igor range
    This function doesn't set the value 0 for Q<minQ!!!

    arguments:
    N: number of atoms - number
    Icoh: cohrent scattering intensity - array
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
    