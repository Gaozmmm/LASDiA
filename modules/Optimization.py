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

For the variables name I used this convention:
if the variable symbolizes a function, its argument is preceded by an underscore: f(x) -> f_x
"""

import matplotlib.pyplot as plt

import sys
import os

import numpy as np
import scipy.constants as sc
from scipy import fftpack
from scipy.integrate import simps

from lmfit import minimize, Parameters

from modules.LiquidStructure import *

def calc_Fintra():
    """To implemente!!!
    eq. 42
    """
    
    Fintra_r = 0
    
    return Fintra_r
    
    
def calc_deltaFr(F_r, Fintra_r, rho0):
    """Function to calculate deltaF(r) (eq. 44, 48)
    
    arguments:
    F(r): - array
    Fintra_r: intramolecular contribution of F(r) - array
    rho0: average atomic density - number
    
    return:
    deltaF_r: - array
    """
    
    deltaF_r = F_r - (Fintra_r - 4*np.pi*rho0)
    
    return deltaF_r
    

def calc_iQi(i_Q, Q, Sinf, J_Q, deltaF_r, r, r_cutoff):
    """Function to calculate the i-th iteration of i(Q) (eq. 49)
    
    arguments:
    i_Q:
    Q:
    Sinf:
    J_Q:
    deltaF_r:
    r_cutoff: - number
    
    returns:
    i_Qi: i-th iteration of i(Q) - array
    """
    
    rint = np.linspace(0, r_cutoff, r.size)
    # print("rint ", rint.size)
    # print("r ", r.size)
    # print("Q ", Q.size)
    # print("deltaF_r ", deltaF_r.size)
    
    integral = simps(deltaF_r * (np.array(np.sin(np.mat(r).T *  np.mat(Q)))).T, rint) #simps(deltaF_r * np.sin(np.mat(r)*Q), rint)
    
    i_Qi = i_Q - ( 1/Q * ( i_Q / (Sinf + J_Q) + 1)) * integral
    
    return i_Qi

    
def calc_deltaAlpha(i_Q, i_Qi, S_Q):
    """
    """
    
    deltaAlpha = (i_Q - i_Qi) / S_Q
    
    return deltaAlpha


def calc_alphai(alpha, deltaAlpha):
    """
    
    """
    
    alphai = alpha * (1 + deltaAlpha)
    
    return alphai
    
    
def calc_SQi():
 
 

def calc_optimize_Fr(iteration, i_Q, F_r, Fintra_r, rho0, Q, Sinf, J_Q, r, r_cutoff, S_Q, alpha, N, Isample_Q, Iincoh_Q, Ztot, fe_Q):
    """Function to calculate the optimization
    
    """
   
    for i in range(iteration):
        # i_QTemp = i_Q
        # S_QTemp = S_Q
        # deltaF = calc_deltaFr(F_r, Fintra_r, rho0)
        # i_Q = calc_iQi(i_Q, Q, Sinf, J_Q, deltaF, r, r_cutoff)
        
        # deltaAlpha = calc_deltaAlpha(i_QTemp, i_Q, S_QTemp)
        # alpha = calc_alphai(alpha, deltaAlpha)
        # Icoh_Q = calc_Icoh(N, alpha, Isample_Q, Iincoh_Q)
        # S_Q = calc_SQ(N, Icoh_Q, Ztot, fe_Q)
        # F_r = calc_Fr(r, Q, i_Q)
        deltaF = calc_deltaFr(F_r, Fintra_r, rho0)
        i_Q = calc_iQi(i_Q, Q, Sinf, J_Q, deltaF, r, r_cutoff)
        F_r = calc_Fr(r, Q, i_Q)
    
    optF_r = F_r
    return optF_r
    
    
    
    
