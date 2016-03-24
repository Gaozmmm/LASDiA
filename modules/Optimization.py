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

"""Module containing the functions used for the F(r) optimization.

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
import time
from scipy import fftpack
from scipy.integrate import simps

from modules.MainFunctions import *
from modules.Utility import *
    
    
def calc_iintra(Q, max_index, elementList, path):
    """Function to calculate the intramolecular contribution of i(Q) (eq. 41)
    
    """
    
    fe_Q, Ztot = calc_eeff(elementList, Q)
    
    numAtoms, element, x, y, z = read_xyz_file(path)
    iintra_Q = np.zeros(Q.size)
    sinpq = np.zeros(Q.size)
    
    for ielem in range(len(element)):
        for jelem in range(len(element)):
            if ielem != jelem:
                KK = calc_Kp(fe_Q, element[ielem], Q) * calc_Kp(fe_Q, element[jelem], Q)
                d = calc_distMol(x[ielem], y[ielem], z[ielem], x[jelem], y[jelem], z[jelem])
                if d != 0.0:
                    iintra_Q += KK * np.sin(d*Q) / (d*Q)
                    iintra_Q[Q==0.0] = KK
    
    iintra_Q[max_index] = 0.0
    iintra_Q /= Ztot**2
    
    return (iintra_Q, fe_Q)
    
    
def calc_Fintra(r, Q, QmaxIntegrate):
    """Function to calculate the intramolecular contribution of F(r) (eq. 42)
    
    To implemente!!! -> For now just for CO2!!!
    """
    
    # Fintra_r = np.zeros(r.size)
    
    dCO = 0.1165 # nm
    dOO = 2 * dCO
    
    elementList = {"C":1,"O":2}
    fe_Q, Ztot = calc_eeff(elementList, Q)
    KC = calc_Kp(fe_Q, "C", Q)
    KO = calc_Kp(fe_Q, "O", Q)
    
    constCO = 4/(np.pi * Ztot**2 * dCO)
    constOO = 2/(np.pi * Ztot**2 * dOO)
    
    Fintra_r_CO = constCO * KC * KO * \
        ((np.sin((r - dCO)*QmaxIntegrate)) / (r - dCO) - (np.sin((r + dCO)*QmaxIntegrate)) / (r + dCO))
    
    Fintra_r_OO = constOO * KO * KO * \
        ((np.sin((r - dOO)*QmaxIntegrate)) / (r - dOO) - (np.sin((r + dOO)*QmaxIntegrate)) / (r + dOO))
    
    Fintra_r = Fintra_r_CO + Fintra_r_OO
    
    return Fintra_r
    
    
def calc_deltaFr(F_r, Fintra_r, r, rho0):
    """Function to calculate deltaF(r) (eq. 44, 48)
    
    arguments:
    F_r: F(r) - array
    Fintra_r: intramolecular contribution of F(r) - array
    r: radius - array
    rho0: average atomic density - number
    
    return:
    deltaF_r: deltaF(r) - array
    """
    
    deltaF_r = F_r - (Fintra_r - 4*np.pi*r*rho0)
    
    return deltaF_r
    

def calc_iQi(i_Q, Q, Sinf, J_Q, deltaF_r, r, rmin):
    """Function to calculate the i-th iteration of i(Q) (eq. 46, 49)
    
    arguments:
    i_Q: i(Q) - array
    Q: momentum transfer - array
    Sinf: Sinf - number
    J_Q: J(Q) - array
    deltaF_r: deltaF(r) - array
    rmin: value of r cutoff - number
    
    returns:
    i_Qi: i-th iteration of i(Q) - array
    """
    
    mask = np.where(r < rmin)
    rInt = r[mask]
    deltaF_rInt = deltaF_r[mask]
    
    Deltar = np.diff(rInt)
    meanDeltar = np.mean(Deltar)
    Qr = np.outer(Q, rInt)
    sinQr = np.sin(Qr)
    integral = np.sum(deltaF_rInt * sinQr, axis=1) * meanDeltar
    
    i_Qi = i_Q - ( 1/Q * ( i_Q / (Sinf + J_Q) + 1)) * integral
    
    return i_Qi

    
def calc_optimize_Fr(iteration, F_r, Fintra_r, rho0, i_Q, Q, Sinf, J_Q, r, rmin):
    """Function to calculate the F(r) optimization (eq 47, 48, 49)
    
    arguments:
    iteration: number of iteration - number
    F_r: initial value of F(r) - array
    rho0: average atomic density - number
    i_Q: i(Q) - array
    Q: momentum transfer - array
    Sinf: Sinf - number
    J_Q: J(Q) - array
    r: radius - array
    rmin: value of r cutoff - number
    
    returns:
    F_r: optimazed F(r) - array
    """
    
    # plt.ion()
    # plt.figure('F_r')
    # plt.plot(r, F_r, label='F(r)')
    # plt.xlabel('r (nm)')
    # plt.ylabel('F(r)')
    # plt.legend()
    # plt.grid()
    
    for i in range(iteration):
        deltaF_r = calc_deltaFr(F_r, Fintra_r, r, rho0)
        i_Q = calc_iQi(i_Q, Q, Sinf, J_Q, deltaF_r, r, rmin)
        i_Q[0] = 0.0
        F_r = calc_Fr(r, Q, Q*i_Q)
        
        # j = i+1
        # plt.figure('F_r')
        # plt.plot(r, F_r, label='%s iteration F(r)' %j)
        # plt.legend()
        # plt.draw()
        
        # time.sleep(1.0)
    
    # plt.ioff()
    # plt.show()
    
    return F_r
    
    
    
    
