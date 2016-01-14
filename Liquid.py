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

"""Environment where to test the new modules
The nomenclature and the procedures follow the article: Eggert et al. 2002 PRB, 65, 174105
"""

from __future__ import (absolute_import, division, print_function, unicode_literals)
import six

import sys
import os

import scipy.constants as sc
from scipy import fftpack
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

from modules.Utility import *
from modules.LiquidStructure import *
from modules.InterpolateData import *
from modules.Optimization import *
from modules.Minimization import *

# import cmath
# from cmath import exp, pi

if __name__ == '__main__':
    N = 1 # sc.N_A
    
    Q, I_Q = read_file("./data/cea_files/HT2_034T++.chi")
    Qbkg, I_Qbkg = read_file("./data/cea_files/HT2_036T++.chi")
    
    minQ = 3
    maxQ = 109
    QmaxIntegrate = 90
    
    operative_index = np.where((Q>minQ) & (Q<=maxQ))
    operative_Q = Q[operative_index]
    operative_I_Q = I_Q[operative_index]
    operative_I_Qbkg = I_Qbkg[operative_index]
    
    min_index = np.where(Q<=minQ)
    max_index = np.where((Q>QmaxIntegrate) & (Q<=maxQ))
    
    integrate_index = np.where((Q>minQ) & (Q<=QmaxIntegrate))
    integrate_Q = Q[integrate_index]
    integrate_I_Q = I_Q[integrate_index]
    integrate_I_Qbkg = I_Qbkg[integrate_index]
    
    elementList = {"Ar":1}
    s =  0.5
    rho0 = 24.00
    
    # range_flag = "operative"
    # range_flag = "integrate"
    range_flag = ""
    if range_flag == "operative":
        used_Q = operative_Q
        used_I_Q = operative_I_Q
        used_I_Qbkg = operative_I_Qbkg
    elif range_flag == "integrate":
        used_Q = integrate_Q
        used_I_Q = integrate_I_Q
        used_I_Qbkg = integrate_I_Qbkg
    else:
        used_Q = Q
        used_I_Q = I_Q
        used_I_Qbkg = I_Qbkg
    
    # remember the electron unit in atomic form factor!!!
    fe_Q, Ztot = calc_eeff(elementList, used_Q)
    Iincoh_Q = calc_Iincoh(elementList, used_Q)
    J_Q = calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf = calc_Sinf(elementList, fe_Q, used_Q, Ztot)
    Isample_Q = used_I_Q - s * used_I_Qbkg
    alpha = calc_alpha(J_Q, Sinf, used_Q, Isample_Q, fe_Q, Ztot, rho0)
    Icoh_Q = calc_Icoh(N, alpha, Isample_Q, Iincoh_Q)

    # plt.figure(1)
    # plt.plot(used_Q, J_Q)
    # plt.grid()
    # plt.show

    
    print("used_Q = ", used_Q.size)
    print("J_Q = ", J_Q.size)

    if range_flag == "operative":
        S_Q = calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, min_index, max_index)
    elif range_flag == "integrate":
        S_Q = calc_SQ2()
    else:
        S_Q = calc_SQ3(N, Icoh_Q, Ztot, fe_Q)
        
    # used_Q = np.concatenate([Q[min_index], used_Q])
    print("S(Q) = ", S_Q.size)
    print("used_Q = ", used_Q.size)
    
    plt.figure(1)
    plt.plot(used_Q, S_Q)
    plt.grid()
    plt.show
    
    i_Q = calc_iQ(S_Q, Sinf)
    Qi_Q = used_Q * i_Q
    
    print("i(Q) = ", i_Q.size)
    print("Qi(Q) = ", Qi_Q.size)
    
    DeltaQ = np.diff(Q)
    meanDeltaQ = np.mean(DeltaQ)
    r = fftpack.fftfreq(used_Q.size, meanDeltaQ)
    mask = np.where(r>0)
    
    F_r3 = calc_Fr(r[mask], used_Q, i_Q)
    
    plt.figure(2)
    plt.plot(r[mask], F_r3)
    plt.grid()
    plt.show
    
    print("F(r) = ", F_r3.size)
    print("J(Q) = ", J_Q.size)
    
    iteration = 4
    r_cutoff = 0.25
    F_rInt = calc_optimize_Fr(iteration, F_r3, rho0, i_Q, used_Q, Sinf, J_Q, r[mask], r_cutoff)
    
    
    plt.figure(2)
    plt.plot(r[mask], F_rInt)
    # plt.grid()
    plt.show()
