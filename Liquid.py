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
    
    range_flag = "operative"
    # range_flag = "integrate"
    # range_flag = ""
    # if range_flag == "operative":
        # used_Q = operative_Q
        # used_I_Q = operative_I_Q
        # used_I_Qbkg = operative_I_Qbkg
    # elif range_flag == "integrate":
        # used_Q = integrate_Q
        # used_I_Q = integrate_I_Q
        # used_I_Qbkg = integrate_I_Qbkg
    # else:
        # used_Q = Q
        # used_I_Q = I_Q
        # used_I_Qbkg = I_Qbkg
    
    # remember the electron unit in atomic form factor!!!
    
    fe_Q1, Ztot1 = calc_eeff(elementList, Q)
    fe_Q2, Ztot2 = calc_eeff(elementList, operative_Q)
    fe_Q3, Ztot3 = calc_eeff(elementList, integrate_Q)
    print(Ztot1, Ztot2, Ztot3)
    plt.figure(1)
    plt.plot(Q, fe_Q1)
    plt.plot(operative_Q, fe_Q2)
    plt.plot(integrate_Q, fe_Q3)
    plt.grid()
    plt.show
    
    
    Iincoh_Q1 = calc_Iincoh(elementList, Q)
    Iincoh_Q2 = calc_Iincoh(elementList, operative_Q)
    Iincoh_Q3 = calc_Iincoh(elementList, integrate_Q)
    plt.figure(2)
    plt.plot(Q, Iincoh_Q1)
    plt.plot(operative_Q, Iincoh_Q2)
    plt.plot(integrate_Q, Iincoh_Q3)
    plt.grid()
    plt.show
    
    J_Q1 = calc_JQ(Iincoh_Q1, Ztot1, fe_Q1)
    J_Q2 = calc_JQ(Iincoh_Q2, Ztot2, fe_Q2)
    J_Q3 = calc_JQ(Iincoh_Q3, Ztot3, fe_Q3)
    plt.figure(3)
    plt.plot(Q, J_Q1)
    plt.plot(operative_Q, J_Q2)
    plt.plot(integrate_Q, J_Q3)
    plt.grid()
    plt.show
    
    
    Sinf1 = calc_Sinf(elementList, fe_Q1, Q, Ztot1)
    Sinf2 = calc_Sinf(elementList, fe_Q2, operative_Q, Ztot2)
    Sinf3 = calc_Sinf(elementList, fe_Q3, integrate_Q, Ztot3)
    print(Sinf1, Sinf2, Sinf3)
    
    Isample_Q1 = I_Q - s * I_Qbkg
    Isample_Q2 = operative_I_Q - s * operative_I_Qbkg
    Isample_Q3 = integrate_I_Q - s * integrate_I_Qbkg
    plt.figure(4)
    plt.plot(Q, Isample_Q1)
    plt.plot(operative_Q, Isample_Q2)
    plt.plot(integrate_Q, Isample_Q3)
    plt.grid()
    plt.show
    
    
    alpha1 = calc_alpha(J_Q1, Sinf1, Q, Isample_Q1, fe_Q1, Ztot1, rho0)
    alpha2 = calc_alpha(J_Q2, Sinf2, operative_Q, Isample_Q2, fe_Q2, Ztot2, rho0)
    alpha3 = calc_alpha(J_Q3, Sinf3, integrate_Q, Isample_Q3, fe_Q3, Ztot3, rho0)
    print(alpha1, alpha2, alpha3)
    
    Icoh_Q1 = calc_Icoh(N, alpha1, Isample_Q1, Iincoh_Q1)
    Icoh_Q2 = calc_Icoh(N, alpha2, Isample_Q2, Iincoh_Q2)
    Icoh_Q3 = calc_Icoh(N, alpha3, Isample_Q3, Iincoh_Q3)
    plt.figure(5)
    plt.plot(Q, Icoh_Q1)
    plt.plot(operative_Q, Icoh_Q2)
    plt.plot(integrate_Q, Icoh_Q3)
    plt.grid()
    plt.show()
    

    # plt.figure(1)
    # plt.plot(used_Q, J_Q)
    # plt.grid()
    # plt.show

    
    # print("used_Q = ", used_Q.size)
    # print("J_Q = ", J_Q.size)

    # if range_flag == "operative":
        # S_Q = calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, min_index, max_index)
    # elif range_flag == "integrate":
        # S_Q = calc_SQ2()
    # else:
        # S_Q = calc_SQ3(N, Icoh_Q, Ztot, fe_Q)
        
    # # used_Q = np.concatenate([Q[min_index], used_Q])
    # print("S(Q) = ", S_Q.size)
    # print("used_Q = ", used_Q.size)
    
    # # plt.figure(1)
    # # plt.plot(used_Q, S_Q)
    # # plt.grid()
    # # plt.show
    
    # i_Q = calc_iQ(S_Q, Sinf)
    # Qi_Q = used_Q * i_Q
    
    # print("i(Q) = ", i_Q.size)
    # print("Qi(Q) = ", Qi_Q.size)
    
    # DeltaQ = np.diff(Q)
    # meanDeltaQ = np.mean(DeltaQ)
    # r = fftpack.fftfreq(used_Q.size, meanDeltaQ)
    # mask = np.where(r>0)
    
    # F_r3 = calc_Fr(r[mask], used_Q, i_Q)
    
    # # plt.figure(2)
    # # plt.plot(r[mask], F_r3)
    # # plt.grid()
    # # plt.show
    
    # print("F(r) = ", F_r3.size)
    # print("J(Q) = ", J_Q.size)
    
    # iteration = 4
    # r_cutoff = 0.25
    # F_rInt = calc_optimize_Fr(iteration, F_r3, rho0, i_Q, used_Q, Sinf, J_Q, r[mask], r_cutoff)
    
    
    # # plt.figure(2)
    # # plt.plot(r[mask], F_rInt)
    # # plt.show()
