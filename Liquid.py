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
import matplotlib.pyplot as plt
import numpy as np

from modules.Utility import *
from modules.LiquidStructure import *
from modules.InterpolateData import *
from modules.Optimization import *
from modules.Minimization import *

if __name__ == '__main__':

    # Q, I_Q = read_file("./data/cea_files/HT2_034.chi")
    # Qbkg, I_Qbkg = read_file("./data/cea_files/HT2_036.chi")

    Q, I_Q = read_file("./data/my_test/HT2_034_20151104.chi")
    Qbkg, I_Qbkg = read_file("./data/my_test/HT2_036_20151104.chi")
    
    elementList = {"Ar":1}
    
    fe_Q, Ztot = calc_eeff(elementList, Q)
    
    numAtoms = sc.N_A
    alpha = 1
    s = 1
    Isample_Q = I_Q - s * I_Qbkg
    Iincoh_Q = calc_Iincoh(elementList, Q)
    Icoh_Q = calc_Icoh(numAtoms, alpha, Isample_Q, Iincoh_Q)
    
    S_Q = calc_SQ(Icoh_Q, Ztot, fe_Q)

    plt.figure(1)
    plt.plot(Q, S_Q)
    plt.grid()
    plt.show
    
    Sinf = calc_Sinf(elementList, fe_Q, Q, Ztot)
    i_Q = calc_iQ(S_Q, Sinf)
    r = calc_spectrum(i_Q)
    F_r = calc_Fr(r, Q, i_Q)
    
    plt.figure(2)
    plt.plot(r, F_r)
    plt.grid()
    plt.show
    
    rho0 = calc_rho0(Q, i_Q)
    g_r = calc_gr(r, F_r, rho0)
    
    # plt.figure(3)
    # plt.plot(r, g_r)
    # plt.grid()
    # plt.show
    
    iteration = 2
    r_cutoff = 0.25
    Fintra_r = calc_Fintra()
    J_Q = calc_JQ(Iincoh_Q, Ztot, fe_Q)
    optF_r = calc_optimize_Fr(iteration, F_r, Fintra_r, rho0, i_Q, Q, Sinf, J_Q, r, r_cutoff)
    
    plt.figure(2)
    plt.plot(r, optF_r)
    #plt.grid()
    plt.show
    
    SQ1 = calc_SQi(optF_r, r, Q, Sinf)
    
    plt.figure(1)
    plt.plot(Q, SQ1)
    #plt.grid()
    plt.show()
    

    
    # #-------------------------------------------------------------
    # Very important for the FFT!!!
    # r = fftpack.fftfreq(S_Q.size, Q[1] - Q[0])    
    # F_r = fftpack.fft(S_Q)
    # F_rI = F_r.imag

    # pidxs = np.where(r > 0)
    # r = r[pidxs]
    # F_rI = F_rI[pidxs]
    
    # print(r.size, F_rI.size)
    
    # plt.figure(2)
    # plt.plot(r, F_rI)
    # plt.grid()
    # plt.show()
    
    # #-------------------------------------------------------------    

    # QiQ = Q*(S_Q - Sinf)
    
    # r, F_r = FFT_QiQ(QiQ)
