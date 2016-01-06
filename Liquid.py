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
    N = 1#sc.N_A
    # Q, I_Q = read_file("./data/cea_files/HT2_034.chi")
    # Qbkg, I_Qbkg = read_file("./data/cea_files/HT2_036.chi")

    # Q, I_Q = read_file("./data/my_test/HT2_034_20151104.chi")
    # Qbkg, I_Qbkg = read_file("./data/my_test/HT2_036_20151104.chi")

    Q, I_Q = read_file("./data/cea_files/HT2_034T++.chi")
    Qbkg, I_Qbkg = read_file("./data/cea_files/HT2_036T++.chi")

    plt.figure(1)
    plt.plot(Q,I_Q)
    plt.grid()
    plt.show

    plt.figure(1)
    plt.plot(Qbkg,I_Qbkg)
    plt.grid()
    plt.show

    print(np.amin(Q))
    print(np.amax(Q))

    minQ = 3
    maxQ = 109
    QmaxIntegrate = 90

    right_index = np.where((Q>minQ) & (Q<QmaxIntegrate))
    inf_index = np.where((Q>QmaxIntegrate) & (Q<maxQ))
    newQ = Q[right_index]
    newI_Q = I_Q[right_index]
    newI_Qbkg = I_Qbkg[right_index]
    
    elementList = {"Ar":1}

    # remember the electronic unit in atomic form factor!!!
    fe_Q, Ztot = calc_eeff(elementList, newQ)
    
    s =  0.8
    rho0 = 25.06#24.00
    
    Iincoh_Q = calc_Iincoh(elementList, newQ)
    J_Q = calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf = calc_Sinf(elementList, fe_Q, newQ, Ztot)
    
    Isample_Q = newI_Q - s * newI_Qbkg

    alpha = calc_alpha(J_Q, Sinf, newQ, Isample_Q, fe_Q, Ztot, rho0)
    
    Icoh_Q = calc_Icoh(N, alpha, Isample_Q, Iincoh_Q)

    S_Q = calc_SQ(N, Icoh_Q, Ztot, fe_Q)

    
    Qinf = Q[inf_index]
    newQinf = np.concatenate([newQ, Qinf])

    SQinf = np.zeros(Qinf.size)
    SQinf.fill(Sinf)
    newSQinf = np.concatenate([S_Q, SQinf])
    
    plt.figure(2)
    plt.plot(newQ,S_Q)
    plt.grid()
    plt.show

    plt.figure(3)
    plt.plot(newQinf,newSQinf)
    plt.grid()
    plt.show

    i_Q = calc_iQ(newSQinf, Sinf)
    Qi_Q = newQinf * i_Q
    #r = calc_spectrum(i_Q)
    r = np.linspace(0.0, 1.5, newQinf.size)
    F_r = calc_Fr(r, newQinf, i_Q)
    # print(r)
    # r_index = np.where(r<1.6)
    # newr = r[r_index]
    # newF_r = F_r[r_index]

    
    # plt.figure(5)
    # plt.plot(newQinf,i_Q)
    # plt.grid()
    # plt.show

    # plt.figure(6)
    # plt.plot(newQinf,Qi_Q)
    # plt.grid()
    # plt.show

    
    
    plt.figure(4)
    plt.plot(r,F_r)
    plt.grid()
    plt.show()

    

    # qmax = np.amax(Q)
    # Qa = np.arange(qmax,20,Q[1]-Q[0])
    # newQ = np.concatenate([Q, Qa])
    # print(newQ.size)
    # SQa = np.zeros(Qa.size)
    # SQa.fill(1)
    # newSQ = np.concatenate([S_Q, SQa])
    # print(newSQ.size)

    # Iincoh_Q = calc_Iincoh(elementList, newQ)
    # fe_Q, Ztot = calc_eeff(elementList, newQ)
    # J_Q = calc_JQ(Iincoh_Q, Ztot, fe_Q)

    


    # i_Q = calc_iQ(S_Q, Sinf)
    # r = calc_spectrum(i_Q)
    # F_r = calc_Fr(r, Q, i_Q)
#    print(F_r)

    # plt.figure(2)
    # plt.plot(r,F_r)
    # plt.grid()
    # plt.show

    # Fintra_r = calc_Fintra()
    # iteration=2
    # r_cutoff=2
    # optF_r = calc_optimize_Fr(iteration, F_r, Fintra_r, rho0, i_Q, newQ, Sinf, J_Q, r, r_cutoff)

    # plt.figure(3)
    # plt.plot(r,optF_r)
    # plt.grid()
    # plt.show

    
    # plt.figure(4)
    # plt.plot(Q, I_Q)
    # plt.show

    # plt.figure(4)
    # plt.plot(Qbkg, I_Qbkg)
    # plt.show()

    
    # Integral1 = simps((J_Q + Sinf) * Q**2, Q)
    # print("int1 ", Integral1)
    # Integral2 = (J_Q + Sinf) * Q**2
    # print("int2 ", Integral2.sum(axis=0))

    # rQ = np.outer(r,Q)
    # sinrQ = np.sin(rQ)
    # FR = (2.0 / np.pi) * np.sum(Q * i_Q * sinrQ, axis=1)
    # FR /=1000

    # iteration = 1
    # r_cutoff = 1.0
    # Fintra_r = calc_Fintra()
    
    # optF_r = calc_optimize_Fr(iteration, F_r, Fintra_r, rho0, i_Q, Q, Sinf, J_Q, r, r_cutoff)

    # plt.figure(1)
    # plt.plot(r, F_r)
    # plt.grid()
    # plt.show

    
    # plt.figure(1)
    # plt.plot(r, optF_r)
    # plt.grid()
    # plt.show()
    
    # SQ_F = calc_SQ_F(optF_r, r, Q, Sinf)
    
    # plt.figure(1)
    # plt.plot(Q, SQ_F)
    # #plt.grid()
    # plt.show
    
    
    # R, FR = FFT_QiQ(i_Q + 1)
    
    # plt.figure(2)
    # plt.plot(r, FR)
    # plt.grid()
    # plt.show()

    
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
