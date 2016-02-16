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

"""Environment for testing the software
The nomenclature and the procedures follow the article: Eggert et al. 2002 PRB, 65, 174105
"""

from __future__ import (absolute_import, division, print_function, unicode_literals)
import six

import sys
import os

import scipy.constants as sc
from scipy import fftpack
from scipy import signal
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np
import time
import math

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
    # Q, I_Q = read_file("./data/cea_files/WO2_007BBin.chi")
    # Qbkg, I_Qbkg = read_file("./data/cea_files/WO2_013BBin.chi")
    
    BinNum = 1
    Num = 2048 #len(Q)
    maxQ = 109
    minQ = 3
    rebQ, rebI_Q = rebinning(Q, I_Q, BinNum, Num, maxQ, minQ)
    
    plt.figure('RawData')
    plt.plot(Q, I_Q)
    plt.plot(rebQ, rebI_Q)
    plt.xlabel('Q')
    plt.ylabel('I(Q)')
    plt.grid()
    plt.show()
    
    # newI_Q = removePeaks(Q, I_Q)
    # newI_Qbkg = removePeaks(Qbkg, I_Qbkg)
    
    # plt.figure('RawData Removed')
    # plt.plot(Q, newI_Q)
    # plt.plot(Qbkg, newI_Qbkg)
    # plt.xlabel('Q')
    # plt.ylabel('I(Q) Removed')
    # plt.grid()
    # plt.show()
    
    # minQ = 4
    # maxQ = 109
    # # minQ = 8
    # # maxQ = 100
    # QmaxIntegrate = 90
    
    # BinNum = 1
    # Num = 2048 # len(Q)
    
    # min_index = np.where(Q<=minQ)
    # max_index = np.where((Q>QmaxIntegrate) & (Q<=maxQ))
    # validation_index = np.where(Q<=maxQ)
    # integration_index = np.where(Q<=QmaxIntegrate)
    # calculation_index = np.where((Q>minQ) & (Q<=QmaxIntegrate))
    
    # rebinQ = rebinning(Q, BinNum, Num, maxQ)
    # interpI_Q = interpolation(Q, I_Q, rebinQ)
    # interpI_Q[min_index] = 0.0
    
    
    
    

    # elementList = {"Ar":1}
    # # elementList = {"C":1,"O":2}
    # # test values
    # # s = np.arange(0.79, 0.83, 0.01)
    # # s = np.arange(0.2, 1.0, 0.1)
    # # rho0 = np.arange(24, 28, 1)
    
    # # real values
    # # s = np.arange(0.2, 1.0, 0.01)
    # # rho0 = np.arange(24.0, 28.0, 0.1)
    
    # # best values
    # # s = np.array([0.81])
    # # rho0 = np.array([26.1])
    # s = np.array([0.57])
    # rho0 = np.array([26.1])
    
    # chi2 = np.zeros((rho0.size, s.size))
    
    # # remember the electron unit in atomic form factor!!!
    # fe_Q, Ztot = calc_eeff(elementList, Q)
    # Iincoh_Q = calc_Iincoh(elementList, Q)
    # J_Q = calc_JQ(Iincoh_Q, Ztot, fe_Q)
    # Sinf = calc_Sinf(elementList, fe_Q, Q, Ztot)
    
    # for i, val_rho0 in enumerate(rho0):
        # for j, val_s in enumerate(s):
            # Isample_Q = calc_IsampleQ(I_Q, s[j], I_Qbkg)
            # alpha = calc_alpha(J_Q, Sinf, Q, Isample_Q, fe_Q, Ztot, rho0[i], integration_index)
            # Icoh_Q = calc_Icoh(N, alpha, Isample_Q, Iincoh_Q)
            
            # S_Q = calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, min_index, max_index, calculation_index, QmaxIntegrate)
            
            # # plt.figure('S_Q')
            # # plt.plot(Q[validation_index], S_Q)
            # # plt.xlabel('Q')
            # # plt.ylabel('S(Q)')
            # # plt.grid()
            # # plt.show
            
            # # damp = calc_damp(Q[validation_index], QmaxIntegrate)
            # # S_Qdamp = (damp * (S_Q - Sinf)) + Sinf
            
            # # easy test!!!
            # # S_Q = S_Qdamp
            
            # i_Q = calc_iQ(S_Q, Sinf)
            # Qi_Q = Q[validation_index] * i_Q
            
            # # plt.figure('Qi_Q')
            # # plt.plot(Q[validation_index], Qi_Q)
            # # plt.xlabel('Q')
            # # plt.ylabel('i(Q)')
            # # plt.grid()
            # # plt.show
            
            # DeltaQ = np.diff(Q)
            # meanDeltaQ = np.mean(DeltaQ)
            # r = fftpack.fftfreq(Q[validation_index].size, meanDeltaQ)
            # mask = np.where(r>0)
            
            # F_r = calc_Fr(r[mask], Q[integration_index], i_Q[integration_index])
            
            # # plt.figure('F_r')
            # # plt.plot(r[mask], F_r)
            # # plt.xlabel('r')
            # # plt.ylabel('F(r)')
            # # plt.grid()
            # # plt.show
            
            # iteration = 2
            # rmin = 0.24
            # # rmin = 0.22
            
            # Fintra_r = calc_Fintra(r[mask], Q, QmaxIntegrate)
            # # F_rIt = calc_optimize_Fr(iteration, F_r, Fintra_r, rho0[i], i_Q[validation_index], Q[validation_index], Sinf, J_Q[validation_index], r[mask], rmin)
            # F_rIt = calc_optimize_Fr(iteration, F_r, Fintra_r, rho0[i], i_Q[integration_index], Q[integration_index], Sinf, J_Q[integration_index], r[mask], rmin)
            
            # S_Qcorr = calc_SQ_F(F_rIt, r[mask], Q, Sinf)
            
            # plt.figure('S_Qcorr')
            # plt.plot(Q, S_Qcorr)
            # plt.xlabel('Q')
            # plt.ylabel('S(Q)corr')
            # plt.grid()
            # plt.show()
            
            # deltaF_r = calc_deltaFr(F_r, Fintra_r, r[mask], rho0[i])
            # i_Q1 = calc_iQi(i_Q[validation_index], Q[validation_index], Sinf, J_Q[validation_index], deltaF_r, r[mask], rmin)
            # print(i_Q1)
            # F_r1 = calc_Fr(r[mask], Q[integration_index], i_Q1[integration_index])
            
            # deltaF_r2 = calc_deltaFr(F_r1, Fintra_r, r[mask], rho0[i])
            # i_Q2 = calc_iQi(i_Q1[validation_index], Q[validation_index], Sinf, J_Q[validation_index], deltaF_r2, r[mask], rmin)
            # F_r2 = calc_Fr(r[mask], Q[integration_index], i_Q2[integration_index])
            
            
            
            # plt.figure(1)
            # plt.plot(Q[validation_index], i_Qi)
            # plt.plot(Q[validation_index], i_Qi1)
            # plt.plot(Q[validation_index], i_Qi - i_Qi1)
            # plt.xlabel('Q')
            # plt.ylabel('i(Q)i')
            # plt.grid()
            # plt.show()

            
            
            
            # i_Q = i_Q[validation_index] - ( 1/Q[validation_index] * ( i_Q[validation_index] / (Sinf + J_Q[validation_index]) + 1))
            # print(i_Qi)
            
            # F_r = (2.0 / np.pi) * simps(Q * i_Q * np.array(np.sin(np.mat(Q).T * np.mat(r))).T, Q)
            
            
            
            # print(Q[integration_index] * i_Q1[integration_index])
            # print(i_Q1[0])
            # print(Q[integration_index] * i_Q1[integration_index] * np.array(np.sin(np.mat(Q[integration_index]).T * np.mat(r[mask]))).T)
            # print(simps(Q[integration_index] * i_Q1[integration_index] * np.array(np.sin(np.mat(Q[integration_index]).T * np.mat(r[mask]))).T,Q[integration_index]))
            # print((2.0 / np.pi)*simps(Q[integration_index] * i_Q1[integration_index] * np.array(np.sin(np.mat(Q[integration_index]).T * np.mat(r[mask]))).T,Q[integration_index]))
            
            # plt.figure(1)
            # plt.plot(Q[integration_index], Q[integration_index]*i_Q1[integration_index])
            # # plt.plot(r[mask], F_r1)
            # # plt.plot(r[mask], F_r2)
            # plt.xlabel('Q')
            # plt.ylabel('Qi(Q)')
            # plt.grid()
            # plt.show()
            












            
           
            
            # maskIt = np.where((r>0) & (r < rmin))
            # rIt = r[maskIt]
            # Fintra_r = calc_Fintra(rIt)
            # deltaF_r = calc_deltaFr(F_rIt[maskIt], Fintra_r, rIt, rho0[i])
            # chi2[i][j] = simps(deltaF_r**2, rIt)
        
    # minIndxRho0, minIndxS = np.unravel_index(chi2.argmin(), chi2.shape)
    # print(chi2[minIndxRho0][minIndxS])
    # print(rho0[minIndxRho0], s[minIndxS])
    
    # print(chi2)
    # print(chi2[0, : ])
    # print(chi2.shape)
    
    # plt.figure(5)
    # plt.plot(s,chi2[0, : ])
    # # plt.plot(rho0,chi2[ : ,0])
    # # plt.xlabel('rho0')
    # plt.xlabel('s')
    # plt.ylabel('chi2')
    # plt.grid()
    # plt.show()
    
    # x, y = np.meshgrid(s, rho0)
    # fig = plt.figure(3)
    # ax = Axes3D(fig)
    # ax.plot_surface(x, y, chi2, rstride=1, cstride=1, cmap='rainbow')
    
    # plt.figure(4)
    # plt.contour(s, rho0, chi2, 200)
    # plt.show()
