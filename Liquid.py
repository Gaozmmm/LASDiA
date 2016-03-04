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
from modules.Formalism import *
from modules.IgorFunctions import *

# import cmath
# from cmath import exp, pi

if __name__ == '__main__':
    N = 3 # sc.N_A
    
    # Q, I_Q = read_file("../data/cea_files/Ar/HT2_034T++.chi")
    # Qbkg, I_Qbkg = read_file("../data/cea_files/Ar/HT2_036T++.chi")
    # Q, I_Q = read_file("../data/cea_files/Ar/HT2_034T++_rem.chi")
    # Qbkg, I_Qbkg = read_file("../data/cea_files/Ar/HT2_036T++_rem.chi")
    Q, I_Q = read_file("../data/cea_files/CO2/WO2_007BBin.chi")
    Qbkg, I_Qbkg = read_file("../data/cea_files/CO2/WO2_013BBin.chi")
    # Q, I_Q = read_file("../data/cea_files/CO2/WO2_007T++.chi")
    # Qbkg, I_Qbkg = read_file("../data/cea_files/CO2/WO2_013T++.chi")
    
    # Ar
    # minQ = 3
    # maxQ = 109
    # CO2
    minQ = 8.005
    maxQ = 100
    QmaxIntegrate = 90
    
    min_index = np.where(Q<=minQ)
    max_index = np.where((Q>QmaxIntegrate) & (Q<=maxQ))
    validation_index = np.where(Q<=maxQ)
    integration_index = np.where(Q<=QmaxIntegrate)
    calculation_index = np.where((Q>minQ) & (Q<=QmaxIntegrate))
    
    # elementList = {"Ar":1}
    elementList = {"C":1,"O":2}
    
    # test values
    # Ar
    # s = np.arange(0.2, 1.0, 0.1)
    # CO2
    # s = np.arange(0.95, 1.05, 0.01)
    # rho0 = np.arange(26, 31, 1)
    
    # # real values
    # # s = np.arange(0.2, 0.8, 0.01)
    # # rho0 = np.arange(24.0, 28.0, 0.1)
    
    # best values
    # Ar
    # s = np.array([0.57])
    # rho0 = np.array([26.1])
    # CO2
    # s = np.array([1.04346])
    # s = np.array([0.79])
    # rho0 = np.array([29.7404])
    
    # s = np.array([1.00114])
    s = np.array([1.0])
    rho0 = np.array([29.7877])
    
    chi2 = np.zeros((rho0.size, s.size))
    
    # remember the electron unit in atomic form factor!!!
    fe_Q, Ztot = calc_eeff(elementList, Q)
    Iincoh_Q = calc_Iincoh(elementList, Q)
    J_Q = calc_JQ(Iincoh_Q, Ztot, fe_Q)
    # J_QIgor = calc_JQIgor(Iincoh_Q, fe_Q)
    Sinf = calc_Sinf(elementList, fe_Q, Q, Ztot)
    
    # write_file("../Results/CO2/J_Q.txt", Q, J_Q)
    # write_file("../Results/CO2/J_Q_formIgor.txt", Q, J_QIgor)
    
    # plt.figure('J_Q')
    # plt.plot(Q, J_Q, label='J(Q)')
    # plt.plot(Q, J_QIgor, label='J(Q) Igor')
    # plt.xlabel('Q ($nm^{-1}$)')
    # plt.ylabel('J(Q)')
    # plt.legend()
    # plt.grid()
    # plt.show()
    
    QIsample, Isample_QIgor = read_file("../data/cea_files/CO2/WO2_007Subt.chi")
    
    for i, val_rho0 in enumerate(rho0):
        for j, val_s in enumerate(s):
            # Isample_Q = calc_IsampleQ(I_Q, s[j], I_Qbkg)
            alpha = calc_alpha(J_Q, Sinf, Q, Isample_QIgor, fe_Q, Ztot, rho0[i], integration_index)
            Icoh_Q = calc_Icoh(N, alpha, Isample_QIgor, Iincoh_Q)
            
            # plt.figure('Isample_Q')
            # plt.plot(Q, Isample_Q, label='Isample(Q)')
            # plt.xlabel('Q ($nm^{-1}$)')
            # plt.ylabel('Isample(Q)')
            # plt.legend()
            # plt.grid()
            # plt.show()
            
            S_Q = calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, min_index, max_index, calculation_index)
            # S_Qdamp = calc_SQdamp(S_Q, Q[validation_index], Sinf, QmaxIntegrate, 1)
            S_Qsmooth = SQsmoothing(Q, S_Q, Sinf, 0.25, min_index, max_index, validation_index)
            S_Qsmooth2 = SQsmoothing2(Q, S_Q, Sinf, 0.25, min_index, max_index, calculation_index)
            newQ, S_Qsmooth3 = SQsmoothing3(Q, S_Q, Sinf, 0.25, minQ, maxQ, QmaxIntegrate)
            
            # write_file("../Results/CO2/S_Q.txt", Q[validation_index], S_Q)
            
            plt.figure('S_Q')
            plt.plot(Q[validation_index], S_Q, label='S(Q)')
            # plt.plot(Q[validation_index], S_Qdamp, label='S(Q) damped')
            plt.plot(Q[validation_index], S_Qsmooth,label='S(Q) smoothed')
            plt.plot(Q[validation_index], S_Qsmooth2,label='S(Q) smoothed2')
            plt.plot(newQ, S_Qsmooth3,label='S(Q) smoothed3')
            plt.xlabel('Q ($nm^{-1}$)')
            plt.ylabel('S(Q)')
            plt.legend()
            plt.grid()
            plt.show()
            
            # i_Q = calc_iQ(S_Qsmooth, Sinf)
            # Qi_Q = Q[validation_index] * i_Q
            
            # # plt.figure('i_Q')
            # # plt.plot(Q[validation_index], i_Q, label='i(Q)')
            # # plt.xlabel('Q ($nm^{-1}$)')
            # # plt.ylabel('i(Q)')
            # # plt.legend()
            # # plt.grid()
            # # plt.show
            
            # # plt.figure('Qi_Q')
            # # plt.plot(Q[validation_index], Qi_Q, label='Qi(Q)')
            # # plt.xlabel('Q ($nm^{-1}$)')
            # # plt.ylabel('Qi(Q)')
            # # plt.legend()
            # # plt.grid()
            # # plt.show()
                        
            # DeltaQ = np.diff(Q)
            # meanDeltaQ = np.mean(DeltaQ)
            # r = fftpack.fftfreq(Q[validation_index].size, meanDeltaQ)
            # mask = np.where(r>0)
            
            # F_r, F_r2 = calc_Fr(r[mask], Q[integration_index], i_Q[integration_index])
            
            # # write_file("../Results/CO2/F_r.txt", r[mask], F_r)
            
            # plt.figure('F_r')
            # plt.plot(r[mask], F_r, label='F(r)')
            # plt.plot(r[mask], F_r, label='F2(r)')
            # plt.xlabel('r ($nm$)')
            # plt.ylabel('F(r)')
            # plt.legend()
            # plt.grid()
            # plt.show()
            
            # iteration = 2
            # rmin = 0.22
            
            # Fintra_r = calc_Fintra(r[mask], Q, QmaxIntegrate)
            # F_rIt = calc_optimize_Fr(iteration, F_r, Fintra_r, rho0[i], i_Q[integration_index], Q[integration_index], Sinf, J_Q[integration_index], r[mask], rmin)
            
            # maskIt = np.where((r>0) & (r < rmin))
            # rIt = r[maskIt]
            # Fintra_r = calc_Fintra(rIt, Q, QmaxIntegrate)
            # deltaF_r = calc_deltaFr(F_rIt[maskIt], Fintra_r, rIt, rho0[i])
            # chi2[i][j] = simps(deltaF_r**2, r[maskIt])
            # print(chi2[i][j])
    
    # take min of chi2
    # minIndxRho0, minIndxS = np.unravel_index(chi2.argmin(), chi2.shape)
    # print(chi2[minIndxRho0][minIndxS])
    # print(rho0[minIndxRho0], s[minIndxS])
    
    # # plot 2d chi2
    # plt.figure('chi2')
    # plt.plot(s,chi2[0, : ])
    # plt.xlabel('s')
    # # plt.plot(rho0,chi2[ : ,0])
    # # plt.xlabel('rho0')
    # plt.ylabel('chi2')
    # plt.grid()
    # plt.show()
    
    # plot the 3d chi2 and its profile
    # x, y = np.meshgrid(s, rho0)
    # fig = plt.figure('chi2 3D')
    # ax = Axes3D(fig)
    # ax.set_xlabel('s')
    # ax.set_ylabel(r'$\rho_0$ (atoms/$nm^3$)')
    # ax.set_zlabel(r'$\chi^2$')
    # ax.plot_surface(x, y, chi2, rstride=1, cstride=1, cmap='rainbow')
    
    # plt.figure('chi2 Profile')
    # plt.contour(s, rho0, chi2, 200)
    # plt.xlabel('s')
    # plt.ylabel(r'$\rho_0$ (atoms/$nm^3$)')
    # plt.show()
