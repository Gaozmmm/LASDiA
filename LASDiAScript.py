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

"""LASDiA main script file.
This script is mainly used for testing the software.
The nomenclature and the procedures follow the article: Eggert et al. 2002 PRB, 65, 174105
"""

from __future__ import (absolute_import, division, print_function, unicode_literals)
import six

import sys
import os

import scipy.constants as sc
from scipy import fftpack
from scipy import signal
from scipy.integrate import simps
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np
import time
import math

from modules.MainFunctions import *
from modules.Utility import *
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
    
    min_index, max_index = calc_indices(Q, minQ, QmaxIntegrate, maxQ)
    validation_index, integration_index, calculation_index = calc_ranges(Q, minQ, QmaxIntegrate, maxQ)
    
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
    Sinf = calc_Sinf(elementList, fe_Q, Q, Ztot)
    
    QIsample, Isample_QIgor = read_file("../data/cea_files/CO2/WO2_007Subt.chi")
    
    for i, val_rho0 in enumerate(rho0):
        for j, val_s in enumerate(s):
            # Isample_Q = calc_IsampleQ(I_Q, s[j], I_Qbkg)
            alpha = calc_alpha(J_Q, Sinf, Q, Isample_QIgor, fe_Q, Ztot, rho0[i], integration_index)
            Icoh_Q = calc_Icoh(N, alpha, Isample_QIgor, Iincoh_Q)
            
            S_Q = calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, max_index, integration_index)
            newQ, S_Qsmoothed = calc_SQsmoothing(Q[validation_index], S_Q[validation_index], Sinf, 0.25, min_index, minQ, QmaxIntegrate, maxQ, 550)
            
            S_QsmoothedDamp = calc_SQdamp(S_Qsmoothed, newQ, Sinf, QmaxIntegrate, 1)
            
            i_Q = calc_iQ(S_QsmoothedDamp, Sinf)
            Qi_Q = calc_QiQ(newQ, S_QsmoothedDamp, Sinf)
            
            print(Qi_Q[newQ>QmaxIntegrate])
            
            newDim = 2*2*2**math.ceil(math.log(5*(QmaxIntegrate+1))/math.log(2))
            oldDim = Qi_Q.size
            addDim = newDim - oldDim
            Qi_Qzero = np.zeros(addDim)
            Qi_Q2 = np.concatenate([Qi_Q, Qi_Qzero])
            # Qi_Q2 = np.resize(Qi_Q, newDim)
            newQ2 = np.linspace(np.amin(Q), maxQ, newDim, endpoint=True)
            # Qi_Q2[newQ>QmaxIntegrate] = 0.0
            
            # write_file("../Results/CO2/Qi_Q.txt", newQ, Qi_Q)
            
            # plt.figure('Qi_Q')
            # # plt.plot(newQ, Qi_Q, label='Qi(Q)')
            # plt.plot(newQ2, Qi_Q2, label='Qi(Q) resized')
            # plt.xlabel('Q ($nm^{-1}$)')
            # plt.ylabel('Qi(Q)')
            # plt.legend()
            # plt.grid()
            # plt.show
            
            # plt.figure('Qi_Q resized')
            # # plt.plot(newQ, Qi_Q, label='Qi(Q)')
            # plt.plot(Qi_Q2, label='Qi(Q) resized')
            # plt.xlabel('Q ($nm^{-1}$)')
            # plt.ylabel('Qi(Q)')
            # plt.legend()
            # plt.grid()
            # plt.show()
                        
            DeltaQ = np.diff(newQ2)
            meanDeltaQ = np.mean(DeltaQ)
            r = fftpack.fftfreq(newQ2.size, meanDeltaQ)
            # r2 = np.fft.fftfreq(newQ.size, meanDeltaQ)
            mask = np.where(r>0)
            
            F_r, F_r2 = calc_Fr(r[mask], newQ, i_Q, Qi_Q2)
            
            # # write_file("../Results/CO2/F_r.txt", r[mask], F_r)
            
            plt.figure('F_r')
            plt.plot(r[mask], F_r, label='F(r)')
            # plt.plot(r[mask], F_r2[mask], label='F(r)')
            plt.xlabel('r ($nm$)')
            plt.ylabel('F(r)')
            plt.legend()
            plt.grid()
            plt.show()
            
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
