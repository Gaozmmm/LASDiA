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
This script is mainly used for testing the software, but it can be used to run LASDiA in text mode.
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
from modules.UtilityAnalysis import *
from modules.Optimization import *
from modules.Minimization import *
from modules.Formalism import *
from modules.IgorFunctions import *

# import cmath
# from cmath import exp, pi

if __name__ == '__main__':
    N = 1 # sc.N_A
    smooth_factor = 0.25
    damp_factor = 1
    iteration = 2
    rmin = 0.22
    # xyz_file = "./xyzFiles/co2.xyz"
    xyz_file = "./xyzFiles/argon.xyz"

    # Q, I_Q = read_file("./HT2_034T++.chi")
    # Qbkg, I_Qbkg = read_file("./HT2_036T++.chi")
    Q, I_Q = read_file("../data/cea_files/Ar/HT2_034T++_rem.chi")
    Qbkg, I_Qbkg = read_file("../data/cea_files/Ar/HT2_036T++_rem.chi")
    # Q, I_Q = read_file("../data/cea_files/CO2/WO2_007BBin.chi")
    # Qbkg, I_Qbkg = read_file("../data/cea_files/CO2/WO2_013BBin.chi")
    # Q, I_Q = read_file("../data/cea_files/CO2/WO2_007T++.chi")
    # Qbkg, I_Qbkg = read_file("../data/cea_files/CO2/WO2_013T++.chi")


    # plt.figure('data')
    # plt.plot(Q, I_Q, label='F(r) Corr')
    # plt.xlabel('r ($nm$)')
    # plt.ylabel('F(r)')
    # plt.legend()
    # plt.grid()
    # plt.show

    # plt.figure('data')
    # plt.plot(Qbkg, I_Qbkg, label='F(r) Corr')
    # plt.xlabel('r ($nm$)')
    # plt.ylabel('F(r)')
    # plt.legend()
    # # plt.grid()
    # plt.show()


    # Ar
    minQ = 3
    maxQ = 109
    # CO2
    # minQ = 8.005
    # maxQ = 100
    # QmaxIntegrate = 90
    QmaxIntegrate = np.arange(60, 100, 2.5)

    # min_index, max_index = calc_indices(Q, minQ, QmaxIntegrate, maxQ)
    # validation_index, integration_index, calculation_index = calc_ranges(Q, minQ, QmaxIntegrate, maxQ)

    elementList = {"Ar":1}
    # elementList = {"C":1,"O":2}

    # test values
    # Ar
    s = np.arange(0.2, 1.0, 0.1)
    # CO2
    # s = np.arange(0.60, 1.05, 0.01)
    rho0 = np.arange(20, 31, 1)

    # # real values
    # # s = np.arange(0.2, 0.8, 0.01)
    # # rho0 = np.arange(24.0, 28.0, 0.1)

    # best values
    # Ar
    # s = np.array([0.57])
    # rho0 = np.array([26.1])
    # CO2
    # s = np.array([0.984228])
    # rho0 = np.array([29.6625])
    chi2 = np.zeros((rho0.size, s.size))

    # remember the electron unit in atomic form factor!!!
    fe_Q, Ztot = calc_eeff(elementList, Q)
    Iincoh_Q = calc_Iincoh(elementList, Q)
    J_Q = calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf = calc_Sinf(elementList, fe_Q, Q, Ztot)

    # QIsample, Isample_QIgor = read_file("../data/cea_files/CO2/WO2_007Subt.chi")
    
    best_rho0_s = np.zeros((len(QmaxIntegrate), 2))
    
    for l, val_QmaxIntegrate in enumerate(QmaxIntegrate):
        min_index, max_index = calc_indices(Q, minQ, QmaxIntegrate[l], maxQ)
        validation_index, integration_index, calculation_index = calc_ranges(Q, minQ, QmaxIntegrate[l], maxQ)

        for i, val_rho0 in enumerate(rho0):
            for j, val_s in enumerate(s):
                Isample_Q = calc_IsampleQ(I_Q, s[j], I_Qbkg)
                alpha = calc_alpha(J_Q[integration_index], Sinf, Q[integration_index], Isample_Q[integration_index], fe_Q[integration_index], Ztot, rho0[i])
                Icoh_Q = calc_Icoh(N, alpha, Isample_Q, Iincoh_Q)
                # alpha = calc_alpha(J_Q[integration_index], Sinf, Q[integration_index], Isample_QIgor[integration_index], fe_Q[integration_index], Ztot, rho0[i])
                # Icoh_Q = calc_Icoh(N, alpha, Isample_QIgor, Iincoh_Q)

                S_Q = calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, max_index, integration_index)
                newQ, S_Qsmoothed = calc_SQsmoothing(Q[validation_index], S_Q[validation_index], Sinf, smooth_factor, min_index, minQ, QmaxIntegrate[l], maxQ, 550)
                S_QsmoothedDamp = calc_SQdamp(S_Qsmoothed, newQ, Sinf, QmaxIntegrate[l], damp_factor)

                Qi_Q = calc_QiQ(newQ, S_QsmoothedDamp, Sinf)
                i_Q = calc_iQ(S_QsmoothedDamp, Sinf)

                validation_indexSmooth, integration_indexSmooth, calculation_indexSmooth = calc_ranges(newQ, minQ, QmaxIntegrate[l], maxQ)
                min_indexSmooth, max_indexSmooth = calc_indices(newQ, minQ, QmaxIntegrate[l], maxQ)

                r = calc_r(newQ)
                F_r = calc_Fr(r, newQ[integration_indexSmooth], Qi_Q[integration_indexSmooth])

                numAtoms, element, x, y, z = read_xyz_file(xyz_file)
                iintra_Q, fe_QSmooth = calc_iintra(newQ, max_indexSmooth, elementList, element, x, y, z)
                Qiintradamp = calc_iintradamp(iintra_Q, newQ, QmaxIntegrate[l], damp_factor)
                Fintra_r = calc_Fr(r, newQ[integration_indexSmooth], Qiintradamp[integration_indexSmooth])
                # Fintra_r = calc_Fintra(r, newQ, QmaxIntegrate[l])

                Iincoh_QSmooth = calc_Iincoh(elementList, newQ)
                J_QSmooth = calc_JQ(Iincoh_QSmooth, Ztot, fe_QSmooth)

                F_rIt = calc_optimize_Fr(iteration, F_r, Fintra_r, rho0[i], i_Q[integration_indexSmooth], newQ[integration_indexSmooth], Sinf, J_QSmooth[integration_indexSmooth], r, rmin)
                # write_file("../Results/CO2/F_rCorr.txt", r, F_rIt)

                maskIt = np.where((r>0) & (r < rmin))
                rIt = r[maskIt]
                deltaF_r = calc_deltaFr(F_rIt[maskIt], Fintra_r[maskIt], rIt, rho0[i])

                chi2[i][j] = simps(deltaF_r**2, r[maskIt])
                # print(chi2[i][j])

        # take min of chi2
        minIndxRho0, minIndxS = np.unravel_index(chi2.argmin(), chi2.shape)
        # print(chi2[minIndxRho0][minIndxS])
        # print(rho0[minIndxRho0], s[minIndxS])
        
        best_rho0_s[l][0] = rho0[minIndxRho0]
        best_rho0_s[l][1] = s[minIndxS]
        
        print(best_rho0_s[l][0], best_rho0_s[l][1])

    plt.figure('QmaxIntegrate rho0')
    plt.plot(QmaxIntegrate, best_rho0_s[ : ,0])
    plt.axis([59, 101, 23.9, 25.1])
    plt.xlabel('QmaxIntegrate($nm^{-1}$)')
    plt.ylabel('rho0')
    plt.grid()
    plt.show
    
    plt.figure('QmaxIntegrate s')
    plt.plot(QmaxIntegrate, best_rho0_s[ : ,1])
    plt.axis([59, 101, 0.59, 0.81])
    plt.xlabel('QmaxIntegrate($nm^{-1}$)')
    plt.ylabel('s')
    plt.grid()
    plt.show
    
    best_rho0 = np.mean(best_rho0_s[ : ,0])
    
    plt.ion()
    for i in range(len(QmaxIntegrate)):
        min_index, max_index = calc_indices(Q, minQ, QmaxIntegrate[i], maxQ)
        validation_index, integration_index, calculation_index = calc_ranges(Q, minQ, QmaxIntegrate[i], maxQ)
        
        Isample_Q = calc_IsampleQ(I_Q, best_rho0_s[i,1], I_Qbkg)
        alpha = calc_alpha(J_Q[integration_index], Sinf, Q[integration_index], Isample_Q[integration_index], fe_Q[integration_index], Ztot, best_rho0)
        Icoh_Q = calc_Icoh(N, alpha, Isample_Q, Iincoh_Q)
        S_Q = calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, max_index, integration_index)
        
        plt.figure('S_Q QmaxIntegrate')
        plt.plot(Q[validation_index], S_Q, label='%s' %QmaxIntegrate[i])
        plt.legend()
        plt.draw()
        
        time.sleep(1.0)
        
    plt.grid()
    plt.ioff()
    plt.show()

    
    # plot 2d chi2
    # plt.figure('chi2s')
    # plt.plot(s,chi2[0, : ])
    # plt.xlabel('s')
    # plt.ylabel('chi2')
    # plt.grid()
    # plt.show
    
    # plt.figure('chi2rho0')
    # plt.plot(rho0,chi2[ : ,0])
    # plt.xlabel('rho0')
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
