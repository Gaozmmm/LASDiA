# The MIT License (MIT)

# Copyright (c) 2015-2016 European Synchrotron Radiation Facility

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
from modules.Geometry import *

# import cmath
# from cmath import exp, pi

if __name__ == '__main__':
    smooth_factor = 0.25
    damp_factor = 1
    iteration = 2
    rmin = 0.22
    # xyz_file = "./xyzFiles/co2.xyz"
    xyz_file = "./xyzFiles/argon.xyz"
    numAtoms, element, x, y, z = read_xyz_file(xyz_file)

    # Q, I_Q = read_file("../data/cea_files/Ar/HT2_034T++.chi")
    # Qbkg, I_Qbkg = read_file("../data/cea_files/Ar/HT2_036T++.chi")

    # Q, I_Q = read_file("../data/cea_files/Ar/HT2_034T++_rem.chi")
    # Qbkg, I_Qbkg = read_file("../data/cea_files/Ar/HT2_036T++_rem.chi")

    Q1, I_Q1 = read_file("../data/cea_files/CO2/WO2_007BBin.chi")
    Qbkg1, I_Qbkg1 = read_file("../data/cea_files/CO2/WO2_013BBin.chi")

    Q, I_Q = read_file("../data/cea_files/CO2/WO2_007T++.chi")
    Qbkg, I_Qbkg = read_file("../data/cea_files/CO2/WO2_013T++.chi")

    # print(Q.shape, I_Q.shape)
    # print(Q.shape, len(I_Q))
    # print(Qbkg.shape, I_Qbkg.shape)



    # Ar
    # minQ = 3
    # maxQ = 109
    # CO2
    minQ = 8.005
    maxQ = 100
    QmaxIntegrate = 90
    # QmaxIntegrate = np.arange(60, 100, 2.5)
    # QmaxIntegrate = np.arange(90)
    
    Q, I_Q, Qbkg, I_Qbkg = check_data_length(Q, I_Q, Qbkg, I_Qbkg, minQ, maxQ)
        
    # cm
    ws1 = 0.005
    ws2 = 0.005
    r1 = 5
    r2 = 20
    d = 1
    sample_thickness = 0.005 # cm (0.04)
    dac_thickness = 0.04 # cm (0.144)
    
    # _2theta = np.degrees(Qto2theta(Q))
    _2theta = Qto2theta(Q)
    
    
    
    sth, phi_angle = calc_phi_matrix(sample_thickness/2, _2theta, ws1, ws2, r1, r2, d, 50)
    T_MCC_sample = calc_T_MCC_sample(phi_angle)
    
    sth_all, phi_angle_dac = calc_phi_matrix(dac_thickness, _2theta, ws1, ws2, r1, r2, d, 500)
    # T_MCC_DAC = calc_T_MCC_sample(phi_angle_dac)
    T_MCC_ALL, T_MCC_DAC = calc_T_MCC_DAC(phi_angle_dac, T_MCC_sample)
    
    plot3d(_2theta, sth, phi_angle, r"2$\vartheta$ (rad)", "x(cm)", r"$\varphi(2\vartheta,x)$ (rad)", "Sample")
    plot3d(_2theta, sth_all, phi_angle_dac, r"2$\vartheta$ (rad)", "x(cm)", r"$\varphi(2\vartheta,x)$ (rad)", "DAC")
    plt.show
    
    plt.figure("simps")
    plt.plot(_2theta, T_MCC_sample, label="sample")
    plt.plot(_2theta, T_MCC_ALL, label="ALL")
    plt.plot(_2theta, T_MCC_DAC, label="DAC")
    plt.xlabel(r"2$\vartheta$ (rad)")
    plt.ylabel(r"$T^{MCC}(2\vartheta, s_{th})$")
    plt.legend()
    plt.grid()
    plt.show()
    
    
    # abs_length = 1.208
    # I_Qcorr, corr_factor_meas = calc_absorption_correction(abs_length, _2theta, dac_thickness, I_Q, 0)
    # I_Qbkgcorr, corr_factor_bkg = calc_absorption_correction(abs_length, _2theta, dac_thickness, I_Qbkg, 0)
    
    # I_Qcorr10, corr_factor_meas10 = calc_absorption_correction(abs_length, _2theta, dac_thickness, I_Q, 10)
    # I_Qbkgcorr10, corr_factor_bkg10 = calc_absorption_correction(abs_length, _2theta, dac_thickness, I_Qbkg, 10)
    
    # I_Qcorr20, corr_factor_meas20 = calc_absorption_correction(abs_length, _2theta, dac_thickness, I_Q, 20)
    # I_Qbkgcorr20, corr_factor_bkg20 = calc_absorption_correction(abs_length, _2theta, dac_thickness, I_Qbkg, 20)
    
    # plt.figure("diamond correction factor")
    # plt.plot(_2theta, corr_factor_meas, label="0")
    # plt.plot(_2theta, corr_factor_meas10, label="10")
    # plt.plot(_2theta, corr_factor_meas20, label="20")
    # plt.xlabel(r"2$\theta$ (rad)")
    # plt.ylabel("Abs correction factor")
    # plt.legend()
    # plt.grid()
    # plt.show()
    
    # plt.figure("diamond correction")
    # plt.plot(Q, I_Q, label="measured")
    # plt.plot(Q, I_Qcorr, label="corr 0")
    # plt.plot(Q, I_Qcorr10, label="corr 10")
    # plt.plot(Q, I_Qcorr20, label="corr 20")
    # plt.plot(Qbkg, I_Qbkg, label="bkg")
    # plt.plot(Qbkg, I_Qbkgcorr, label="bkg corr")
    # plt.plot(Qbkg, I_Qbkgcorr20, label="bkg corr 20")
    # plt.xlabel('Q ($nm^{-1}$)')
    # plt.ylabel('I(Q)')
    # plt.legend()
    # plt.grid()
    # plt.show
    
    # plt.figure("diamond correction diff")
    # plt.plot(Q, I_Qcorr - I_Qcorr20, label="measured")
    # plt.plot(Q, I_Qcorr - I_Qcorr20, label="measured")
    # # plt.plot(Qbkg, I_Qbkgcorr - I_Qbkgcorr20, label="bkg")
    # plt.xlabel('Q ($nm^{-1}$)')
    # plt.ylabel('$\Delta$(I(Q))')
    # plt.legend()
    # plt.grid()
    # plt.show()
    
    # min_index, max_index = calc_indices(Q, minQ, QmaxIntegrate, maxQ)
    # validation_index, integration_index, calculation_index = calc_ranges(Q, minQ, QmaxIntegrate, maxQ)

    # elementList = {"Ar":1}
    # elementList = {"C":1,"O":2}

    # test values
    # Ar
    # s = np.arange(0.2, 1.0, 0.1)
    # CO2
    # s = np.arange(0.60, 1.05, 0.01)
    # rho0 = np.arange(20, 31, 1)

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

    # chi2 = np.zeros((rho0.size, s.size))

    # # remember the electron unit in atomic form factor!!!
    # fe_Q, Ztot = calc_eeff(elementList, Q)
    # Iincoh_Q = calc_Iincoh(elementList, Q)
    # J_Q = calc_JQ(Iincoh_Q, Ztot, fe_Q)
    # Sinf = calc_Sinf(elementList, fe_Q, Q, Ztot)

    # # QIsample, Isample_QIgor = read_file("../data/cea_files/CO2/WO2_007Subt.chi")

    # # best_rho0_s = np.zeros((len(QmaxIntegrate), 2))

    # # for l, val_QmaxIntegrate in enumerate(QmaxIntegrate):
    # min_index, max_index = calc_indices(Q, minQ, QmaxIntegrate, maxQ)
    # validation_index, integration_index, calculation_index = calc_ranges(Q, minQ, QmaxIntegrate, maxQ)

    # for i, val_rho0 in enumerate(rho0):
        # for j, val_s in enumerate(s):
            # Isample_Q = calc_IsampleQ(I_Q, s[j], I_Qbkg)
            # alpha = calc_alpha(J_Q[integration_index], Sinf, Q[integration_index], Isample_Q[integration_index], fe_Q[integration_index], Ztot, rho0[i])
            # Icoh_Q = calc_Icoh(numAtoms, alpha, Isample_Q, Iincoh_Q)
            # # alpha = calc_alpha(J_Q[integration_index], Sinf, Q[integration_index], Isample_QIgor[integration_index], fe_Q[integration_index], Ztot, rho0[i])
            # # Icoh_Q = calc_Icoh(numAtoms, alpha, Isample_QIgor, Iincoh_Q)

            # S_Q = calc_SQ(numAtoms, Icoh_Q, Ztot, fe_Q, Sinf, Q, max_index, integration_index)
            # newQ, S_Qsmoothed = calc_SQsmoothing(Q[validation_index], S_Q[validation_index], Sinf, smooth_factor, min_index, minQ, QmaxIntegrate, maxQ, 550)
            # S_QsmoothedDamp = calc_SQdamp(S_Qsmoothed, newQ, Sinf, QmaxIntegrate, damp_factor)

            # Qi_Q = calc_QiQ(newQ, S_QsmoothedDamp, Sinf)
            # i_Q = calc_iQ(S_QsmoothedDamp, Sinf)

            # validation_indexSmooth, integration_indexSmooth, calculation_indexSmooth = calc_ranges(newQ, minQ, QmaxIntegrate, maxQ)
            # min_indexSmooth, max_indexSmooth = calc_indices(newQ, minQ, QmaxIntegrate, maxQ)

            # r = calc_r(newQ)
            # F_r = calc_Fr(r, newQ[integration_indexSmooth], Qi_Q[integration_indexSmooth])

            # iintra_Q, fe_QSmooth = calc_iintra(newQ, max_indexSmooth, elementList, element, x, y, z)
            # Qiintradamp = calc_iintradamp(iintra_Q, newQ, QmaxIntegrate, damp_factor)
            # Fintra_r = calc_Fr(r, newQ[integration_indexSmooth], Qiintradamp[integration_indexSmooth])
            # # Fintra_r = calc_Fintra(r, newQ, QmaxIntegrate[l])


            # Iincoh_QSmooth = calc_Iincoh(elementList, newQ)
            # J_QSmooth = calc_JQ(Iincoh_QSmooth, Ztot, fe_QSmooth)

            # F_rIt = calc_optimize_Fr(iteration, F_r, Fintra_r, rho0[i], i_Q[integration_indexSmooth], newQ[integration_indexSmooth], Sinf, J_QSmooth[integration_indexSmooth], r, rmin, "n")

            # # plt.figure("F_rit")
            # # plt.plot(r, F_rit)
            # # plt.xlabel('s')
            # # plt.ylabel('chi2')
            # # plt.grid()
            # # plt.show()

            # maskIt = np.where((r>0) & (r < rmin))
            # rIt = r[maskIt]
            # deltaF_r = calc_deltaFr(F_rIt[maskIt], Fintra_r[maskIt], rIt, rho0[i])

            # chi2[i][j] = simps(deltaF_r**2, r[maskIt])
            # print(i, j, val_rho0, val_s)
            # print("chi2 ", chi2[i][j])


    # sBest, rho0Best = calc_min_chi2(s, rho0, chi2)

    # plot_chi2(s, rho0, chi2)
    # plot3d(s, rho0, chi2)
    # plt.show()












    # for QmaxIntegrate
    # best_rho0_s[l][0] = rho0[minIndxRho0]
    # best_rho0_s[l][1] = s[minIndxS]

    # print(best_rho0_s[l][0], best_rho0_s[l][1])

# plt.figure('QmaxIntegrate rho0')
# plt.plot(QmaxIntegrate, best_rho0_s[ : ,0])
# plt.axis([59, 101, 23.9, 25.1])
# plt.xlabel('QmaxIntegrate($nm^{-1}$)')
# plt.ylabel('rho0')
# plt.grid()
# plt.show

# plt.figure('QmaxIntegrate s')
# plt.plot(QmaxIntegrate, best_rho0_s[ : ,1])
# plt.axis([59, 101, 0.59, 0.81])
# plt.xlabel('QmaxIntegrate($nm^{-1}$)')
# plt.ylabel('s')
# plt.grid()
# plt.show

# best_rho0 = np.mean(best_rho0_s[ : ,0])

# plt.ion()
# for i in range(len(QmaxIntegrate)):
    # min_index, max_index = calc_indices(Q, minQ, QmaxIntegrate[i], maxQ)
    # validation_index, integration_index, calculation_index = calc_ranges(Q, minQ, QmaxIntegrate[i], maxQ)

    # Isample_Q = calc_IsampleQ(I_Q, best_rho0_s[i,1], I_Qbkg)
    # alpha = calc_alpha(J_Q[integration_index], Sinf, Q[integration_index], Isample_Q[integration_index], fe_Q[integration_index], Ztot, best_rho0)
    # Icoh_Q = calc_Icoh(N, alpha, Isample_Q, Iincoh_Q)
    # S_Q = calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, max_index, integration_index)

    # plt.figure('S_Q QmaxIntegrate')
    # plt.plot(Q[validation_index], S_Q, label='%s' %QmaxIntegrate[i])
    # plt.legend()
    # plt.draw()

    # time.sleep(1.0)

# plt.grid()
# plt.ioff()
# plt.show()
