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
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np
import time
import math

# from modules.MainFunctions import *
# from modules.Utility import *
# from modules.UtilityAnalysis import *
# from modules.Optimization import *
# from modules.Minimization import *
# from modules.Formalism import *
# from modules.IgorFunctions import *
# from modules.Geometry import *

from modules import MainFunctions
from modules import Utility
from modules import UtilityAnalysis
from modules import Optimization
from modules import Minimization
from modules import Formalism
from modules import IgorFunctions
from modules import Geometry


# import cmath
# from cmath import exp, pi

if __name__ == '__main__':
    molecule = "CO2"
    elementList = Utility.molToelemList(molecule)
    elementParameters = Utility.read_parameters(elementList, "./elementParameters.txt")

    smooth_factor = 0.25
    damp_factor = 1
    iteration = 2
    rmin = 0.22
    numAtoms, element, x, y, z = Utility.read_xyz_file("./xyzFiles/co2.xyz")
    aff_params = "./affParamCEA.txt"
    incoh_params = "./incohParamCEA.txt"

    Q, I_Q = Utility.read_file("../data/cea_files/CO2/WO2_007T++.chi")
    Qbkg, Ibkg_Q  = Utility.read_file("../data/cea_files/CO2/WO2_013T++.chi")

    # Ar
    # minQ = 3
    # maxQ = 109
    # CO2
    minQ = 8.005
    maxQ = 100
    QmaxIntegrate = 90
    # QmaxIntegrate = np.arange(60, 100, 2.5)
    # QmaxIntegrate = np.arange(90)
    
    Q, I_Q, Qbkg, Ibkg_Q  = UtilityAnalysis.check_data_length(Q, I_Q, Qbkg, Ibkg_Q , minQ, maxQ)
    
    min_index, max_index = UtilityAnalysis.calc_indices(Q, minQ, QmaxIntegrate, maxQ)
    validation_index, integration_index, calculation_index = UtilityAnalysis.calc_ranges(Q, minQ, QmaxIntegrate, maxQ)
    
    
    # cm
    ws1 = 0.005
    ws2 = 0.005
    r1 = 5
    r2 = 20
    d = 1
    sth = np.arange(0.02, 0.06, 0.01)
    dac_thickness = 0.144 # 0.04 # cm
    phi_matrix_thickness = 0.17 #0.04 # cm

    # test values
    # Ar
    # s = np.arange(0.2, 1.0, 0.1)
    # CO2
    # s = np.arange(0.7, 1.0, 0.5)
    # rho0 = np.arange(25, 30, 1)
    # sth = np.arange(0.002, 0.004, 0.001)

    # best values
    # Ar
    # s = np.array([0.57])
    # rho0 = np.array([26.1])
    # CO2
    s = np.array([0.984228])
    rho0 = np.array([29.6625])
    
    s_value = 0.98
    rho0_value = 29

    chi2 = np.zeros((rho0.size, s.size))

    fe_Q, Ztot = MainFunctions.calc_eeff(elementList, Q, elementParameters)
    Iincoh_Q = MainFunctions.calc_Iincoh(elementList, Q, elementParameters)
    J_Q = MainFunctions.calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf = MainFunctions.calc_Sinf(elementList, fe_Q, Q, Ztot, elementParameters)
    
    # iintra_Q2 = Optimization.calc_iintra2(Q, fe_Q, Ztot, QmaxIntegrate, maxQ, elementList, element, x, y, z, elementParameters)
    
    
    # while True:
        # Isample_Q = MainFunctions.calc_IsampleQ(I_Q, s[j], Ibkg_Q )
        # alpha = MainFunctions.calc_alpha(J_Q[Q<=QmaxIntegrate], Sinf, Q[Q<=QmaxIntegrate], \
            # Isample_Q[Q<=QmaxIntegrate], fe_Q[Q<=QmaxIntegrate], Ztot, rho0[i])
        # Icoh_Q = MainFunctions.calc_Icoh(numAtoms, alpha, Isample_Q, Iincoh_Q)
        # S_Q = MainFunctions.calc_SQ(numAtoms, Icoh_Q, Ztot, fe_Q, Sinf, Q, minQ, QmaxIntegrate, maxQ)
        
        # if ...:
            # break
    
    
    
    
    
    for i, val_rho0 in enumerate(rho0):
        for j, val_s in enumerate(s):
            Isample_Q = MainFunctions.calc_IsampleQ(I_Q, s[j], Ibkg_Q )
            alpha = MainFunctions.calc_alpha(J_Q[Q<=QmaxIntegrate], Sinf, Q[Q<=QmaxIntegrate], \
                Isample_Q[Q<=QmaxIntegrate], fe_Q[Q<=QmaxIntegrate], Ztot, rho0[i])
            Icoh_Q = MainFunctions.calc_Icoh(numAtoms, alpha, Isample_Q, Iincoh_Q)
            S_Q = MainFunctions.calc_SQ(numAtoms, Icoh_Q, Ztot, fe_Q, Sinf, Q, minQ, QmaxIntegrate, maxQ)
            
            
            alphaFZ = Formalism.calc_alphaFZ(numAtoms, Q, Isample_Q, Iincoh_Q, rho0, elementParameters)
            IcohFZ_Q = MainFunctions.calc_Icoh(numAtoms, alphaFZ, Isample_Q, Iincoh_Q)
            SFZ_Q = Formalism.calc_S_QFZ(numAtoms, IcohFZ_Q, Ztot, Q, elementParameters)
            
            Utility.plot_data(Q, S_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S(Q)$", "y")
            Utility.plot_data(Q, SFZ_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S_{FZ}(Q)$", "y")
            plt.show()
            
            
            # newQ, S_Qsmoothed = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, smooth_factor, minQ, QmaxIntegrate, maxQ, 550)
            # S_QsmoothedDamp = UtilityAnalysis.calc_SQdamp(S_Qsmoothed, newQ, Sinf, QmaxIntegrate, damp_factor)
            
            # newfe_Q = UtilityAnalysis.interpolation_after_smoothing(Q, newQ, fe_Q)
            
            # i_Q = MainFunctions.calc_iQ(S_QsmoothedDamp, Sinf)
            
            # r = MainFunctions.calc_r(newQ)
            # F_r = MainFunctions.calc_Fr(r, newQ[newQ<=QmaxIntegrate], i_Q[newQ<=QmaxIntegrate])
            # # Utility.plot_data(r, F_r, "F_r", r"$r(nm)$", r"F(r)", "F(r)", "n")
            # # plt.show()
            
            # min_indexSmooth, max_indexSmooth = UtilityAnalysis.calc_indices(newQ, minQ, QmaxIntegrate, maxQ)
            # iintra_Q, fe_QSmooth = Optimization.calc_iintra(newQ, max_indexSmooth, elementList, element, x, y, z, incoh_params, aff_params)
            # # Qiintradamp = calc_iintradamp(iintra_Q, newQ, QmaxIntegrate, damp_factor)
            # # Fintra_r = calc_Fr(r, newQ[integration_indexSmooth], Qiintradamp[integration_indexSmooth])
            # # Fintra_r = calc_Fintra(r, newQ, QmaxIntegrate)
            
            # print(iintra_Q.size)
            # print(iintra_Q2.size)
            
            # newiintra_Q2 = UtilityAnalysis.interpolation_after_smoothing(Q, newQ, iintra_Q2)
            # print(newiintra_Q2.size)
            
            # Utility.plot_data(newQ, iintra_Q, "iintra_Q", r"$Q(nm^{-1})$", r"$i_{intra}(Q)$", r"$i_{intra}(Q)$", "y")
            # Utility.plot_data(Q, iintra_Q2, "iintra_Q", r"$Q(nm^{-1})$", r"$i_{intra}(Q)$", r"$i_{intra}(Q)2$", "y")
            # Utility.plot_data(newQ, newiintra_Q2, "iintra_Q", r"$Q(nm^{-1})$", r"$i_{intra}(Q)$", r"new$i_{intra}(Q)2$", "y")
            # plt.show()
            
            # Iincoh_QSmooth = calc_Iincoh(elementList, newQ, incoh_params, aff_params)
            # J_QSmooth = calc_JQ(Iincoh_QSmooth, Ztot, fe_QSmooth)
            #
            # F_rIt = calc_optimize_Fr(iteration, F_r, Fintra_r, rho0[i], i_Q[integration_indexSmooth], newQ[integration_indexSmooth], Sinf, J_QSmooth[integration_indexSmooth], r, rmin, "n")
            #
            # maskIt = np.where((r>0) & (r < rmin))
            # rIt = r[maskIt]
            # deltaF_r = calc_deltaFr(F_rIt[maskIt], Fintra_r[maskIt], rIt, rho0[i])

            # # # chi2[i][j][k] = simps(deltaF_r**2, r[maskIt])
            # print(i, j, val_rho0, val_s)
            # print("chi2 ", chi2[i][j])


    
    # # # # chi2Min, sBest, sBestIdx, rho0Best, rho0BestIdx, sthBest, sthBestIdx = calc_min_chi2(s, rho0, sth, chi2)

    # # # # print(sBest, rho0Best, sthBest)




















    # Geometrical correction
    # two_theta = UtilityAnalysis.Qto2theta(Q) # rad
    
    # abs_length = 1.208
    # corr_factor_meas0 = Geometry.calc_absorption_correction(abs_length, two_theta, dac_thickness, 0)
    # I_Q = I_Q / corr_factor_bkg
    # Ibkg_Q  = Ibkg_Q  / corr_factor_bkg

    # all_thickness_sampling, phi_matrix = Geometry.calc_phi_matrix(phi_matrix_thickness, two_theta, ws1, ws2, r1, r2, d, 1000)

    # Utility.plot3d(two_theta, all_thickness_sampling, phi_matrix, "phi_matrix", "2theta", "x(cm)", "phi")

    # T_MCC_sample1, T_MCC_DAC1, T_MCC_ALL1 = Geometry.calc_T_MCC(0.001, all_thickness_sampling, phi_matrix, "y")
    # T_MCC_sample2, T_MCC_DAC2, T_MCC_ALL2 = Geometry.calc_T_MCC(0.002, all_thickness_sampling, phi_matrix, "y")
    # T_MCC_sample3, T_MCC_DAC3, T_MCC_ALL3 = Geometry.calc_T_MCC(0.003, all_thickness_sampling, phi_matrix, "y")
    # T_MCC_sample4, T_MCC_DAC4, T_MCC_ALL4 = Geometry.calc_T_MCC(0.004, all_thickness_sampling, phi_matrix, "y")

    # Utility.write_file("./T_MCC_sample4a.txt", Q, T_MCC_sample4, "Q", "T_MCC_sample4")
    # Utility.write_file("./T_MCC_DAC4a.txt", Q, T_MCC_DAC4, "Q", "T_MCC_DAC4")

    # Utility.plot_data(Q, T_MCC_ALL1, "T_MCC_all", r"$Q (nm^{-1})$", r"$T^{MCC}_{ALL}(2\vartheta, s_{th})$", "0.001 cm", "y")
    # Utility.plot_data(Q, T_MCC_ALL2, "T_MCC_all", r"$Q (nm^{-1})$", r"$T^{MCC}_{ALL}(2\vartheta, s_{th})$", "0.002 cm", "y")
    # Utility.plot_data(Q, T_MCC_ALL3, "T_MCC_all", r"$Q (nm^{-1})$", r"$T^{MCC}_{ALL}(2\vartheta, s_{th})$", "0.003 cm", "y")
    # Utility.plot_data(Q, T_MCC_ALL4, "T_MCC_all", r"$Q (nm^{-1})$", r"$T^{MCC}_{ALL}(2\vartheta, s_{th})$", "0.004 cm", "y")

    # Utility.plot_data(Q, T_MCC_sample1, "T_MCC_sample", r"$Q (nm^{-1})$", r"$T^{MCC}_{sample}(2\vartheta, s_{th})$", "0.001 cm", "y")
    # Utility.plot_data(Q, T_MCC_sample2, "T_MCC_sample", r"$Q (nm^{-1})$", r"$T^{MCC}_{sample}(2\vartheta, s_{th})$", "0.002 cm", "y")
    # Utility.plot_data(Q, T_MCC_sample3, "T_MCC_sample", r"$Q (nm^{-1})$", r"$T^{MCC}_{sample}(2\vartheta, s_{th})$", "0.003 cm", "y")
    # Utility.plot_data(Q, T_MCC_sample4, "T_MCC_sample", r"$Q (nm^{-1})$", r"$T^{MCC}_{sample}(2\vartheta, s_{th})$", "0.004 cm", "y")

    # Utility.plot_data(Q, T_MCC_DAC1, "T_MCC_DAC", r"$Q (nm^{-1})$", r"$T^{MCC}_{DAC}(2\vartheta, s_{th})$", "0.001 cm", "y")
    # Utility.plot_data(Q, T_MCC_DAC2, "T_MCC_DAC", r"$Q (nm^{-1})$", r"$T^{MCC}_{DAC}(2\vartheta, s_{th})$", "0.002 cm", "y")
    # Utility.plot_data(Q, T_MCC_DAC3, "T_MCC_DAC", r"$Q (nm^{-1})$", r"$T^{MCC}_{DAC}(2\vartheta, s_{th})$", "0.003 cm", "y")
    # Utility.plot_data(Q, T_MCC_DAC4, "T_MCC_DAC", r"$Q (nm^{-1})$", r"$T^{MCC}_{DAC}(2\vartheta, s_{th})$", "0.004 cm", "y")

    # Utility.plot_data(Q, T_MCC_ALL4,    "T_MCC", r"$Q (nm^{-1})$", r"$T^{MCC}(2\vartheta, s_{th})$", "ALL", "y")
    # Utility.plot_data(Q, T_MCC_sample4, "T_MCC", r"$Q (nm^{-1})$", r"$T^{MCC}(2\vartheta, s_{th})$", "Sample", "y")
    # Utility.plot_data(Q, T_MCC_DAC4,    "T_MCC", r"$Q (nm^{-1})$", r"$T^{MCC}{2\vartheta, s_{th})$", "DAC", "y")

    # T_MCC_corr_factor_bkg1 = Geometry.calc_T_DAC_MCC_bkg_corr(T_MCC_DAC1, T_MCC_DAC3)
    # T_MCC_corr_factor_bkg2 = Geometry.calc_T_DAC_MCC_bkg_corr(T_MCC_DAC2, T_MCC_DAC3)
    # T_MCC_corr_factor_bkg3 = Geometry.calc_T_DAC_MCC_bkg_corr(T_MCC_DAC3, T_MCC_DAC3)
    # T_MCC_corr_factor_bkg4 = Geometry.calc_T_DAC_MCC_bkg_corr(T_MCC_DAC4, T_MCC_DAC3)

    # Utility.write_file("./T_MCC_corr_factor_bkg4a.txt", Q, T_MCC_corr_factor_bkg4, "Q", "T_MCC_corr_factor_bkg4")

    # Utility.plot_data(Q, T_MCC_corr_factor_bkg1, "T_MCC_bkg", r"$Q (nm^{-1})$", r"$T^{MCC}(2\vartheta, s_{th})/T^{MCC}(2\vartheta, s_{th})$", "0.001 cm", "y")
    # Utility.plot_data(Q, T_MCC_corr_factor_bkg2, "T_MCC_bkg", r"$Q (nm^{-1})$", r"$T^{MCC}(2\vartheta, s_{th})/T^{MCC}(2\vartheta, s_{th})$", "0.002 cm", "y")
    # Utility.plot_data(Q, T_MCC_corr_factor_bkg3, "T_MCC_bkg", r"$Q (nm^{-1})$", r"$T^{MCC}(2\vartheta, s_{th})/T^{MCC}(2\vartheta, s_{th})$", "0.003 cm", "y")
    # Utility.plot_data(Q, T_MCC_corr_factor_bkg4, "T_MCC_bkg", r"$Q (nm^{-1})$", r"$T^{MCC}(2\vartheta, s_{th})/T^{MCC}(2\vartheta, s_{th})$", "0.004 cm", "y")

    # plt.show()

    

















    # # QIsample, Isample_QIgor = read_file("../data/cea_files/CO2/WO2_007Subt.chi")

    # # best_rho0_s = np.zeros((len(QmaxIntegrate), 2))

    # # for l, val_QmaxIntegrate in enumerate(QmaxIntegrate):
    # min_index, max_index = calc_indices(Q, minQ, QmaxIntegrate, maxQ)
    # validation_index, integration_index, calculation_index = calc_ranges(Q, minQ, QmaxIntegrate, maxQ)

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

    # Isample_Q = calc_IsampleQ(I_Q, best_rho0_s[i,1], Ibkg_Q )
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
