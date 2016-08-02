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
This script is mainly used for testing the software, but it can be used to run LASDiA 
in text mode.

The nomenclature and the procedures follow the article: Eggert et al. 2002 PRB, 65, 174105.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by
an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""

from __future__ import (absolute_import, division, print_function, unicode_literals)
import six

import matplotlib.pyplot as plt

from modules import MainFunctions
from modules import Minimization
from modules import Optimization
from modules import Utility
from modules import UtilityAnalysis
from modules import Formalism
from modules import Geometry
from modules import IgorFunctions
from modules import KaplowMethod

if __name__ == '__main__':
    variables = Utility.read_inputFile("./inputFile.txt")
    
    elementList = Utility.molToelemList(variables.molecule)
    elementParameters = Utility.read_parameters(elementList, variables.element_params)
    
    path = Utility.path_xyz_file(variables.molecule)
    numAtoms, element, x, y, z = Utility.read_xyz_file(path)
    
    Q, I_Q = Utility.read_file(variables.data_file)
    Qbkg, Ibkg_Q  = Utility.read_file(variables.bkg_file)
    
    Q, I_Q, Qbkg, Ibkg_Q = UtilityAnalysis.check_data_length(Q, I_Q, Qbkg, Ibkg_Q , \
        variables.minQ, variables.maxQ)
    
    # chi2 = np.zeros((rho0.size, s.size))
    
    aff_path = "C:/Users/devoto/work/ID27/LASDiA/affParamCEA.txt"
    incoh_path = "C:/Users/devoto/work/ID27/LASDiA/incohParamCEA.txt" 
    
    fe_Q, Ztot = MainFunctions.calc_eeff(elementList, Q, elementParameters)
    fe_Q2, Ztot2 = MainFunctions.calc_eeff2(elementList, Q, incoh_path, aff_path)
    
    # Utility.plot_data(Q, fe_Q, "fe_Q", r"$Q(nm^{-1})$", r"$fe(Q)$", r"$fe(Q)$", "y")
    Utility.plot_data(Q, fe_Q-fe_Q2, "fe_Q", r"$Q(nm^{-1})$", r"$fe(Q)$", r"$fe(Q)2$", "y")
    
    print(MainFunctions.calc_Kp(fe_Q, "C", Q, elementParameters))
    print(MainFunctions.calc_Kp2(fe_Q, "C", Q, aff_path))
    
    
    
    
    Iincoh_Q = MainFunctions.calc_Iincoh(elementList, Q, elementParameters)
    J_Q = MainFunctions.calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf = MainFunctions.calc_Sinf(elementList, fe_Q, Q, Ztot, elementParameters)
    
    s = variables.s_value
    rho0 = variables.rho0_value
    
    s_array = UtilityAnalysis.make_array(variables.s_value, 20)
    rho0_array = UtilityAnalysis.make_array(variables.rho0_value, 20)
    
    
    # for rho0 in rho0_array:
        # for s in s_array:
    # S_Q, r, F_r = KaplowMethod.Kaplow_method(numAtoms, variables, Q, I_Q, \
        # Ibkg_Q, J_Q, fe_Q, Iincoh_Q, Sinf, Ztot, s, rho0)
        
    # Utility.plot_data(Q, S_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S(Q)$", "y")
    # Utility.plot_data(r, F_r, "F_r", r"$r(nm)$", r"$F(r)$", r"$F(r)$", "y")
    
    plt.show()
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
         
            
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
