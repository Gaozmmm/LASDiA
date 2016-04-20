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

import matplotlib.pyplot as plt

from modules.MainFunctions import *
# from modules.Utility import *
# from modules.UtilityAnalysis import *
from modules.Optimization import *
from modules.Minimization import *
# from modules.Formalism import *
# from modules.IgorFunctions import *


if __name__ == "__main__":
    
    variables = read_inputFile("./inputFile.txt")
    
    elementList = molToelemList(variables.molecule)
    numAtoms, element, x, y, z = read_xyz_file(variables.xyz_file)
    
    Q, I_Q = read_file(variables.data_file)
    Qbkg, I_Qbkg = read_file(variables.bkg_file)
    
    plt.ion()
    if variables.pw_raw_data[0].lower() == "y":
        plot_raw_data(Q, I_Q, "raw_data", "Q($nm^{-1}$)", "I(Q)", "I(Q)")
        plot_raw_data(Qbkg, I_Qbkg, "raw_data", "Q($nm^{-1}$)", "I(Q)", "I(Q) bkg")
    
    fe_Q, Ztot = calc_eeff(elementList, Q, variables.incoh_params, variables.aff_params)
    Iincoh_Q = calc_Iincoh(elementList, Q, variables.incoh_params, variables.aff_params)
    J_Q = calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf = calc_Sinf(elementList, fe_Q, Q, Ztot, variables.aff_params)

    min_index, max_index = calc_indices(Q, variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    validation_index, integration_index, calculation_index = calc_ranges(Q, variables.minQ, variables.QmaxIntegrate, variables.maxQ)

    scale_factor = setArray(variables.s_min, variables.s_max, variables.s_step)
    rho0 = setArray(variables.rho0_min, variables.rho0_max, variables.rho0_step)

    chi2 = np.zeros((rho0.size, scale_factor.size))
    valChi2 = 100000.0

    for i, val_rho0 in enumerate(rho0):
        for j, val_s in enumerate(scale_factor):
            # print(i, j, val_rho0, val_s)
            Isample_Q = calc_IsampleQ(I_Q, scale_factor[j], I_Qbkg)
            alpha = calc_alpha(J_Q[integration_index], Sinf, Q[integration_index], Isample_Q[integration_index], fe_Q[integration_index], Ztot, rho0[i])
            Icoh_Q = calc_Icoh(numAtoms, alpha, Isample_Q, Iincoh_Q)

            S_Q = calc_SQ(numAtoms, Icoh_Q, Ztot, fe_Q, Sinf, Q, max_index, integration_index)
            newQ, S_Qsmoothed = calc_SQsmoothing(Q[validation_index], S_Q[validation_index], Sinf, variables.smooth_factor, min_index, variables.minQ, variables.QmaxIntegrate, variables.maxQ, 550)
            S_QsmoothedDamp = calc_SQdamp(S_Qsmoothed, newQ, Sinf, variables.QmaxIntegrate, variables.damp_factor)
            
            if variables.pw_S_Q[0].lower() == "y":
                plot_data(Q[validation_index], S_Q, "S_Q", "Q($nm^{-1}$)", "S(Q)", "S(Q)")
            # if variables.pw_S_Q[1].lower() == "y":
                # write_file(variables.pw_S_Q[2], Q[validation_index], S_Q, "Q($nm^{-1}$)", "S(Q)")
            # if variables.pw_S_Qsmoothed_damped[0].lower() == "y":
                # plot_data(newQ, S_QsmoothedDamp, "S_Q", "Q($nm^{-1}$)", "S(Q)", "S(Q) smooth-damp")
            # if variables.pw_S_Qsmoothed_damped[1].lower() == "y":
                # write_file(variables.pw_S_Qsmoothed_damped[2], newQ, S_QsmoothedDamp, "Q($nm^{-1}$)", "S(Q) smooth-damp")
            
            Qi_Q = calc_QiQ(newQ, S_QsmoothedDamp, Sinf)
            i_Q = calc_iQ(S_QsmoothedDamp, Sinf)
            
            validation_indexSmooth, integration_indexSmooth, calculation_indexSmooth = calc_ranges(newQ, variables.minQ, variables.QmaxIntegrate, variables.maxQ)
            min_indexSmooth, max_indexSmooth = calc_indices(newQ, variables.minQ, variables.QmaxIntegrate, variables.maxQ)
            
            r = calc_r(newQ)
            F_r = calc_Fr(r, newQ[integration_indexSmooth], Qi_Q[integration_indexSmooth])
            
            if variables.pw_F_r[0].lower() == "y":
                plot_data(r, F_r, "F_r", "r($nm$)", "F(r)", "F(r)")
            # if variables.pw_F_r[1].lower() == "y":
                # write_file(variables.pw_F_r[2], r, F_r, "r($nm$)", "F(r)")
            
            iintra_Q, fe_QSmooth = calc_iintra(newQ, max_indexSmooth, elementList, element, x, y, z, variables.incoh_params, variables.aff_params)
            Qiintradamp = calc_iintradamp(iintra_Q, newQ, variables.QmaxIntegrate, variables.damp_factor)
            Fintra_r = calc_Fr(r, newQ[integration_indexSmooth], Qiintradamp[integration_indexSmooth])
            # Fintra_r = calc_Fintra(r, newQ, QmaxIntegrate[l], variables.aff_params)
            
            Iincoh_QSmooth = calc_Iincoh(elementList, newQ, variables.incoh_params, variables.aff_params)
            J_QSmooth = calc_JQ(Iincoh_QSmooth, Ztot, fe_QSmooth)

            F_rIt = calc_optimize_Fr(variables.iteration, F_r, Fintra_r, rho0[i], i_Q[integration_indexSmooth], newQ[integration_indexSmooth], Sinf, J_QSmooth[integration_indexSmooth], r, variables.rmin, variables.pw_F_r_iter[0])
            
            chi2[i][j] = calc_chi2(r, variables.rmin, F_rIt, Fintra_r, rho0[i])
            if chi2[i][j] < valChi2:
                valChi2 = chi2[i][j]
                F_rOpt = F_rIt
                # print(i, j, val_rho0, val_s)
                # print("chi2 ", chi2[i][j])
            
            # if variables.pw_F_rIt[0].lower() == "y":
                # plot_data(r, F_rIt, "F_rIt", "r($nm$)", "F(r)", "F(r)")
            # if variables.pw_F_rIt[1].lower() == "y":
                # write_file(variables.pw_F_rIt[2], r, F_rIt, "r($nm$)", "F(r)")
    
    
    if variables.pw_F_rOpt[0].lower() == "y":
        plot_data(r, F_rOpt, "F_rOpt", "r($nm$)", "F(r)", "F(r) Optimal")
    if variables.pw_F_rOpt[0].lower() == "y":
        write_file(variables.pw_F_rOpt[2], r, F_rOpt, "r($nm$)", "F(r)")
    
    S_QCorr = calc_SQCorr(F_rOpt, r, newQ, Sinf)
    
    if variables.pw_S_QCorr[0].lower() == "y":
        plot_data(newQ, S_QCorr, "S_QCorr", "Q($nm^{-1}$)", "S(Q)", "S(Q) Corrected")
    if variables.pw_S_QCorr[0].lower() == "y":
        write_file(variables.pw_S_QCorr[2], newQ, S_QCorr, "Q($nm^{-1}$))", "S(Q)")
    
    chi2Min, sBest, rho0Best = calc_min_chi2(scale_factor, rho0, chi2)
    if variables.pw_results[0].lower() == "y":
        write_results(variables.pw_results[1], variables.molecule, sBest, rho0Best)
        print("chi2 min ", chi2Min)
        print("Best scale factor ", sBest)
        print("Best rho0 density ", rho0Best)
    else:
        print("chi2 min ", chi2Min)
        print("Best scale factor ", sBest)
        print("Best rho0 density ", rho0Best)
    
    if variables.pw_chi2[0].lower() == "y":
        plot_chi2(scale_factor, rho0, chi2)
        plot3d_chi2(scale_factor, rho0, chi2)
    
    plt.ioff()
    plt.show()