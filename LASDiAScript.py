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
import numpy as np
from itertools import product
from timeit import default_timer as timer
from scipy.integrate import simps

from modules import Formalism
from modules import Geometry
from modules import IgorFunctions
from modules import KaplowMethod
from modules import MainFunctions
from modules import Minimization
from modules import Optimization
from modules import Utility
from modules import UtilityAnalysis


if __name__ == '__main__':
    variables = Utility.read_inputFile("./inputFile.txt")
    
    elementList = Utility.molToelemList(variables.molecule)
    elementParameters = Utility.read_parameters(elementList, variables.element_params)
    
    path = Utility.path_xyz_file(variables.molecule)
    numAtoms, element, x, y, z = Utility.read_xyz_file(path)
    
    Q, I_Q = Utility.read_file(variables.data_file)
    Qbkg, Ibkg_Q  = Utility.read_file(variables.bkg_file)
    
    Q, I_Q, Qbkg, Ibkg_Q = UtilityAnalysis.check_data_length(Q, I_Q, Qbkg, Ibkg_Q, \
        variables.minQ, variables.maxQ)
    
    I_Q, Ibkg_Q = Geometry.geometry_correction(Q, I_Q, Qbkg, Ibkg_Q, variables, \
        variables.phi_matrix_flag)
    
    fe_Q, Ztot = MainFunctions.calc_eeff(elementList, Q, elementParameters)
    Iincoh_Q = MainFunctions.calc_Iincoh(elementList, Q, elementParameters)
    J_Q = MainFunctions.calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf, Sinf_Q = MainFunctions.calc_Sinf(elementList, fe_Q, Q, Ztot, elementParameters)
    
    # Utility.plot_data(Q, Sinf_Q, "Sinf_Q", r"$Q(nm^{-1})$", r"$S_\infty(Q)$", r"$S_\infty(Q)$", "y")
    
    r, Fintra_r = Optimization.calc_intraComponent(Q, fe_Q, Ztot, variables.QmaxIntegrate, \
        variables.maxQ, elementList, element, x, y, z, elementParameters, variables.damping_factor)
    
    scale_factor = variables.sf_value
    density = variables.rho0_value
    percentage = 20
    
    #-----------------------------------------------------
    # Test Kaplow method for AL Morard
    
    chi2_E, S_Q_E, F_r_E = KaplowMethod.Kaplow_method(numAtoms, variables, \
                Q, I_Q, Ibkg_Q, J_Q, fe_Q, Iincoh_Q, Sinf, Ztot, scale_factor, density, Fintra_r, r)
    
    aff_sq_mean = Formalism.calc_aff_squared_mean(numAtoms, elementList, Q, elementParameters)
    aff_mean_sq = Formalism.calc_aff_mean_squared(numAtoms, elementList, Q, elementParameters)
    
    chi2AL, S_Q_AL, F_r_AL = KaplowMethod.Kaplow_methodWAL(numAtoms, variables, \
                Q, I_Q, Ibkg_Q, aff_sq_mean, aff_mean_sq, Iincoh_Q, scale_factor, density, Fintra_r, r)
    
    # Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scale_factor, Ibkg_Q)
    # alpha = Formalism.calc_alphaM(Q, Isample_Q, Iincoh_Q, density, aff_sq_mean, aff_mean_sq)
    # Icoh_Q = MainFunctions.calc_Icoh(numAtoms, alpha, Isample_Q, Iincoh_Q)
    
    # S_Q = Formalism.calc_SAL_Q(Q, Icoh_Q, aff_sq_mean, variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    
    # Utility.plot_data(Q, S_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{AL}_M(Q)$", "y")
    
    # S_Qsmoothed = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, 1, 25, \
        # variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    # S_QsmoothedDamp = UtilityAnalysis.calc_SQdamp(S_Qsmoothed, Q, 1, \
        # variables.QmaxIntegrate, variables.damping_factor)
    
    # Utility.plot_data(Q, S_QsmoothedDamp, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{AL}_M(Q)SD$", "y")
    
    # i_Q = MainFunctions.calc_iQ(S_QsmoothedDamp, 1)
    # F_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], \
        # i_Q[Q<=variables.QmaxIntegrate])
    
    # Utility.plot_data(r, F_r, "F_r", r"$r(nm)$", r"$F(r)$", r"$F^{AL}_M(r)$", "y")
    
    # J_Q = Iincoh_Q/aff_sq_mean
    
    # F_rIt, deltaF_rIt = Optimization.calc_optimize_Fr(variables.iteration, F_r, \
        # Fintra_r, density, i_Q[Q<=variables.QmaxIntegrate], Q[Q<=variables.QmaxIntegrate], \
        # 1, J_Q[Q<=variables.QmaxIntegrate], r, variables.rmin, variables.plot_iter)
    
    # Utility.plot_data(r, F_rIt, "F_r", r"$r(nm)$", r"$F(r)$", r"$F^{AL}_M(r)O$", "y")
    
    # chi2 = simps(deltaF_rIt[r < variables.rmin]**2, r[r < variables.rmin])
    
    print(chi2_E, chi2)
    
    #-----------------------------------------------------
    # Automatic loop for rho0 and sf
    
    # This is good, but I still have to test several parts and I need more flexibility!
    
    jj = 0
    # FOpt_r = []
    # Sbest_Q = []
    
    # FOpt_rAL = []
    # Sbest_QAL = []
    
    # FOpt_rFZ = []
    # Sbest_QFZ = []
    
    # sf_loop = "n"
    # rho0_loop = "n"
    
    # aff_sq_mean = Formalism.calc_aff_squared_mean(numAtoms, elementList, Q, elementParameters)
    # aff_mean_sq = Formalism.calc_aff_mean_squared(numAtoms, elementList, Q, elementParameters)
    
    # while True:
        # # print(jj)
        
        # if variables.sf_loop == "y":
            # sf_array = UtilityAnalysis.make_array_loop(scale_factor, percentage)
        # else:
            # sf_array = np.array([scale_factor])
        # if variables.rho0_loop == "y":
            # rho0_array = UtilityAnalysis.make_array_loop(density, percentage)
        # else:
            # rho0_array = np.array([density])
        
        # chi2 = np.zeros((rho0_array.size, sf_array.size))
        # chi2AL = np.zeros((rho0_array.size, sf_array.size))
        # chi2FZ = np.zeros((rho0_array.size, sf_array.size))
        # chi2ALM = np.zeros((rho0_array.size, sf_array.size))
        # chi2FZM = np.zeros((rho0_array.size, sf_array.size))
        
        # chi2_min = 1000000
        
        # for (idx_rho0, rho0), (idx_sf, sf) in product(enumerate(rho0_array), enumerate(sf_array)):
            # chi2[idx_rho0][idx_sf], S_Q, F_r = KaplowMethod.Kaplow_method(numAtoms, variables, \
                # Q, I_Q, Ibkg_Q, J_Q, fe_Q, Iincoh_Q, Sinf, Ztot, sf, rho0, Fintra_r, r)
            
            # if chi2[idx_rho0][idx_sf] < chi2_min:
                # chi2_min = chi2[idx_rho0][idx_sf]
                # best_sf = sf
                # best_rho0 = rho0
                # FOpt_r = F_r
                # Sbest_Q = S_Q
                # # FOpt_r.append(F_r)
                # # Sbest_Q.append(S_Q)
        
        # scale_factor, density = Minimization.calc_min_chi2(sf_array, rho0_array, chi2)
        
        # chi2_min = 1000000
        
        # for (idx_rho0, rho0), (idx_sf, sf) in product(enumerate(rho0_array), enumerate(sf_array)):
            # chi2AL[idx_rho0][idx_sf], S_Q, F_r = KaplowMethod.Kaplow_methodWAL(numAtoms, variables, \
                # Q, I_Q, Ibkg_Q, aff_sq_mean, aff_mean_sq, Iincoh_Q, sf, rho0, Fintra_r, r)
            # if chi2AL[idx_rho0][idx_sf] < chi2_min:
                # chi2_min = chi2AL[idx_rho0][idx_sf]
                # best_sfAL = sf
                # best_rho0AL = rho0
                # FOpt_rAL = F_r
                # Sbest_QAL = S_Q
        
        # scale_factorAL, densityAL = Minimization.calc_min_chi2(sf_array, rho0_array, chi2AL)
        
        # chi2_min = 1000000
        
        # for (idx_rho0, rho0), (idx_sf, sf) in product(enumerate(rho0_array), enumerate(sf_array)):
            # chi2FZ[idx_rho0][idx_sf], S_Q, F_r = KaplowMethod.Kaplow_methodWFZ(numAtoms, variables, \
                # Q, I_Q, Ibkg_Q, aff_sq_mean, aff_mean_sq, Iincoh_Q, sf, rho0, Fintra_r, r)
            # if chi2FZ[idx_rho0][idx_sf] < chi2_min:
                # chi2_min = chi2FZ[idx_rho0][idx_sf]
                # best_sfFZ = sf
                # best_rho0FZ = rho0
                # FOpt_rFZ = F_r
                # Sbest_QFZ = S_Q
        
        # scale_factorFZ, densityFZ = Minimization.calc_min_chi2(sf_array, rho0_array, chi2FZ)
        
        # chi2_min = 1000000
        
        # for (idx_rho0, rho0), (idx_sf, sf) in product(enumerate(rho0_array), enumerate(sf_array)):
            # chi2ALM[idx_rho0][idx_sf], S_Q, F_r = KaplowMethod.Kaplow_methodMAL(numAtoms, variables, \
                # Q, I_Q, Ibkg_Q, aff_sq_mean, aff_mean_sq, Iincoh_Q, sf, rho0, Fintra_r, r)
            # print(chi2ALM[idx_rho0][idx_sf], chi2_min)
            # if chi2ALM[idx_rho0][idx_sf] < chi2_min:
                # chi2_min = chi2ALM[idx_rho0][idx_sf]
                # best_sfALM = sf
                # best_rho0ALM = rho0
                # FOpt_rALM = F_r
                # Sbest_QALM = S_Q
        
        # print(Sbest_QALM)
        
        # scale_factorALM, densityALM = Minimization.calc_min_chi2(sf_array, rho0_array, chi2ALM)
        
        # chi2_min = 1000000
        
        # for (idx_rho0, rho0), (idx_sf, sf) in product(enumerate(rho0_array), enumerate(sf_array)):
            # chi2FZM[idx_rho0][idx_sf], S_Q, F_r = KaplowMethod.Kaplow_methodMFZ(numAtoms, variables, \
                # Q, I_Q, Ibkg_Q, aff_sq_mean, aff_mean_sq, Iincoh_Q, sf, rho0, Fintra_r, r)
            # if chi2FZM[idx_rho0][idx_sf] < chi2_min:
                # chi2_min = chi2FZM[idx_rho0][idx_sf]
                # best_sfFZM = sf
                # best_rho0FZM = rho0
                # FOpt_rFZM = F_r
                # Sbest_QFZM = S_Q
        
        # scale_factorFZM, densityFZM = Minimization.calc_min_chi2(sf_array, rho0_array, chi2FZM)
        
        # percentage /= 2
        # jj += 1
        # if jj == 2:
            # break
    
    # Utility.plot_data(Q, Sbest_Q, "Sbest_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{E}(Q)$", "y")
    # Utility.plot_data(Q, Sbest_QAL, "Sbest_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{AL}(Q)$", "y")
    # Utility.plot_data(Q, Sbest_QFZ, "Sbest_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{FZ}(Q)$", "y")
    # Utility.plot_data(Q, Sbest_QALM, "Sbest_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{AL}_M(Q)$", "y")
    # Utility.plot_data(Q, Sbest_QFZM, "Sbest_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{FZ}_M(Q)$", "y")
    
    # Utility.plot_data(r, FOpt_r, "FOpt_r", r"$r(nm)$", r"$F(r)$", r"$F^{E}(r)$", "y")
    # Utility.plot_data(r, FOpt_rAL, "FOpt_r", r"$r(nm)$", r"$F(r)$", r"$F^{AL}(r)$", "y")
    # Utility.plot_data(r, FOpt_rFZ, "FOpt_r", r"$r(nm)$", r"$F(r)$", r"$F^{FZ}(r)$", "y")
    # Utility.plot_data(r, FOpt_rALM, "FOpt_r", r"$r(nm)$", r"$F(r)$", r"$F^{AL}_M(r)$", "y")
    # Utility.plot_data(r, FOpt_rFZM, "FOpt_r", r"$r(nm)$", r"$F(r)$", r"$F^{FZ}_M(r)$", "y")
    
    # Scorr_Q = MainFunctions.calc_SQCorr(FOpt_r, r, Q, Sinf)
    # Scorr_QAL = MainFunctions.calc_SQCorr(FOpt_rAL, r, Q, 1)
    # Scorr_QFZ = MainFunctions.calc_SQCorr(FOpt_rFZ, r, Q, 1)
    # Scorr_QALM = MainFunctions.calc_SQCorr(FOpt_rALM, r, Q, 1)
    # Scorr_QFZM = MainFunctions.calc_SQCorr(FOpt_rFZM, r, Q, 1)
    
    # Utility.plot_data(Q, Scorr_Q, "Sbest_Qc", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{E}(Q)$", "y")
    # Utility.plot_data(Q, Scorr_QAL, "Sbest_Qc", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{AL}(Q)$", "y")
    # Utility.plot_data(Q, Scorr_QFZ, "Sbest_Qc", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{FZ}(Q)$", "y")
    # Utility.plot_data(Q, Scorr_QALM, "Sbest_Qc", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{AL}_M(Q)$", "y")
    # Utility.plot_data(Q, Scorr_QFZM, "Sbest_Qc", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{FZ}_M(Q)$", "y")
    
    # end loop
    #-----------------------------------------------------
    
    # print("Eggert", scale_factor, best_sf)
    # print("Eggert", density, best_rho0)
    # print("----------------------")
    # print("AL", scale_factorAL, best_sfAL)
    # print("AL", densityAL, best_rho0AL)
    # print("----------------------")
    # print("FZ", scale_factorFZ, best_sfFZ)
    # print("FZ", densityFZ, best_rho0FZ)
    # print("----------------------")
    
    # Utility.plot_data(Q, Sbest_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^E(Q)$", "y")
    # Utility.plot_data(Q, Sbest_QAL, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{AL}(Q)$", "y")
    # Utility.plot_data(Q, Sbest_QFZ, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{FZ}(Q)$", "y")
    
    # Utility.plot_data(r, FOpt_r, "F_r", r"$r(nm)$", r"$F(r)$", r"$F^E(r)$", "y")
    # Utility.plot_data(r, FOpt_rAL, "F_r", r"$r(nm)$", r"$F(r)$", r"$F^{AL}(r)$", "y")
    # Utility.plot_data(r, FOpt_rFZ, "F_r", r"$r(nm)$", r"$F(r)$", r"$F^{FZ}(r)$", "y")
    
    # chi2, S_Q, F_r = KaplowMethod.Kaplow_method(numAtoms, variables, \
        # Q, I_Q, Ibkg_Q, J_Q, fe_Q, Iincoh_Q, Sinf, Ztot, scale_factor, density, Fintra_r, r)
    
    # SOpt_Q0 = MainFunctions.calc_SQCorr(FOpt_r[0], r, Q, Sinf)
    # SOpt_Q1 = MainFunctions.calc_SQCorr(FOpt_r[1], r, Q, Sinf)
    
    # Utility.plot_data(Q, Sbest_Q[0], "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S(Q)0$", "y")
    # Utility.plot_data(Q, Sbest_Q[1], "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S(Q)1$", "y")
    # Utility.plot_data(Q, SOpt_Q0, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S_{Opt}(Q)0$", "y")
    # Utility.plot_data(Q, SOpt_Q1, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S_{Opt}(Q)1$", "y")
    
    # plt.show()
    
    #-----------------------------------------------------
    # test formalism
    
    # Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scale_factor, Ibkg_Q)
    # aff_sq_mean = Formalism.calc_aff_squared_mean(numAtoms, elementList, Q, elementParameters)
    # aff_mean_sq = Formalism.calc_aff_mean_squared(numAtoms, elementList, Q, elementParameters)
    
    # alpha = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, \
        # Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], \
        # fe_Q[Q<=variables.QmaxIntegrate], Ztot, density)
    
    # alphaW = Formalism.calc_alphaW(Q, Isample_Q, Iincoh_Q, density, aff_sq_mean, aff_mean_sq, 0.001)
    # alphaM = Formalism.calc_alphaM(Q, Isample_Q, Iincoh_Q, density, aff_sq_mean, aff_mean_sq)
    
    # print("Eggert ", alpha)
    # print("Waseda ", alphaW)
    # print("Morard ", alphaM)
   
    # Icoh_Q = MainFunctions.calc_Icoh(numAtoms, alpha, Isample_Q, Iincoh_Q)
    # Icoh_W_Q = MainFunctions.calc_Icoh(numAtoms, alphaW, Isample_Q, Iincoh_Q)
    # Icoh_M_Q = MainFunctions.calc_Icoh(numAtoms, alphaM, Isample_Q, Iincoh_Q)
    
    # S_Q = MainFunctions.calc_SQ(numAtoms, Icoh_Q, Ztot, fe_Q, Sinf, Q, variables.minQ, \
        # variables.QmaxIntegrate, variables.maxQ)
    # S_Qsmoothed = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, variables.smooth_factor, \
        # variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    # S_QsmoothedDamp = UtilityAnalysis.calc_SQdamp(S_Qsmoothed, Q, Sinf, \
        # variables.QmaxIntegrate, variables.damping_factor)
    
    # i_Q = MainFunctions.calc_iQ(S_QsmoothedDamp, Sinf)
    # F_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], \
        # i_Q[Q<=variables.QmaxIntegrate])
    
    # SAL_Q = Formalism.calc_SAL_Q(Q, Icoh_W_Q, aff_sq_mean, variables.minQ, variables.QmaxIntegrate, \
        # variables.maxQ)
    # SAL_Qsmoothed = UtilityAnalysis.calc_SQsmoothing(Q, SAL_Q, 1, variables.smooth_factor, \
        # variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    # SAL_QsmoothedDamp = UtilityAnalysis.calc_SQdamp(SAL_Qsmoothed, Q, 1, \
        # variables.QmaxIntegrate, variables.damping_factor)
    
    # iAL_Q = MainFunctions.calc_iQ(SAL_QsmoothedDamp, 1)
    # FAL_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], \
        # iAL_Q[Q<=variables.QmaxIntegrate])
    
    # SFZ_Q = Formalism.calc_SFZ_Q(Q, Icoh_W_Q, aff_sq_mean, aff_mean_sq, variables.minQ, \
        # variables.QmaxIntegrate, variables.maxQ)
    # SFZ_Qsmoothed = UtilityAnalysis.calc_SQsmoothing(Q, SFZ_Q, 1, variables.smooth_factor, \
        # variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    # SFZ_QsmoothedDamp = UtilityAnalysis.calc_SQdamp(SFZ_Qsmoothed, Q, 1, \
        # variables.QmaxIntegrate, variables.damping_factor)
    
    # iFZ_Q = MainFunctions.calc_iQ(SFZ_QsmoothedDamp, 1)
    # FFZ_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], \
        # iFZ_Q[Q<=variables.QmaxIntegrate])
    
    # SAL_M_Q = Formalism.calc_SAL_Q(Q, Icoh_M_Q, aff_sq_mean, variables.minQ, variables.QmaxIntegrate, \
        # variables.maxQ)
    # SAL_M_Qsmoothed = UtilityAnalysis.calc_SQsmoothing(Q, SAL_M_Q, 1, variables.smooth_factor, \
        # variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    # SAL_M_QsmoothedDamp = UtilityAnalysis.calc_SQdamp(SAL_M_Qsmoothed, Q, 1, \
        # variables.QmaxIntegrate, variables.damping_factor)
    
    # iAL_M_Q = MainFunctions.calc_iQ(SAL_M_QsmoothedDamp, 1)
    # FAL_M_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], \
        # iAL_M_Q[Q<=variables.QmaxIntegrate])
    
    # SFZ_M_Q = Formalism.calc_SFZ_Q(Q, Icoh_M_Q, aff_sq_mean, aff_mean_sq, variables.minQ, \
        # variables.QmaxIntegrate, variables.maxQ)
    # SFZ_M_Qsmoothed = UtilityAnalysis.calc_SQsmoothing(Q, SFZ_M_Q, 1, variables.smooth_factor, \
        # variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    # SFZ_M_QsmoothedDamp = UtilityAnalysis.calc_SQdamp(SFZ_M_Qsmoothed, Q, 1, \
        # variables.QmaxIntegrate, variables.damping_factor)
    
    # iFZ_M_Q = MainFunctions.calc_iQ(SFZ_M_QsmoothedDamp, 1)
    # FFZ_M_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], \
        # iFZ_M_Q[Q<=variables.QmaxIntegrate])
    
    # Utility.plot_data(Q, S_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{E}(Q)$", "y")
    # Utility.plot_data(Q, SAL_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{AL}(Q)$", "y")
    # Utility.plot_data(Q, SFZ_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{FZ}(Q)$", "y")
    # Utility.plot_data(Q, SAL_M_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{AL}_M(Q)$", "y")
    # Utility.plot_data(Q, SFZ_M_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{FZ}_M(Q)$", "y")
    
    # Utility.plot_data(Q, S_QsmoothedDamp, "S_QsmoothedDamp", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{E}(Q)$", "y")
    # Utility.plot_data(Q, SAL_QsmoothedDamp, "S_QsmoothedDamp", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{AL}(Q)$", "y")
    # Utility.plot_data(Q, SFZ_QsmoothedDamp, "S_QsmoothedDamp", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{FZ}(Q)$", "y")
    # Utility.plot_data(Q, SAL_M_QsmoothedDamp, "S_QsmoothedDamp", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{AL}_M(Q)$", "y")
    # Utility.plot_data(Q, SFZ_M_QsmoothedDamp, "S_QsmoothedDamp", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{FZ}_M(Q)$", "y")
    
    # Utility.plot_data(r, F_r, "F_r", r"$r(nm)$", r"$F(r)$", r"$F^{E}(r)$", "y")
    # Utility.plot_data(r, FAL_r, "F_r", r"$r(nm)$", r"$F(r)$", r"$F^{AL}(r)$", "y")
    # Utility.plot_data(r, FFZ_r, "F_r", r"$r(nm)$", r"$F(r)$", r"$F^{FZ}(r)$", "y")
    # Utility.plot_data(r, FAL_M_r, "F_r", r"$r(nm)$", r"$F(r)$", r"$F^{AL}_M(r)$", "y")
    # Utility.plot_data(r, FFZ_M_r, "F_r", r"$r(nm)$", r"$F(r)$", r"$F^{FZ}_M(r)$", "y")
    
    # end formalism test
    #-----------------------------------------------------
    
    
    #-----------------------------------------------------
    # Test alpha
    
    # Isample_Q = MainFunctions.calc_IsampleQ(I_Q, variables.sf_value, Ibkg_Q)
    # val_alphaE = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, \
        # Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], \
        # fe_Q[Q<=variables.QmaxIntegrate], Ztot, variables.rho0_value)
    
    # val_alphaE2 = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf_Q[Q<=variables.QmaxIntegrate], \
        # Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], \
        # fe_Q[Q<=variables.QmaxIntegrate], Ztot, variables.rho0_value)
    
    # aff_sq_mean = Formalism.calc_aff_squared_mean(numAtoms, elementList, Q, elementParameters)
    # aff_mean_sq = Formalism.calc_aff_mean_squared(numAtoms, elementList, Q, elementParameters)
    # gamma = np.linspace(0, 0.01, 1000)
    
    # alphaE = np.empty(gamma.size)
    # alphaE.fill(val_alphaE)
    # alphaE2 = np.empty(gamma.size)
    # alphaE2.fill(val_alphaE2)
    
    # alphaW = np.array([])
    # for g_val in gamma:
        # # print(g_val)
        # val_alphaW = Formalism.calc_alphaW(Q, Isample_Q, Iincoh_Q, variables.rho0_value, aff_sq_mean, aff_mean_sq, g_val)
        # alphaW = np.append(alphaW, val_alphaW)
    
    # alphaM = np.empty(gamma.size)
    # val_alphaM = Formalism.calc_alphaM(Q, Isample_Q, Iincoh_Q, variables.rho0_value, aff_sq_mean, aff_mean_sq)
    # alphaM.fill(val_alphaM)
    
    # print("E", val_alphaE)
    # print("E2", val_alphaE2)
    # # print("W", val_alphaW)
    # print("M", val_alphaM)
    
    # Utility.plot_data(gamma, alphaE, "alpha", r"$\gamma$", r"$\alpha$", r"$\alpha^E$", "y")
    # Utility.plot_data(gamma, alphaE2, "alpha", r"$\gamma$", r"$\alpha$", r"$\alpha^E$ $S_\infty(Q)$", "y")
    # # Utility.plot_data(gamma, alphaW, "alpha", r"$\gamma$", r"$\alpha$", r"$\alpha^W$", "y")
    # Utility.plot_data(gamma, alphaM, "alpha", r"$\gamma$", r"$\alpha$", r"$\alpha^M$", "y")
    
    #-----------------------------------------------------
    
    #-----------------------------------------------------
    # Kaplow method
    # Isample_Q = MainFunctions.calc_IsampleQ(I_Q, sf, Ibkg_Q)
    # alpha = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, \
        # Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], \
        # fe_Q[Q<=variables.QmaxIntegrate], Ztot, rho0)
    # Icoh_Q = MainFunctions.calc_Icoh(numAtoms, alpha, Isample_Q, Iincoh_Q)
    
    # S_Q = MainFunctions.calc_SQ(numAtoms, Icoh_Q, Ztot, fe_Q, Sinf, Q, variables.minQ, \
        # variables.QmaxIntegrate, variables.maxQ)
    # S_Qsmoothed = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, variables.smooth_factor, \
        # variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    # S_QsmoothedDamp = UtilityAnalysis.calc_SQdamp(S_Qsmoothed, Q, Sinf, \
        # variables.QmaxIntegrate, variables.damping_factor)
    
    # i_Q = MainFunctions.calc_iQ(S_QsmoothedDamp, Sinf)
    # F_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], \
        # i_Q[Q<=variables.QmaxIntegrate])
    
    # F_rIt, deltaF_rIt = Optimization.calc_optimize_Fr(variables.iteration, F_r, \
        # Fintra_r, rho0, i_Q[Q<=variables.QmaxIntegrate], Q[Q<=variables.QmaxIntegrate], \
        # Sinf, J_Q[Q<=variables.QmaxIntegrate], r, variables.rmin, variables.plot_iter)
    
    # chi2[idx_rho0][idx_sf] = simps(deltaF_rIt[r < variables.rmin]**2, r[r < variables.rmin])
    #-----------------------------------------------------
    
    plt.show()