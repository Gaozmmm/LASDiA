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
    
    # r, Fintra_r = Optimization.calc_intraComponent(Q, fe_Q, Ztot, variables.QmaxIntegrate, \
        # variables.maxQ, elementList, element, x, y, z, elementParameters, variables.damping_factor)
    
    # scale_factor = variables.sf_value
    # density = variables.rho0_value
    # percentage = 20
    
    #-----------------------------------------------------
    # Test alpha
    
    Isample_Q = MainFunctions.calc_IsampleQ(I_Q, variables.sf_value, Ibkg_Q)
    val_alphaE = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, \
        Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], \
        fe_Q[Q<=variables.QmaxIntegrate], Ztot, variables.rho0_value)
    
    val_alphaE2 = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf_Q[Q<=variables.QmaxIntegrate], \
        Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], \
        fe_Q[Q<=variables.QmaxIntegrate], Ztot, variables.rho0_value)
    
    aff_sq_mean = Formalism.calc_aff_squared_mean(numAtoms, elementList, Q, elementParameters)
    aff_mean_sq = Formalism.calc_aff_mean_squared(numAtoms, elementList, Q, elementParameters)
    gamma = np.linspace(0, 0.01, 1000)
    
    alphaE = np.empty(gamma.size)
    alphaE.fill(val_alphaE)
    alphaE2 = np.empty(gamma.size)
    alphaE2.fill(val_alphaE2)
    
    alphaW = np.array([])
    for g_val in gamma:
        # print(g_val)
        val_alphaW = Formalism.calc_alphaW(Q, Isample_Q, Iincoh_Q, variables.rho0_value, aff_sq_mean, aff_mean_sq, g_val)
        alphaW = np.append(alphaW, val_alphaW)
    
    alphaM = np.empty(gamma.size)
    val_alphaM = Formalism.calc_alphaM(Q, Isample_Q, Iincoh_Q, variables.rho0_value, aff_sq_mean, aff_mean_sq)
    alphaM.fill(val_alphaM)
    
    print("E", val_alphaE)
    print("E2", val_alphaE2)
    # print("W", val_alphaW)
    print("M", val_alphaM)
    
    Utility.plot_data(gamma, alphaE, "alpha", r"$\gamma$", r"$\alpha$", r"$\alpha^E$", "y")
    Utility.plot_data(gamma, alphaE2, "alpha", r"$\gamma$", r"$\alpha$", r"$\alpha^E$ $S_\infty(Q)$", "y")
    # Utility.plot_data(gamma, alphaW, "alpha", r"$\gamma$", r"$\alpha$", r"$\alpha^W$", "y")
    Utility.plot_data(gamma, alphaM, "alpha", r"$\gamma$", r"$\alpha$", r"$\alpha^M$", "y")
    
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
        # percentage /= 2
        
        # print(jj, scale_factor, best_sf)
        # print(jj, density, best_rho0)
        # print("----------------------")
        
        # chi2_min = 1000000
        
        # for (idx_rho0, rho0), (idx_sf, sf) in product(enumerate(rho0_array), enumerate(sf_array)):
            # chi2AL[idx_rho0][idx_sf], S_Q, F_r = KaplowMethod.Kaplow_methodWAL(numAtoms, variables, \
                # Q, I_Q, Ibkg_Q, aff_sq_mean, aff_mean_sq, Iincoh_Q, sf, rho0, Fintra_r, r)
            # print(sf, rho0)
            # if chi2AL[idx_rho0][idx_sf] < chi2_min:
                # chi2_min = chi2AL[idx_rho0][idx_sf]
                # best_sfAL = sf
                # best_rho0AL = rho0
                # FOpt_rAL = F_r
                # Sbest_QAL = S_Q
        
        # Utility.plot_data(Q, Sbest_QAL, "S_QAL", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S(Q)$", "y")
        
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
        
        
        # jj += 1
        # if jj == 1:
            # break
        
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
        # fe_Q[Q<=variables.QmaxIntegrate], Ztot, rho0)
    
    # alphaW = Formalism.calc_alphaW(Q, Isample_Q, Iincoh_Q, density, aff_sq_mean, aff_mean_sq)
    # alphaM = Formalism.calc_alphaM(Q, Isample_Q, Iincoh_Q, density, aff_sq_mean, aff_mean_sq)
    
    # print("Waseda ", alphaW)
    # print("Morard ", alphaM)
   
    # Icoh_Q = MainFunctions.calc_Icoh(numAtoms, alpha, Isample_Q, Iincoh_Q)
    # Icoh_W_Q = MainFunctions.calc_Icoh(numAtoms, alphaW, Isample_Q, Iincoh_Q)
    
   
    # SAL_Q = Formalism.calc_SAL_Q(Q, Icoh_W_Q, aff_sq_mean, variables.minQ, variables.QmaxIntegrate, \
        # variables.maxQ)
    # SFZ_Q = Formalism.calc_SFZ_Q(Q, Icoh_W_Q, aff_sq_mean, aff_mean_sq, variables.minQ, \
        # variables.QmaxIntegrate, variables.maxQ)
    
    # S_Q = MainFunctions.calc_SQ(numAtoms, Icoh_Q, Ztot, fe_Q, Sinf, Q, variables.minQ, \
        # variables.QmaxIntegrate, variables.maxQ)
    
    # Utility.plot_data(Q, SAL_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{AL}(Q)$", "y")
    # Utility.plot_data(Q, SFZ_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{FZ}(Q)$", "y")
    # Utility.plot_data(Q, S_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{E}(Q)$", "y")
    
    # end formalism test
    #-----------------------------------------------------
    
    plt.show()