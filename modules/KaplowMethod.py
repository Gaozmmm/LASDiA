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

"""Module containing the function for the Kaplow method used to calculate the
chi2 minimization.

The nomenclature and the procedure follow the article:
Eggert et al. 2002 PRB, 65, 174105.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by
an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""

from modules import MainFunctions
from modules import Utility
from modules import UtilityAnalysis

def Kaplow_method(numAtoms, variables, Q, I_Q, Ibkg_Q, J_Q, fe_Q, Iincoh_Q, iintra_Q, \
    Sinf, Ztot, s, rho0, Fintra_r2, r2):
    """Function to apply the Kaplow method
    """
    
    Isample_Q = MainFunctions.calc_IsampleQ(I_Q, s, Ibkg_Q)
    alpha = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, \
        Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], \
        fe_Q[Q<=variables.QmaxIntegrate], Ztot, rho0)
    Icoh_Q = MainFunctions.calc_Icoh(numAtoms, alpha, Isample_Q, Iincoh_Q)
    S_Q = MainFunctions.calc_SQ(numAtoms, Icoh_Q, Ztot, fe_Q, Sinf, Q, variables.minQ, \
        variables.QmaxIntegrate, variables.maxQ)
    
    
    newQ, S_Qsmoothed = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, variables.smooth_factor, \
        variables.minQ, variables.QmaxIntegrate, variables.maxQ, 550)
    S_QsmoothedDamp = UtilityAnalysis.calc_SQdamp(S_Qsmoothed, newQ, Sinf, \
        variables.QmaxIntegrate, variables.damp_factor)
    
    # newfe_Q = UtilityAnalysis.interpolation_after_smoothing(Q, newQ, fe_Q)
    
    r = MainFunctions.calc_r(newQ)
    i_Q = MainFunctions.calc_iQ(S_QsmoothedDamp, Sinf)
    F_r = MainFunctions.calc_Fr(r, newQ[newQ<=variables.QmaxIntegrate], \
        i_Q[newQ<=variables.QmaxIntegrate])
    
    iintra_Q = UtilityAnalysis.interpolation_after_smoothing(Q, newQ, iintra_Q)
    # Fintra_r2 = UtilityAnalysis.interpolation_after_smoothing(r2, r, Fintra_r2)
    
    Qiintradamp = UtilityAnalysis.calc_Qiintradamp(iintra_Q, newQ, variables.QmaxIntegrate, \
        variables.damp_factor)
    Fintra_r = MainFunctions.calc_Fr(r, newQ[newQ<=variables.QmaxIntegrate], \
        Qiintradamp[newQ<=variables.QmaxIntegrate])
    
    print(len(Fintra_r))
    print(len(Fintra_r2))
    print(len(r))
    print(len(r2))
    
    
    # Iincoh_QSmooth = calc_Iincoh(elementList, newQ)
    # J_QSmooth = calc_JQ(Iincoh_QSmooth, Ztot, fe_QSmooth)

    # F_rIt = calc_optimize_Fr(variables.iteration, F_r, Fintra_r, rho0, i_Q[integration_indexSmooth], \
        # newQ[integration_indexSmooth], Sinf, J_QSmooth[integration_indexSmooth], r, rmin, "n")
    
    # Utility.plot_data(Q, S_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S(Q)$", "y")
    # Utility.plot_data(newQ, S_QsmoothedDamp, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S(Q)SD$", "y")
    # Utility.plot_data(Q, fe_Q, "fe_Q", r"$Q(nm^{-1})$", r"$fe(Q)$", r"$fe(Q)$", "y")
    # Utility.plot_data(newQ, newfe_Q, "fe_Q", r"$Q(nm^{-1})$", r"$fe(Q)$", r"$newfe(Q)$", "y")
    # Utility.plot_data(Q, iintra_Q, "iintra_Q", r"$Q(nm^{-1})$", r"$fe(Q)$", r"$iintra(Q)$", "y")
    # Utility.plot_data(newQ, iintra_Q, "iintra_Q", r"$Q(nm^{-1})$", r"$iintra(Q)$", r"$iintra(Q)$", "y")
    # Utility.plot_data(newQ, newQ*iintra_Q, "iintra_Q", r"$Q(nm^{-1})$", r"$iintra(Q)$", r"$Qiintra(Q)$", "y")
    # Utility.plot_data(newQ, Qiintradamp, "iintra_Q", r"$Q(nm^{-1})$", r"$iintra(Q)$", r"$Qiintra(Q)D$", "y")
    Utility.plot_data(r, Fintra_r, "Fintra_r", r"$r(nm)$", r"$F_{intra}(r)$", r"$F_{intra}(r)$", "y")
    Utility.plot_data(r2, Fintra_r2, "Fintra_r", r"$r(nm)$", r"$F_{intra}(r)$", r"$F_{intra}2(r)$", "y")
    
    return (S_Q, r, F_r)