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


from scipy.integrate import simps

from modules import MainFunctions
from modules import Optimization
from modules import Utility
from modules import UtilityAnalysis


def Kaplow_method(numAtoms, variables, Q, I_Q, Ibkg_Q, J_Q, fe_Q, Iincoh_Q, \
    Sinf, Ztot, s, rho0, Fintra_r, r):
    """Function to apply the Kaplow method
    """
    
    Isample_Q = MainFunctions.calc_IsampleQ(I_Q, s, Ibkg_Q)
    alpha = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, \
        Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], \
        fe_Q[Q<=variables.QmaxIntegrate], Ztot, rho0)
    Icoh_Q = MainFunctions.calc_Icoh(numAtoms, alpha, Isample_Q, Iincoh_Q)
    S_Q = MainFunctions.calc_SQ(numAtoms, Icoh_Q, Ztot, fe_Q, Sinf, Q, variables.minQ, \
        variables.QmaxIntegrate, variables.maxQ)
    
    S_Qsmoothed = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, variables.smooth_factor, \
        variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    S_QsmoothedDamp = UtilityAnalysis.calc_SQdamp(S_Qsmoothed, Q, Sinf, \
        variables.QmaxIntegrate, variables.damping_factor)
    
    i_Q = MainFunctions.calc_iQ(S_QsmoothedDamp, Sinf)
    F_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], \
        i_Q[Q<=variables.QmaxIntegrate])
    
    F_rIt, deltaF_rIt = Optimization.calc_optimize_Fr(variables.iteration, F_r, \
        Fintra_r, rho0, i_Q[Q<=variables.QmaxIntegrate], Q[Q<=variables.QmaxIntegrate], \
        Sinf, J_Q[Q<=variables.QmaxIntegrate], r, variables.rmin, "n")
    
    chi2 = simps(deltaF_rIt[r < variables.rmin]**2, r[r < variables.rmin])
    
    return (chi2)