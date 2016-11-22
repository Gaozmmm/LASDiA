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
from scipy import interpolate

from modules import Formalism
from modules import MainFunctions
from modules import Optimization
from modules import Utility
from modules import UtilityAnalysis


def Kaplow_method(variables, Q, I_Q, Ibkg_Q, J_Q, fe_Q, Iincoh_Q, \
    Sinf, Ztot, sf, rho0, Fintra_r, r):
    """Function to apply the Kaplow method.

    Parameters
    ----------
    variables     : module
                    input variables setted by the user
    Q             : numpy array
                    momentum transfer (nm^-1)
    I_Q           : numpy array
                    measured scattering intensity
    Ibkg_Q        : numpy array
                    background scattering intensity
    J_Q           : numpy array
                    J(Q)
    fe_Q          : numpy array
                    effective electric form factor
    Iincoh_Q      : numpy array
                    incoherent scattering intensity
    Sinf          : float
                    value of S(Q) for Q->inf
    Ztot          : int
                    total Z number
    sf            : float
                    scale factor
    rho0          : float
                    average atomic density
    Fintra_r      : numpy array
                    intramolecular contribution of F(r)
    r             : numpy array
                    atomic distance (nm)

    Returns
    -------
    chi2          : float
                    chi2 value
    SsmoothDamp_Q : numpy array
                    smoothed and damped structure factor
    Fopt_r        : numpy array
                    optimized F(r)
    """

    Isample_Q = MainFunctions.calc_IsampleQ(I_Q, sf, Ibkg_Q)
    alpha = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, \
        Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], \
        fe_Q[Q<=variables.QmaxIntegrate], Ztot, rho0)
    Icoh_Q = MainFunctions.calc_Icoh(alpha, Isample_Q, Iincoh_Q)

    S_Q = MainFunctions.calc_SQ(Icoh_Q, Ztot, fe_Q, Sinf, Q, variables.minQ, \
        variables.QmaxIntegrate, variables.maxQ)
    Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, variables.smoothFactor, \
        variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Q, Sinf, \
        variables.QmaxIntegrate, variables.dampFactor)

    i_Q = MainFunctions.calc_iQ(SsmoothDamp_Q, Sinf)
    F_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], \
        i_Q[Q<=variables.QmaxIntegrate])

    Fopt_r, deltaF_rIt = Optimization.calc_optimize_Fr(variables.iteration, F_r, \
        Fintra_r, rho0, i_Q[Q<=variables.QmaxIntegrate], Q[Q<=variables.QmaxIntegrate], \
        Sinf, J_Q[Q<=variables.QmaxIntegrate], r, variables.rmin, variables.plot_iter)

    chi2 = simps(deltaF_rIt[r < variables.rmin]**2, r[r < variables.rmin])

    return (chi2, SsmoothDamp_Q, F_r, Fopt_r)


def Kaplow_methodFZ(numAtoms, variables, Q, I_Q, Ibkg_Q, aff_squared_mean, \
    aff_mean_squared, Iincoh_Q, sf, rho0, Fintra_r, r):
    """Function to apply the Kaplow method with Waseda FZ formalism.

    Parameters
    ----------
    numAtoms  : int
                number of atoms in the molecule
    variables : module
                input variables setted by the user
    Q         : numpy array
                momentum transfer (nm^-1)
    I_Q       : numpy array
                measured scattering intensity
    Ibkg_Q    : numpy array
                background scattering intensity
    aff_squared_mean  : numpy array
                        mean of the squared form factor: <f^2>
    aff_mean_squared : numpy array
                       squared of the mean form factor: <f>^2
    Iincoh_Q : numpy array
               incoherent scattering intensity
    sf       : float
               scale factor
    rho0     : float
               average atomic density
    Fintra_r : numpy array
               intramolecular contribution of F(r)
    r        : numpy array
               atomic distance (nm)

    Returns
    -------
    chi2     : float
               chi2 value
    """

    Isample_Q = MainFunctions.calc_IsampleQ(I_Q, sf, Ibkg_Q)
    alpha = Formalism.calc_alphaFZ(Q, Isample_Q, Iincoh_Q, rho0, aff_squared_mean, aff_mean_squared)
    Icoh_Q = MainFunctions.calc_Icoh(alpha, Isample_Q, Iincoh_Q)

    S_Q = Formalism.calc_SFZ_Q(Q, Icoh_Q, aff_squared_mean, aff_mean_squared, variables.minQ, \
        variables.QmaxIntegrate, variables.maxQ)
    S_Qsmoothed = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, 1, variables.smooth_factor, \
        variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    S_QsmoothedDamp = UtilityAnalysis.calc_SQdamp(S_Qsmoothed, Q, 1, \
        variables.QmaxIntegrate, variables.damping_factor)

    i_Q = MainFunctions.calc_iQ(S_QsmoothedDamp, 1)
    F_r = MainFunctions.calc_Fr(r, Q[Q<=variables.QmaxIntegrate], \
        i_Q[Q<=variables.QmaxIntegrate])

    # J_Q = Iincoh_Q/aff_mean_squared

    F_rIt, deltaF_rIt = Formalism.calc_optimize_FrFZ(variables.iteration, F_r, \
        Fintra_r, rho0, i_Q[Q<=variables.QmaxIntegrate], Q[Q<=variables.QmaxIntegrate], \
        Iincoh_Q[Q<=variables.QmaxIntegrate], r, variables.rmin)

    chi2 = simps(deltaF_rIt[r < variables.rmin]**2, r[r < variables.rmin])

    return (chi2, S_QsmoothedDamp, F_rIt)
