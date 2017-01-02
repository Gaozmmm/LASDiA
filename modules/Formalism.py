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

"""Set of modules used in LASDiA to calculate the relations among the different set of Structure Factors.

The formalisms used are Ashcroft-Langreth and Faber-Ziman.

The nomenclature and the procedure follow the book:
Waseda - The Structure of Non-Crystalline Materials.

For the functions arguments and the returns I followed this convetion for the notes:
argument: description - type
return: description - type.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by an underscore: f(x) -> f_x
otherwise it is just the name.
"""


import sys
import os
import math
import numpy as np
import scipy.constants as sc
from scipy import fftpack
from scipy.integrate import simps
from scipy.stats import chisquare
import matplotlib.pyplot as plt

from modules import MainFunctions
from modules import Optimization
from modules import Utility
from modules import UtilityAnalysis


def calc_alphaFZ(Q, Isample_Q, Iincoh_Q, rho0, aff_squared_mean, aff_mean_squared):
    """Function to calcultate alpha with Morard formula.
    Morard et al. E&PSL, Vol 263, Is 1-2, p. 128-139 (eq. 15).

    Parameters
    ----------
    Q                : numpy array
                       momentum transfer (nm^-1)
    Isample_Q        : numpy array
                       sample scattering intensity
    Iincoh_Q         : numpy array
                       incoherent scattering intensity
    rho0             : float
                       average atomic density
    aff_squared_mean : numpy array
                        mean of the squared form factor: <f^2>
    aff_mean_squared : numpy array
                       squared of the mean form factor: <f>^2

    Returns
    -------
    alpha            : float
                       normalization factor
    """

    Integral1 = simps((Iincoh_Q + aff_squared_mean)/aff_mean_squared * Q**2, Q)
    Integral2 = simps((Isample_Q/aff_mean_squared) * Q**2,Q)
    alpha = ((-2*np.pi**2*rho0) + Integral1) / Integral2

    return alpha


def calc_aff_squared_mean(numAtoms, elementList, Q, elementParameters):
    """Function to calculate the squared mean effective electron Form Factor, <f^2> (eq. FZ-4).

    Parameters
    ----------
    numAtoms          : int
                        number of atoms in the molecule
    elementList       : dictionary("element": multiplicity)
                        chemical elements of the sample with their multiplicity
                        element      : string
                                       chemical element
                        multiplicity : int
                                       chemical element multiplicity
    Q                 : numpy array
                        momentum transfer (nm^-1)
    elementParameters : dictionary("element": parameters)
                        chemical elements of the sample with their parameters
                        element    : string
                                     chemical element
                        parameters : list
                                     list of the parameters
                                     (Z, a1, b1, a2, b2, a3, b3, a4, b4, c, M, K, L)

    Returns
    -------
    aff_squared_mean  : numpy array
                        mean of the squared form factor: <f^2>
    """

    aff_squared_mean = np.zeros(Q.size)

    for element, multiplicity in elementList.items():
        aff_squared_mean += multiplicity * MainFunctions.calc_aff(element, Q, elementParameters)**2

    aff_squared_mean /= numAtoms

    return aff_squared_mean


def calc_aff_mean_squared(numAtoms, elementList, Q, elementParameters):
    """Function to calculate the mean squared effective electron Form Factor, <f>^2 (eq. FZ-5).

    Parameters
    ----------
    numAtoms          : int
                        number of atoms in the molecule
    elementList       : dictionary("element": multiplicity)
                        chemical elements of the sample with their multiplicity
                        element      : string
                                       chemical element
                        multiplicity : int
                                       chemical element multiplicity
    Q                 : numpy array
                        momentum transfer (nm^-1)
    elementParameters : dictionary("element": parameters)
                        chemical elements of the sample with their parameters
                        element    : string
                                     chemical element
                        parameters : list
                                     list of the parameters
                                     (Z, a1, b1, a2, b2, a3, b3, a4, b4, c, M, K, L)

    Returns
    -------
    aff_mean_squared : numpy array
                       squared of the mean form factor: <f>^2
    """

    aff_mean_squared = np.zeros(Q.size)

    for element, multiplicity in elementList.items():
        aff_mean_squared += multiplicity * MainFunctions.calc_aff(element, Q, elementParameters)

    aff_mean_squared /= numAtoms

    aff_mean_squared = aff_mean_squared**2

    return aff_mean_squared


def calc_SFZ_Q(Q, Icoh_Q, aff_squared_mean, aff_mean_squared, minQ, QmaxIntegrate, maxQ):
    """Function to calculate S(Q) with Faber-Ziman formalism.

    Parameters
    ----------
    Q                 : numpy array
                        momentum transfer (nm)
    Icoh_Q            : numpy array (nm)
                        coherent scattering intensity


    Returns
    -------
    S_Q               : numpy array
                        structure factor in FZ formalism
    """

    S_Q = np.zeros(Q.size)
    S_Q[(Q>minQ) & (Q<=QmaxIntegrate)] = (Icoh_Q[(Q>minQ) & (Q<=QmaxIntegrate)] - \
        (aff_squared_mean[(Q>minQ) & (Q<=QmaxIntegrate)] - aff_mean_squared[(Q>minQ) & (Q<=QmaxIntegrate)])) / (aff_mean_squared[(Q>minQ) & (Q<=QmaxIntegrate)])
    # S_Q[(Q>QmaxIntegrate) & (Q<=maxQ)] = 1
    S_Q[(Q>QmaxIntegrate)] = 1

    return S_Q


def calc_iintraFZ(Q, QmaxIntegrate, maxQ, elementList, element, x, y, z, elementParameters, aff_mean_squared):
    """Function to calculate the intramolecular contribution of i(Q) (eq. 41).

    Parameters
    ----------
    Q                 : numpy array
                        momentum transfer (nm^-1)
    fe_Q              : numpy array
                        effective electric form factor
    Ztot              : int
                        total Z number
    QmaxIntegrate     : float
                        maximum Q value for the intagrations
    maxQ              : float
                        maximum Q value
    elementList       : dictionary("element": multiplicity)
                        chemical elements of the sample with their multiplicity
                        element      : string
                                       chemical element
                        multiplicity : int
                                       chemical element multiplicity
    element           : string array
                        array with the elements in the xyz_file
    x, y, z           : float array
                        atomic coordinate in the xyz_file (nm)
    elementParameters : dictionary("element": parameters)
                        chemical elements of the sample with their parameters
                        element    : string
                                     chemical element
                        parameters : list
                                     list of the parameters
                                     (Z, a1, b1, a2, b2, a3, b3, a4, b4, c, M, K, L)

    Returns
    -------
    iintra_Q          : numpy array
                        intramolecular contribution of i(Q)
    """

    iintra_Q = np.zeros(Q.size)
    sinpq = np.zeros(Q.size)

    for ielem in range(len(element)):
        for jelem in range(len(element)):
            if ielem != jelem:
                # print(ielem, jelem)
                # print(type(element[ielem]))
                # print(element[ielem], elementList[element[ielem]], element[jelem], elementList[element[jelem]])
                f_Qi = MainFunctions.calc_aff(element[ielem], Q, elementParameters)
                f_Qj = MainFunctions.calc_aff(element[jelem], Q, elementParameters)
                f_i = np.mean(elementList[element[ielem]] * f_Qi / 3)
                f_j = np.mean(elementList[element[jelem]] * f_Qj / 3)
                # f_i = np.mean(f_Qi)
                # f_j = np.mean(f_Qj)
                ff = f_i * f_j
                d = Utility.calc_distMol(x[ielem], y[ielem], z[ielem], x[jelem], y[jelem], z[jelem])
                if d != 0.0:
                    iintra_Q += ff * np.sin(d*Q) / (d*Q)
                    iintra_Q[Q==0.0] = ff

    iintra_Q[(Q>QmaxIntegrate) & (Q<=maxQ)] = 0.0
    iintra_Q /= np.mean(aff_mean_squared)
    # iintra_Q /= 3

    return iintra_Q


def calc_intraComponentFZ(Q, fe_Q, Ztot, QmaxIntegrate, maxQ, elementList, element, \
    x, y, z, elementParameters, damping_factor, aff_mean_squared):
    """Function to calculate the intra-molecular components.

    Parameters
    ----------
    Q                 : numpy array
                        momentum transfer (nm^-1)
    fe_Q              : numpy array
                        effective electric form factor
    Ztot              : int
                        total Z number
    QmaxIntegrate     : float
                        maximum Q value for the intagrations
    maxQ              : float
                        maximum Q value
    elementList       : dictionary("element": multiplicity)
                        chemical elements of the sample with their multiplicity
                        element      : string
                                       chemical element
                        multiplicity : int
                                       chemical element multiplicity
    element           : string array
                        array with the elements in the xyz_file
    x, y, z           : float array
                        atomic coordinate in the xyz_file (nm)
    elementParameters : dictionary("element": parameters)
                        chemical elements of the sample with their parameters
                        element    : string
                                     chemical element
                        parameters : list
                                     list of the parameters
                                     (Z, a1, b1, a2, b2, a3, b3, a4, b4, c, M, K, L)
    damping_factor    : float
                        damp factor

    Returns
    -------
    r                 : numpy array
                        atomic distance (nm)
    Fintra_r          : numpy array
                        intramolecular contribution of F(r)
    """

    iintra_Q = calc_iintraFZ(Q, QmaxIntegrate, maxQ, elementList, element, x, y, z, elementParameters, aff_mean_squared)
    iintradamp_Q = UtilityAnalysis.calc_iintradamp(iintra_Q, Q, QmaxIntegrate, damping_factor)
    r = MainFunctions.calc_r(Q)
    Fintra_r = MainFunctions.calc_Fr(r, Q[Q<=QmaxIntegrate], iintradamp_Q[Q<=QmaxIntegrate])

    return (r, iintradamp_Q, Fintra_r)


def calc_iQiFZ(i_Q, Q, Iincoh_Q, deltaF_r, r, rmin):
    """Function to calculate the i-th iteration of i(Q) (eq. 46, 49).

    Parameters
    ----------
    i_Q      : numpy array
               i(Q)
    Q        : numpy array
               momentum transfer (nm^-1)
    Sinf     : float
               value of S(Q) for Q->inf
    J_Q      : numpy array
               J(Q)
    deltaF_r : numpy array
               difference between F(r) and its theoretical value
    r        : numpy array
               atomic distance (nm)
    rmin     : float
               r cut-off value (nm)

    Returns
    -------
    i_Qi     : numpy array
               i-th iteration of i(Q)
    """

    mask = np.where(r < rmin)
    rInt = r[mask]
    deltaF_rInt = deltaF_r[mask]

    Deltar = np.diff(rInt)
    meanDeltar = np.mean(Deltar)
    Qr = np.outer(Q, rInt)
    sinQr = np.sin(Qr)
    integral = np.sum(deltaF_rInt * sinQr, axis=1) * meanDeltar

    i_Qi = i_Q - ( 1/Q * (i_Q/Iincoh_Q + 1)) * integral
    # i_Qi = i_Q - ( 1/Q * (i_Q + 1)) * integral

    return i_Qi


def calc_optimize_FrFZ(iteration, F_r, Fintra_r, rho0, i_Q, Q, Iincoh_Q, r, rmin):
    """Function to calculate the F(r) optimization (eq 47, 48, 49).

    Parameters
    ----------
    iteration : int
                number of iterations
    F_r       : numpy array
                F(r)
    rho0      : float
                atomic density
    i_Q       : numpy array
                i(Q)
    Q         : numpy array
                momentum transfer (nm^-1)
    Sinf      : float
                value of S(Q) for Q->inf
    J_Q       : numpy array
                J(Q)
    r         : numpy array
                atomic distance (nm)
    rmin      : float
                r cut-off value (nm)
    plot_iter : string
                flag to plot the F(r) iterations

    Returns
    -------
    F_r       : numpy array
                optimazed F(r)
    deltaF_r  : numpy array
                difference between the last F(r) and its theoretical value
    """

    for i in range(iteration):
        deltaF_r = Optimization.calc_deltaFr(F_r, Fintra_r, r, rho0)
        i_Q = calc_iQiFZ(i_Q, Q, Iincoh_Q, deltaF_r, r, rmin)
        i_Q[0] = 0.0
        F_r = MainFunctions.calc_Fr(r, Q, i_Q)

    return (F_r, deltaF_r)
