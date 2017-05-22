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

"""Module containing the functions used for the F(r) optimization.

The nomenclature and the procedure follow the article:
Eggert et al. 2002 PRB, 65, 174105.

For the functions arguments and the returns I followed this convetion for the notes:
arguments: description - type
returns: description - type.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by
an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""


import matplotlib.pyplot as plt
import numpy as np
import time

from modules import MainFunctions
from modules import Utility
from modules import UtilityAnalysis


def calc_iintra(Q, fe_Q, Ztot, QmaxIntegrate, maxQ, elementList, element, x, y, z, elementParameters):
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
                Kpi = MainFunctions.calc_Kp(fe_Q, element[ielem], Q, elementParameters)
                Kpj = MainFunctions.calc_Kp(fe_Q, element[jelem], Q, elementParameters)
                KK = Kpi * Kpj
                d = Utility.calc_distMol(x[ielem], y[ielem], z[ielem], x[jelem], y[jelem], z[jelem])
                if d != 0.0:
                    iintra_Q += KK * np.sin(d*Q) / (d*Q)
                    iintra_Q[Q==0.0] = KK

    iintra_Q[(Q>QmaxIntegrate) & (Q<=maxQ)] = 0.0
    iintra_Q /= Ztot**2
    # iintra_Q /= 3

    return iintra_Q


def calc_intraComponent(Q, fe_Q, Ztot, QmaxIntegrate, maxQ, elementList, element, \
    x, y, z, elementParameters, dampingFunction):
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
    dampingFunction   : numpy array
                        damping function

    Returns
    -------
    r                 : numpy array
                        atomic distance (nm)
    Fintra_r          : numpy array
                        intramolecular contribution of F(r)
    """

    iintra_Q = calc_iintra(Q, fe_Q, Ztot, QmaxIntegrate, maxQ, elementList, element, \
        x, y, z, elementParameters)
    iintradamp_Q = UtilityAnalysis.calc_iintradamp(iintra_Q, Q, QmaxIntegrate, dampingFunction)
    r = MainFunctions.calc_r(Q)
    Fintra_r = MainFunctions.calc_Fr(r, Q[Q<=QmaxIntegrate], iintradamp_Q[Q<=QmaxIntegrate])

    return (r, iintradamp_Q, Fintra_r)


def calc_Ftheor(Fintra_r, r, density):
    """Function to calculate the F(r) theoretical behavior.
    """

    Fth_r = Fintra_r - 4*np.pi*r*density

    return Fth_r


def calc_deltaFr(F_r, Fintra_r, r, density):
    """Function to calculate deltaF(r) (eq. 44, 48).

    Parameters
    ----------
    F_r      : numpy array
               F(r)
    Fintra_r : numpy array
               intramolecular contribution of F(r)
    r        : numpy array
               atomic distance (nm)
    density     : float
               atomic density

    Returns
    -------
    deltaF_r : numpy array
               difference between F(r) and its theoretical value
    """

    deltaF_r = F_r - (Fintra_r - 4*np.pi*r*density)

    return deltaF_r


def calc_iQi(i_Q, Q, Sinf, J_Q, deltaF_r, r, rmin):
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

    i_Qi = i_Q - ( 1/Q * ( i_Q / (Sinf + J_Q) + 1)) * integral

    return i_Qi


def calc_optimize_Fr(iterations, F_r, Fintra_r, density, i_Q, Q, Sinf, J_Q, r,
    rmin, plot_iter):
    """Function to calculate the F(r) optimization (eq 47, 48, 49).

    Parameters
    ----------
    iterations : int
                 number of iterations
    F_r        : numpy array
                 F(r)
    density       : float
                 atomic density
    i_Q        : numpy array
                 i(Q)
    Q          : numpy array
                 momentum transfer (nm^-1)
    Sinf       : float
                 value of S(Q) for Q->inf
    J_Q        : numpy array
                 J(Q)
    r          : numpy array
                 atomic distance (nm)
    rmin       : float
                 r cut-off value (nm)
    plot_iter  : string
                 flag to plot the F(r) iterations

    Returns
    -------
    F_r        : numpy array
                 optimazed F(r)
    deltaF_r   : numpy array
                 difference between the last F(r) and its theoretical value
    """

    # if plot_iter.lower() == "y":
        # plt.ion()
        # plt.figure("F_rIt")
        # plt.plot(r, F_r, label="F(r)")
        # plt.xlabel("r (nm)")
        # plt.ylabel("F(r)")
        # plt.legend()
        # plt.grid(True)

    for i in range(iterations):
        deltaF_r = calc_deltaFr(F_r, Fintra_r, r, density)
        i_Q[0] = 0.0
        i_Q[1:] = calc_iQi(i_Q[1:], Q[1:], Sinf, J_Q[1:], deltaF_r, r, rmin)
        
        r, F_r = MainFunctions.calc_Fr(Q, Q*i_Q)
        # if plot_iter.lower() == "y":
            # j = i+1
            # plt.figure("F_rIt")
            # plt.plot(r, F_r, label="%s iteration F(r)" %j)
            # plt.legend()
            # plt.draw()

            # time.sleep(1.0)

    # if plot_iter.lower() == "y":
        # plt.ioff()

    return (F_r, deltaF_r)
