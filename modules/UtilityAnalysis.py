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

"""Module containing useful functions for the analysis used in LASDiA.

The nomenclature and the procedure follow the article:
Eggert et al. 2002 PRB, 65, 174105.

For the functions arguments and the returns I followed this convetion for the notes:
arguments: description - type
returns: description - type.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""

import matplotlib.pyplot as plt

import sys
import os

import numpy as np
import scipy.constants as sc
from scipy import fftpack
from scipy.integrate import simps
from scipy import interpolate
from scipy import signal
import math
import random

from modules.Utility import *

def calc_indices(Q, minQ, QmaxIntegrate, maxQ):
    """Function to calculate the Q ranges where S(Q) is constant.
    
    Parameters
    ----------
    Q             : numpy array 
                    momentum transfer (nm^-1)
    minQ          : float
                    minimum Q value
    maxQ          : float
                    maximum Q value
    QmaxIntegrate : float
                    maximum Q value for the intagrations
    
    Returns
    -------
    min_index     : numpy array
                    indices of elements with Q<=minQ
    max_index     : numpy array
                    indices of elements with Q>QmaxIntegrate & Q<=maxQ
    """

    min_index = np.where(Q<=minQ)
    max_index = np.where((Q>QmaxIntegrate) & (Q<=maxQ))

    return (min_index, max_index)


def calc_ranges(Q, minQ, QmaxIntegrate, maxQ):
    """Function to calculate the Q ranges used in the program.
    
    Parameters
    ----------
    Q                 : numpy array 
                        momentum transfer (nm^-1)
    minQ              : float
                        minimum Q value
    maxQ              : float
                        maximum Q value
    QmaxIntegrate     : float
                        maximum Q value for the intagrations
    
    Returns
    -------
    validation_index  : numpy array
                        range of valide Q
    integration_index : numpy array
                        range where the integration is calculated
    calculation_index : numpy array
                        range where S(Q) is calculated
    """

    validation_index = np.where(Q<=maxQ)
    integration_index = np.where(Q<=QmaxIntegrate)
    calculation_index = np.where((Q>minQ) & (Q<=QmaxIntegrate))

    return(validation_index, integration_index, calculation_index)


def calc_SQsmoothing(Q, S_Q, Sinf, smooth_factor, min_index, minQ, QmaxIntegrate, maxQ, NumPoints):
    """Function for smoothing S(Q).
    This function smooths S(Q) and resets the number of points for the variable Q.
    
    Parameters
    ----------
    Q             : numpy array 
                    momentum transfer (nm^-1)
    S_Q           : numpy array 
                    structure factor
    Sinf          : float
                    Sinf
    smooth_factor : float
                    smoothing factor
    min_index     : numpy array 
                    indices of elements with Q<=minQ
    minQ          : float
                    minimum Q value
    maxQ          : float
                    maximum Q value
    QmaxIntegrate : float
                    maximum Q value for the intagrations
    NumPoints     : int
                    number of points in the smoothed S(Q)
    
    Returns
    -------
    newQ          : numpy array 
                    new set of Q with NumPoints dimension 
    S_Qsmoothed   : numpy array
                    smoothed S(Q) with NumPoints dimension
    """

    mask_smooth = np.where((Q>minQ) & (Q<=maxQ))
    smooth = interpolate.UnivariateSpline(Q[mask_smooth], S_Q[mask_smooth], k=3, s=smooth_factor)
    newQ = np.linspace(np.amin(Q), maxQ, NumPoints, endpoint=True)
    S_Qsmoothed = smooth(newQ)

    # mask_low = np.where(Q<=minQ)
    num_low = S_Qsmoothed[newQ<minQ].size
    smooth = interpolate.UnivariateSpline(Q[min_index], S_Q[min_index], k=3, s=smooth_factor)
    newQLow = np.linspace(np.amin(newQ), minQ, num_low, endpoint=True)
    S_QsmoothLow = smooth(newQLow)

    S_Qsmoothed[newQ<minQ] = S_QsmoothLow
    S_Qsmoothed[(newQ>QmaxIntegrate) & (newQ<=maxQ)] = Sinf

    return (newQ, S_Qsmoothed)


def fitline(Q, I_Q, index1, element1, index2, element2):
    """Function to calculate the fit line between two points.
    This function is used to flat the peaks in the raw data with a first order polynomial.
    
    Parameters
    ----------
    Q        : numpy array 
               momentum transfer (nm^-1)
    I_Q      : numpy array
               measured scattering intensity
    index1   : int
               first point index
    element1 : float
               first point value
    index2   : int
               second point index
    element2 : float
               second point value
    
    Returns
    -------
    y_axis   : numpy array
               ordinate values of fitted curve
    """
    
    xpoints = [element1, element2]
    ypoints = [I_Q[index1], I_Q[index2]]

    coefficients = np.polyfit(xpoints, ypoints, 1)
    polynomial = np.poly1d(coefficients)
    y_axis = polynomial(Q)

    return y_axis


def find_nearest(array, value):
    """Function to find the nearest array element of a given value.
    
    Parameters
    ----------
    array   : numpy array
              array of which it wants to find the nearest element
    value   : float
              comparing element
    
    Returns
    -------
    index   : int
              index of the nearest element
    element : float
              nearest element
    """

    index = (np.abs(array-value)).argmin()
    element = array[index]

    return (index, element)


def remove_peaks(Q, I_Q):
    """Function to remove Bragg's peaks.
    
    Parameters
    ----------
    Q   : numpy array
          momentum transfer (nm^-1)
    I_Q : numpy array
          measured scattering intensity
    
    Returns
    -------
    I_Q : numpy array 
          measured scattering intensity without peaks
    """

    plt.figure('Remove Peaks')
    plt.plot(Q, I_Q)
    plt.grid()
    plt.xlabel('Q')
    plt.ylabel('I(Q)')

    points = np.array(plt.ginput(n=0, timeout=0, show_clicks=True, mouse_add=1, mouse_pop=3, mouse_stop=2))

    plt.close()

    idx_elem = np.zeros(shape=(len(points),2))

    for i in range(0, len(points)):
        idx_elem[i] = find_nearest(Q, points[i,0])

    zipped_idx = np.array(list(zip(*[iter(idx_elem[:,0])]*2)))
    zipped_elem = np.array(list(zip(*[iter(idx_elem[:,1])]*2)))

    for i in range(0, len(zipped_elem)):
        mask = np.where((Q>=zipped_elem[i,0]) & (Q<=zipped_elem[i,1]))
        I_Q1 = fitline(Q[mask], I_Q, zipped_idx[i,0], zipped_elem[i,0], zipped_idx[i,1], zipped_elem[i,1])
        I_Q[mask] = I_Q1

    return I_Q


def calc_iintradamp(iintra_Q, Q, QmaxIntegrate, damping_factor):
    """Function to apply the damping factor to iintra(Q) function.
    
    Parameters
    ----------
    iintra_Q       : numpy array
                     intramolecular contribution of i(Q)
    Q              : numpy array
                     momentum transfer (nm^-1)
    QmaxIntegrate  : float
                     maximum Q value for the intagrations
    damping_factor : float
                     damping factor
    
    Returns
    -------
    Qiintra_Q      : numpy array
                     intramolecular contribution of i(Q) multiplied by Q and the damping function
    """

    # damping_factor = 0.5 # np.log(10)
    exponent_factor = damping_factor / QmaxIntegrate**2
    damp_Q = np.exp(-exponent_factor * Q**2)

    Qiintra_Q = Q*iintra_Q
    Qiintra_Q = damp_Q*Qiintra_Q

    return Qiintra_Q


def calc_SQdamp(S_Q, Q, Sinf, QmaxIntegrate, damping_factor):
    """Function to apply the damping factor to the structure factor.
    
    Parameters
    ----------
    S_Q            : numpy array 
                     structure factor
    Q              : numpy array
                     momentum transfer (nm^-1)
    Sinf           : float
                     Sinf
    QmaxIntegrate  : float
                     maximum Q value for the intagrations
    damping_factor : float
                     damping factor
    
    Returns
    -------
    S_Qdamp        : numpy array
                     damped structure factor
    """

    exponent_factor = damping_factor / QmaxIntegrate**2
    damp_Q = np.exp(-exponent_factor * Q**2)

    S_Qdamp = (damp_Q * (S_Q - Sinf)) + Sinf

    return S_Qdamp


def Qto2theta(Q):
    """Function to convert Q into 2theta.
    
    Parameters
    ----------
    Q       : numpy array
              momentum transfer (nm^-1)
    
    
    Returns
    -------
    _2theta : numpy array
              2theta angle (rad)
    """

    wavelenght = 0.03738
    theta = np.arcsin((wavelenght*Q) / (4*np.pi))
    _2theta = 2*theta

    # return np.degrees(theta2)
    return _2theta


def check_data_length(Q, I_Q, Qbkg, I_Qbkg, minQ, maxQ):
    """Function to check if the raw data of measured and the background have the same numeber of points.
    If the number is different the function rebin them.
    
    Parameters
    ----------
    Q      : numpy array
             momentum transfer (nm^-1)
    I_Q    : numpy array
             measured scattering intensity
    Qbkg   : numpy array
             background momentum transfer (nm^-1)
    I_Qbkg : numpy array
             background scattering intensity
    minQ   : float
             minimum Q value
    maxQ   : float
             maximum Q value
    
    Returns
    -------
    Q      : numpy array
             rebinned momentum transfer (nm^-1)
    I_Q    : numpy array
             rebinned measured scattering intensity
    Qbkg   : numpy array
             rebinned background momentum transfer (nm^-1)
    I_Qbkg : numpy array
             rebinned background scattering intensity
    """
    
    if len(Q) != len(Qbkg):
        min_len = len(Q) if len(Q)<len(Qbkg) else len(Qbkg)
        Q, I_Q = rebinning(Q, I_Q, 1, min_len, maxQ, minQ)
        Qbkg, I_Qbkg = rebinning(Qbkg, I_Qbkg, 1, min_len, maxQ, minQ)
    
    return (Q, I_Q, Qbkg, I_Qbkg)