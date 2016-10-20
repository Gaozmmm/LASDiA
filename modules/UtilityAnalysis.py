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

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by
an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""


import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.integrate import simps

from modules import Utility


def check_data_length(Q, I_Q, Qbkg, I_Qbkg, minQ, maxQ):
    """Function to check if the measured and the background raw data have the
    same number of points.
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


def rebinning(X, f_X, BinNum, Num, maxQ, minQ):
    """Function for the rebinning.
    
    Parameters
    ----------
    X      : numpy array
             abscissa to rebin
    f_X    : numpy array
             ordinate to interpolate
    BinNum : int
             number of points to bin together
    Num    : int
             number of points in data interpolation
    maxQ   : float
             maximum Q value
    minQ   : float
             minimum Q value
    
    Returns
    -------
    BinX   : numpy array
             rebinned abscissa
    BinY   : numpy array
             rebinned ordinate
    """

    newf_X = interpolate.interp1d(X, f_X)
    ShiftX = np.linspace(np.amin(X), maxQ, BinNum*Num, endpoint=True)
    ShiftY = newf_X(ShiftX)

    min = (BinNum - 1)/2 * maxQ /(BinNum * Num - 1)
    max = maxQ - (BinNum - 1)/2 * maxQ / (BinNum*Num - 1)
    BinX = np.linspace(min, max, Num, endpoint=True)
    BinY = np.zeros(Num)

    for i in range(BinNum):
        for j in range(0, Num):
            BinY[j] += ShiftY[j*BinNum+i]

    BinY /= BinNum

    mask = np.where(BinX<=minQ)
    BinY[mask] = 0.0

    # lenX = len(X)
    # numX = 2**int(math.log(lenX,2))
    # rebinnedX = np.linspace(np.amin(X), maxQ, numX, endpoint=True)
    # if min < np.amin(X):
        # min = np.amin(X)

    return (BinX, BinY)


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

    points = np.array(plt.ginput(n=0, timeout=0, show_clicks=True, mouse_add=1, \
        mouse_pop=3, mouse_stop=2))

    plt.close()

    idx_elem = np.zeros(shape=(len(points),2))

    for i in range(0, len(points)):
        idx_elem[i] = find_nearest(Q, points[i,0])

    zipped_idx = np.array(list(zip(*[iter(idx_elem[:,0])]*2)))
    zipped_elem = np.array(list(zip(*[iter(idx_elem[:,1])]*2)))

    for i in range(0, len(zipped_elem)):
        mask = np.where((Q>=zipped_elem[i,0]) & (Q<=zipped_elem[i,1]))
        I_Q1 = fitline(Q[mask], I_Q, zipped_idx[i,0], zipped_elem[i,0], zipped_idx[i,1], \
            zipped_elem[i,1])
        I_Q[mask] = I_Q1

    return I_Q


def Qto2theta(Q):
    """Function to convert Q into 2theta.
    
    Parameters
    ----------
    Q         : numpy array
                momentum transfer (nm^-1)
    
    
    Returns
    -------
    two_theta : numpy array
                2theta angle (rad)
    """

    wavelenght = 0.03738
    theta = np.arcsin((wavelenght*Q) / (4*np.pi))
    two_theta = 2*theta

    # return np.degrees(theta2)
    return two_theta


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
                     intramolecular contribution of i(Q) multiplied by Q and the
                     damping function
    """

    # damping_factor = 0.5 # np.log(10)
    exponent_factor = damping_factor / QmaxIntegrate**2
    damp_Q = np.exp(-exponent_factor * Q**2)

    # Qiintra_Q = Q*iintra_Q
    iintra_Q = damp_Q*iintra_Q

    return iintra_Q


def calc_SQsmoothing(Q, S_Q, Sinf, smooth_factor, minQ, QmaxIntegrate, maxQ):
    """Function for smoothing S(Q).
    This function smooths S(Q) and resets the number of points for the variable Q.
    
    Parameters
    ----------
    Q             : numpy array 
                    momentum transfer (nm^-1)
    S_Q           : numpy array 
                    structure factor
    Sinf          : float
                    value of S(Q) for Q->inf
    smooth_factor : float
                    smoothing factor
    minQ          : float
                    minimum Q value
    QmaxIntegrate : float
                    maximum Q value for the intagrations
    maxQ          : float
                    maximum Q value
    NumPoints     : int
                    number of points in the smoothed S(Q)
    
    Returns
    -------
    S_Qsmoothed   : numpy array
                    smoothed S(Q) with NumPoints dimension
    """
    
    smooth = interpolate.UnivariateSpline(Q[(Q>minQ) & (Q<=QmaxIntegrate)], \
        S_Q[(Q>minQ) & (Q<=QmaxIntegrate)], k=3, s=smooth_factor)
    S_Qsmoothed = smooth(Q)
    
    S_Qsmoothed[Q<minQ] = 0
    S_Qsmoothed[(Q>QmaxIntegrate)] = Sinf
    # S_Qsmoothed[(Q>QmaxIntegrate) & (Q<=maxQ)] = Sinf
    
    return S_Qsmoothed


def calc_SQdamp(S_Q, Q, Sinf, QmaxIntegrate, damping_factor):
    """Function to apply the damping factor to the structure factor.
    
    Parameters
    ----------
    S_Q            : numpy array 
                     structure factor
    Q              : numpy array
                     momentum transfer (nm^-1)
    Sinf           : float
                     value of S(Q) for Q->inf
    QmaxIntegrate  : float
                     maximum Q value for the intagrations
    damping_factor : float
                     damp factor
    
    Returns
    -------
    S_Qdamp        : numpy array
                     damped structure factor
    """

    exponent_factor = damping_factor / QmaxIntegrate**2
    damp_Q = np.exp(-exponent_factor * Q**2)

    S_Qdamp = (damp_Q * (S_Q - Sinf)) + Sinf

    return S_Qdamp


def interpolation_after_smoothing(X, newX, f_X):
    """Function to interpolate the data.
    
    Parameters
    ----------
    X         : numpy array
                old abscissa
    interpf_X : numpy array
                new abscissa
    f_X       : numpy array
                old ordinate
    
    Returns
    -------
    newf_X    : numpy array
                new ordinate
    """
    
    interpf_X = interpolate.interp1d(X, f_X)
    newf_X = interpf_X(newX)
    
    return newf_X


def make_array_loop(var_value, percentage):
    """Function to create an array given its middle value and the percentage
    of the extreme.
    
    Parameters
    ----------
    var_value  : float
                 variable's value to generate the array
    percentage : float
                 percentage to calculate the extreme of the interval
    
    Returns
    -------
    var_array  : numpy array
                 final array
    """
    
    low_extreme = var_value - var_value/percentage
    high_extreme = var_value + var_value/percentage
    
    var_array = np.linspace(low_extreme, high_extreme, 10)
    
    return var_array


def make_array(variables, scale_factor, density, percentage):
    """
    """
    
    if variables.sf_loop == "y":
        sf_array = make_array_loop(scale_factor, percentage)
    else:
        sf_array = np.array([scale_factor])
    if variables.rho0_loop == "y":
        rho0_array = make_array_loop(density, percentage)
    else:
        rho0_array = np.array([density])
    
    return (sf_array, rho0_array)


def normalize_to_1(var_y):
    """Function to normalize a variable to unitary area.
    
    Parameters
    ----------
    var_y : numpy array
            variable to normalize
    
    Returns
    -------
    var_y : numpy array
            normalized variable
    
    """
    
    area = abs(simps(var_y))
    var_y_norm = var_y / area
    
    return var_y_norm
    