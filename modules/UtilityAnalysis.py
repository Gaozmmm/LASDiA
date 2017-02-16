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
#import math
#from scipy import fftpack
from scipy import interpolate
from scipy.integrate import simps
#from scipy import signal
from scipy import constants

from modules import MainFunctions
#from modules import Utility
from modules import UtilityAnalysis


def data_interpolation(Q, I_Q, minQ, maxQ, numPoints):
    """Function to interpolate the data between.
    
    Parameters
    ----------
    Q         : numpy array
                momentum transfer (nm^-1)
    I_Q       : numpy array
                scattering intensity
    minQ      : float
                minimum Q value
    maxQ      : float
                maximum Q value
    numPoints : int
                number of point for the interpolation
    
    Returns
    -------
    Q         : numpy array
                interpolated momentum transfer (nm^-1)
    I_Q       : numpy array
                interpolated scattering intensity
    """
    
    Q, I_Q = rebinning(Q, I_Q, minQ, maxQ, numPoints)
    
    I_Q[Q<=minQ] = 0.0
    
    return (Q, I_Q)
    
    
def rebinning(Q, I_Q, minQ, maxQ, numPoints):
    """Function to rebin the data.
    
    Parameters
    ----------
    Q         : numpy array
                momentum transfer (nm^-1)
    I_Q       : numpy array
                scattering intensity
    minQ      : float
                minimum Q value
    maxQ      : float
                maximum Q value
    numPoints : int
                number of point for the interpolation
    
    Returns
    -------
    Q         : numpy array
                rebinned momentum transfer (nm^-1)
    I_Q       : numpy array
                rebinned scattering intensity
    """
    
    if numPoints == 0:
        numPoints = Q.size
    
    interI_Q = interpolate.interp1d(Q, I_Q)
    Q = np.linspace(np.amin(Q), maxQ, numPoints, endpoint=True)
    I_Q = interI_Q(Q)
    
    return (Q, I_Q)


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
        Q, I_Q = rebinning(Q, I_Q, minQ, maxQ, min_len)
        Qbkg, I_Qbkg = rebinning(Qbkg, I_Qbkg, minQ, maxQ, min_len)
    
    return (Q, I_Q, Qbkg, I_Qbkg)


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
    """Function to find the nearest value to a given number in array.
    
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

    # plt.close()

    print(type(points))
    print(points)
    
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


def calc_iintradamp(iintra_Q, Q, QmaxIntegrate, dampingFunction):
    """Function to apply the damping factor to iintra(Q) function.
    
    Parameters
    ----------
    iintra_Q        : numpy array
                      intramolecular contribution of i(Q)
    Q               : numpy array
                      momentum transfer (nm^-1)
    QmaxIntegrate   : float
                     maximum Q value for the integrations
    dampingFunction : numpy array
                      damping function
    
    Returns
    -------
    Qiintra_Q      : numpy array
                     intramolecular contribution of i(Q) multiplied by the
                     damping function
    """

    iintra_Q = dampingFunction*iintra_Q

    return iintra_Q


def calc_SQsmoothing(Q, S_Q, Sinf, smoothingFactor, minQ, QmaxIntegrate, maxQ):
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
                    maximum Q value for the integrations
    maxQ          : float
                    maximum Q value
    NumPoints     : int
                    number of points in the smoothed S(Q)
    
    Returns
    -------
    S_Qsmoothed   : numpy array
                    smoothed S(Q) with NumPoints dimension
    """
    
    # smooth = interpolate.UnivariateSpline(Q[(Q>minQ) & (Q<=QmaxIntegrate)],
        # S_Q[(Q>minQ) & (Q<=QmaxIntegrate)], k=3, s=smooth_factor)
    # S_Qsmoothed = signal.savgol_filter(S_Q, 51, 3)
    # StdDevWave=smoothingFactor * 10**(-6) * Q**3
    
    smooth = interpolate.UnivariateSpline(Q, S_Q, k=3, s=smoothingFactor)
    S_Qsmoothed = smooth(Q)
    
    S_Qsmoothed[Q<minQ] = 0
    S_Qsmoothed[(Q>=QmaxIntegrate)] = Sinf
    
    return (S_Qsmoothed)


def calc_dampingFunction(Q, dampingFactor, QmaxIntegrate, typeFunction):
    """Function to calculate the damping function.

    Parameters
    ----------
    Q               : numpy array 
                      momentum transfer (nm^-1)
    dampingFactor   : float
                      damping factor
    QmaxIntegrate   : float
                      maximum Q value for the integrations
    typeFunction    : string
                      type of function to use for the damping

    Returns
    -------
    dampingFunction : numpy array
                      damping function
    """
    
    if typeFunction == "Exponential":
        exponentFactor = dampingFactor / QmaxIntegrate**2
        dampingFunction = np.exp(-exponentFactor * Q**2)
    elif typeFunction == "Lorch Function":
        delta = np.pi/QmaxIntegrate
        dampingFunction = np.sin(Q*delta)/Q*delta

    return dampingFunction


def calc_SQdamp(S_Q, Sinf, dampingFunction):
    """Function to apply the damping factor to the structure factor.
    
    Parameters
    ----------
    S_Q             : numpy array 
                      structure factor
    Sinf            : float
                      value of S(Q) for Q->inf
    dampingFunction : float
                      damp factor
    
    Returns
    -------
    S_Qdamp         : numpy array
                      damped structure factor
    """

    S_Qdamp = (dampingFunction * (S_Q - Sinf)) + Sinf

    return S_Qdamp


def S_QCalculation(Q, I_Q, Ibkg_Q, scaleFactor, J_Q, Sinf, fe_Q, Ztot, density, Iincoh_Q, 
    minQ, QmaxIntegrate, maxQ, smoothFactor, dampFactor):
    """Function for the S(Q) calculation.
    """
    
    Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scaleFactor, Ibkg_Q)
    alpha = MainFunctions.calc_alpha(J_Q[Q<=QmaxIntegrate], Sinf,
        Q[Q<=QmaxIntegrate], Isample_Q[Q<=QmaxIntegrate],
        fe_Q[Q<=QmaxIntegrate], Ztot, density)
    Icoh_Q = MainFunctions.calc_Icoh(alpha, Isample_Q, Iincoh_Q)
    
    S_Q = MainFunctions.calc_SQ(Icoh_Q, Ztot, fe_Q, Sinf, Q,
        minQ, QmaxIntegrate, maxQ)
    Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf,
        smoothFactor, minQ, QmaxIntegrate, maxQ)
    SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Q, Sinf,
        QmaxIntegrate, dampFactor)
    
    return SsmoothDamp_Q


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


def make_array_loop(varValue, step, numSample):
    """Function to create an array given its middle value and the percentage
    of the extreme.
    
    Parameters
    ----------
    varValue  : float
                variable's value to generate the array
    step      : float
                array step
    numSample : int
                number of sample
    
    Returns
    -------
    varArray  : numpy array
                variable final array
    """
    
    # lowExtreme = varValue-step*11
    # highExtreme = varValue+step*11
    lowExtreme = varValue
    highExtreme = varValue+step*22
    
    varArray = np.linspace(lowExtreme, highExtreme, numSample)
    
    return varArray


def normalize_to_1(var_y):
    """Function to normalize a distribution to unitary area.
    
    Parameters
    ----------
    var_y      : numpy array
                 distribution to normalize
    
    Returns
    -------
    var_y_norm : numpy array
                 normalized distribution
    
    """
    
    area = abs(simps(var_y))
    var_y_norm = var_y / area
    
    return var_y_norm


def conv_gcm3_to_atnm3(val, atomMass):
    """Function to convert the sample density from g/cm3 to atoms/nm3.
    
    Parameters
    ----------
    val      : float
               density value in g/cm3
    atomMass : float
               atomic mass in g/mol
    
    Returns
    -------
    val_conv : float
               density value in atoms/nm3
    
    """
    
    val_conv = val/atomMass*constants.N_A/10**21
    
    return val_conv


def conv_atnm3_to_gcm3(val, atomMass):
    """Function to convert the sample density from atoms/nm3 to g/cm3.
    
    Parameters
    ----------
    val      : float
               density value in g/cm3
    atomMass : float
               atomic mass in g/mol
    
    Returns
    -------
    val_conv : float
               density value in atoms/nm3
    
    """
    
    val_conv = val*atomMass/constants.N_A*10**21
    
    return val_conv
