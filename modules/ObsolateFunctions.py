# The MIT License (MIT)

# Copyright (c) 2016 Francesco Devoto

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

"""Module containing the obsolate functions used in LASDiA.

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
from scipy.interpolate import UnivariateSpline


def calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, min_index, max_index, calculation_index):
    """Function to calculate the structure factor S(Q) (eq. 18)

    arguments:
    N: number of atoms - number
    Icoh: cohrent scattering intensity - array
    Ztot: total Z number - number
    fe_Q: effective electric form factor - array
    Sinf: Sinf - number
    Q: momentum transfer - array
    min_index: array index of element with Q<minQ - array
    max_index: array index of element with Q>QmaxIntegrate & Q<=maxQ - array
    calculation_index: array index of element in the calculation range Q>minQ & Q<=QmaxIntegrate - array
    
    returns:
    S_Q: structure factor - array
    """
    
    S_Q = Icoh_Q[calculation_index] / (N * Ztot**2 * fe_Q[calculation_index]**2)
        
    S_Qmin = np.zeros(Q[min_index].size)
    S_Q = np.concatenate([S_Qmin, S_Q])
    
    S_Qmax = np.zeros(Q[max_index].size)
    S_Qmax.fill(Sinf)
    S_Q = np.concatenate([S_Q, S_Qmax])
    
    return S_Q
    
    
def interpolation(X, f_X, rebinnedX):
    """Function for the interpolation
    """
    
    interpolatedf_X = interpolate.interp1d(X, f_X)
    newf_X = interpolatedf_X(rebinnedX)
    
    return newf_X
    
    
def smoothing(X, f_X, smoothfactor):
    """Function for smoothing
    """
    
    smooth = interpolate.UnivariateSpline(X, f_X, k=3, s=smoothfactor)
    smoothedf_X = smooth(X)
    
    return smoothedf_X
    
    
def SQsmoothing(Q, S_Q, Sinf, smoothfactor, min_index, max_index, validation_index):
    """Function for smoothing S(Q)
    """
    
    S_Qsmooth = smoothing(Q[validation_index], S_Q[validation_index], smoothfactor)
    
    S_Qsmooth[min_index] = 0.0
    S_Qmax = np.zeros(Q[max_index].size)
    S_Qmax.fill(Sinf)
    S_Qsmooth = np.concatenate([S_Qsmooth, S_Qmax])
    # S_Qsmooth[max_index] = Sinf
    
    return S_Qsmooth
    
    
def SQsmoothing2(Q, S_Q, Sinf, smoothfactor, min_index, max_index, calculation_index):
    """Function for smoothing S(Q)
    """
    
    S_Qsmooth = smoothing(Q[calculation_index], S_Q[calculation_index], smoothfactor)
    
    S_Qmin = np.zeros(Q[min_index].size)
    S_Qsmooth = np.concatenate([S_Qmin, S_Qsmooth])
    
    S_Qmax = np.zeros(Q[max_index].size)
    S_Qmax.fill(Sinf)
    S_Qsmooth = np.concatenate([S_Qsmooth, S_Qmax])
    
    return S_Qsmooth
    
    
def fitcurve(X, f_X, mask):
    """Function to flat the peak
    """
    
    xpoints = X[mask]
    ypoints = f_X[mask]

    coefficients = np.polyfit(xpoints, ypoints, 2)
    polynomial = np.poly1d(coefficients)
    y_axis = polynomial(xpoints)

    return y_axis
    
    
def calc_SQdamp(S_Q, Q, Sinf, QmaxIntegrate, damping_factor):
    """
    """
    
    # damping_factor = 0.5
    # damping_factor = np.log(10)
    exponent_factor = damping_factor / QmaxIntegrate**2
    damp_Q = np.exp(-exponent_factor * Q**2)
    
    S_Qdamp = (damp_Q * (S_Q - Sinf)) + Sinf
    
    return S_Qdamp
    
    
def calc_damp(Q, QmaxIntegrate, damping_factor):
    """Function to calculate the damping function
    """
    
    # damping_factor = 0.5 # np.log(10)
    exponent_factor = damping_factor / QmaxIntegrate**2
    damp_Q = np.exp(-exponent_factor * Q**2)
    
    return damp_Q
    
    
def calc_Fr(r, Q, i_Q):
    """Function to calculate F(r) (eq. 20)
    
    arguments:
    r: radius - array
    Q: momentum transfer - array
    i_Q: i(Q) - array
    
    returns:
    F_r: F(r) - array
    """
    
    F_r = (2.0 / np.pi) * simps(Q * i_Q * np.array(np.sin(np.mat(Q).T * np.mat(r))).T, Q)
    
    DeltaQ = np.diff(Q)
    meanDeltaQ = np.mean(DeltaQ)
    rQ = np.outer(r,Q)
    sinrQ = np.sin(rQ)
    F_r2 = (2.0 / np.pi) * np.sum(Q * i_Q * sinrQ, axis=1) * meanDeltaQ
    
    return (F_r, F_r2)
    
    
def calc_NewDimFFT():
    """Function to redefine the array dimension to use the FFT
    """
    
    idx, elem = find_nearest(newQ, QmaxIntegrate)
    newDim = 2*2*2**math.ceil(math.log(5*(idx+1))/math.log(2))
    Qi_Q2 = np.resize(Qi_Q, newDim)
    Qi_Q2[idx:] = 0.0
    
    newQ2 = np.linspace(np.amin(Q), maxQ, newDim, endpoint=True)
    
    DeltaQ = np.diff(newQ)
    deltaR = 2*np.pi/(DeltaQ[0]*newDim)
    
    
def rebinning(X, f_X, BinNum, Num, maxQ, minQ):
    """Function for the rebinning
    """
    
    newf_X = interpolate.interp1d(X, f_X)
    ShitX = np.linspace(np.amin(X), maxQ, BinNum*Num, endpoint=True)
    ShitY = newf_X(ShitX)
    
    min = (BinNum - 1)/2 * maxQ /(BinNum * Num - 1)
    max = maxQ - (BinNum - 1)/2 * maxQ / (BinNum*Num - 1)
    BinX = np.linspace(min, max, Num, endpoint=True)
    BinY = np.zeros(Num)
    
    for i in range(BinNum):
        for j in range(0, Num):
            BinY[j] += ShitY[j*BinNum+i]
    
    BinY /= BinNum
    
    mask = np.where(X<=minQ)
    BinY[mask] = 0.0
    
    # lenX = len(X)
    # numX = 2**int(math.log(lenX,2))
    # rebinnedX = np.linspace(np.amin(X), maxQ, numX, endpoint=True)
    # if min < np.amin(X):
        # min = np.amin(X)
    
    return (BinX, BinY)
    
    
def calc_gr(r, F_r, rho0):
    """Function to calculate g(r)
    
    """
    
    g_r = 1 + F_r / (4 * np.pi * r * rho0)
    
    return g_r
    
    
def calc_alpha2(J_Q, Sinf, Q, Isample_Q, fe_Q, Ztot, rho0, index):
    """Function to calculate the normalization factor alpha (eq. 34)
    
    arguments:
    J_Q: J(Q) - array
    Sinf: Sinf - number
    Q: momentum transfer - array
    Isample_Q: sample scattering intensity - array
    fe_Q: effective electric form factor - array
    Ztot: total Z number - number
    rho0: average atomic density - number
    index: array index of element in the calculation range - array
    
    returns:
    alpha: normalization factor - number
    """
    
    Integral1 = simps((J_Q[index] + Sinf) * Q[index]**2, Q[index])
    Integral2 = simps((Isample_Q[index]/fe_Q[index]**2) * Q[index]**2,Q[index])
    alpha = Ztot**2 * (((-2*np.pi**2*rho0) + Integral1) / Integral2)
    
    # DeltaQ = np.diff(Q)
    # meanDeltaQ = np.mean(DeltaQ)
    # Int1 = np.sum((J_Q[index] + Sinf) * Q[index]**2) * meanDeltaQ
    # Int2 = np.sum( (Isample_Q[index]/fe_Q[index]**2) * Q[index]**2  ) * meanDeltaQ
    # alpha = Ztot**2 * (((-2*np.pi**2*rho0) + Int1) / Int2)

    return alpha