# The MIT License (MIT)

# Copyright (c) 2015 Francesco Devoto

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

"""Set of modules used in the Liquid program

The nomenclature and the procedures follow the article:
Eggert et al. 2002 PRB, 65, 174105

If the fuctions have arguments I followed this convetion for the notes:
argument: description - type
"""

import matplotlib.pyplot as plt

import sys
import os

import numpy as np
import scipy.constants as sc
from scipy import fftpack
from scipy.integrate import simps

def read_file(path):
    """Read the file and return x and y as numpy vectors
    
    arguments:
    path: file path to read - string
    
    return:
    xVect: abscissa values - array
    yVect: ordinate values - array
    """
    
    file = open(path, "r")
    file_name = file.name
    ext = os.path.splitext(file_name)[1]
    
    dimension = "nano"
    scale_factor = 1
    
    if ext == ".chi":
        # Read and ignore the first 4 lines:
        header1 = file.readline()
        header2 = file.readline()
        header3 = file.readline()
        header4 = file.readline()
        if dimension.lower() in header2.lower():
            scale_factor = 10
    elif ext == ".xy":
        # Read and ignore the first 17 lines:
        header1 = file.readline()
        header2 = file.readline()
        header3 = file.readline()
        header4 = file.readline()
        header5 = file.readline()
        header6 = file.readline()
        header7 = file.readline()
        header8 = file.readline()
        header9 = file.readline()
        header10 = file.readline()
        header11 = file.readline()
        header12 = file.readline()
        header13 = file.readline()
        header14 = file.readline()
        header15 = file.readline()
        header16 = file.readline()
        header17 = file.readline()


    # Read all the file in one time and close the file:
    lines = file.readlines()
    file.close()

    # Set the variables:
    # in our case the variables are: Q and I(Q)
    x = []
    y = []

    # Scan the rows of the file stored in lines, and put the values into the variables:
    for line in lines:
        columns = line.split()
        x.append(float(columns[0]))
        y.append(float(columns[1]))

    # Modify the variables as numpy array:
    xVect = np.array(x)
    xVect /= scale_factor
    yVect = np.array(y)
    
#    plot_data(xVect, yVect,  "Q", "I(Q)", "Data")
    
    return (xVect, yVect)

    
# def plot_data(nFigure, xSample, ySample, xLabel, yLabel, style, dataLabel, overlap):
    # """Plot the data read with the function read_file
    #
    # xSample: x sample - numpy array
    # ySample: y sample - numpy array
    # xLabel: x axis label - string
    # yLabel: y axis label - string
    # style: set the line style - string
    # dataLabel: name of the plot in the legend - string
    # overlap: set if to overlap on an existing plot - boolean (to implement)
    # """
    # plt.figure(nFigure) #.add_subplot(311)
    # if plt.fignum_exists(nFigure):
        # plt.figure(nFigure).hold(True)
    
    # plt.plot(xSample, ySample, style, label=dataLabel)
    # plt.legend()
    # plt.xlabel(xLabel)
    # plt.ylabel(yLabel)
    # plt.show()

    
def calc_aff(element, Q):
    """Function to calculate the Atomic Form Factor
    The atomic form factor is calculated with the formula from:
    http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
    
    and the parameters from the article:
    Hajdu Acta Cryst. (1972). A28, 250
    
    arguments:
    element: chemical element - string
    Q: transferse moment - array
    
    return:
    f_Q: atomic form factor - array
    """
    
    # open, read and close the parameters file
    # the first line contain the header and it is useless for the calculation
    file = open("./affParamCEA.txt", "r")
    header1 = file.readline()
    lines = file.readlines()
    file.close()

    # scan the lines and when it find the right element, save the parameters in variables
    for line in lines:
        columns = line.split()
        if columns[0] == element:
            a1 = float(columns[1])
            b1 = float(columns[2])
            a2 = float(columns[3])
            b2 = float(columns[4])
            a3 = float(columns[5])
            b3 = float(columns[6])
            a4 = float(columns[7])
            b4 = float(columns[8])
            c = float(columns[9])
            break
    
    # Calculate the atomic form factor as:
    # f(Q) = f1(Q) + f2(Q) + f3(Q) + f4(Q) + c
    # fi(Q) = ai * exp(-bi * (Q/4pi)^2)
    
    f1_Q = a1 * np.exp(-b1 * (Q/(4*np.pi))**2)
    f2_Q = a2 * np.exp(-b2 * (Q/(4*np.pi))**2)
    f3_Q = a3 * np.exp(-b3 * (Q/(4*np.pi))**2)
    f4_Q = a4 * np.exp(-b4 * (Q/(4*np.pi))**2)
    
    f_Q = f1_Q + f2_Q + f3_Q + f4_Q + c
    
    return f_Q
      
    
def calc_eeff(elementList, Q, calc_aff):
    """Function to calculate the effective electron Form Factor, fe
    
    arguments:
    elementList: contains the chemical elements of the sample with their multiplicity - dictionary("element": multiplicity)
        element - string; multiplicity - number
    Q: transferse moment - numpy array
    calc_aff: function to calculate the atomic form factor - function
    
    return:
    fe_Q: effective electric form factor - array
    Ztot: total Z number - number
    """
    
    fe_Q = 0
    Ztot = 0

    file = open("./incohParamCEA.txt", "r")
    header1 = file.readline()
    lines = file.readlines()
    file.close()

    # scan the lines and when it find the right element take the Z number
    for element, multiplicity in elementList.items():
        # print (element, multiplicity)
        for line in lines:
            columns = line.split()
            if columns[0] == element:
                Ztot += multiplicity * float(columns[1])
                break
    
    # print (Ztot)
    
    for element, multiplicity in elementList.items():
        fe_Q += multiplicity * calc_aff(element, Q)
    
    fe_Q /= Ztot

    return (fe_Q, Ztot)
    
    
def calc_Iincoh(elementList, Q, calc_aff):
    """ Function to calculate the incoherent scattering intensity Iincoh(Q)
    The incoherent scattering intensity is calculated with the formula and parameters from the article:
    Hajdu Acta Cryst. (1972). A28, 250
    
    arguments:
    elementList: contains the chemical elements of the sample with their multiplicity - dictionary("element": multiplicity)
        element - string; multiplicity - number
    Q: transferse moment - numpy array
    calc_aff: function to calculate the atomic form factor - function
    
    return:
    Iincoh_Q: incoherent scattering intensity - array
    """
    
    file = open("./incohParamCEA.txt", "r")
    header1 = file.readline()
    lines = file.readlines()
    file.close()

    Iincoh_Q = 0
    
    # scan the lines and when it find the right element take the Z number
    for element, multiplicity in elementList.items():
        aff = calc_aff(element, Q)
        for line in lines:
            columns = line.split()
            if columns[0] == element:
                Z = float(columns[1])
                M = float(columns[2])
                K = float(columns[3])
                L = float(columns[4])
                break

        Iincoh_Q += multiplicity * ((Z - aff**2/Z ) * (1 - M * (np.exp(-K*Q/(4*np.pi)) - np.exp(-L*Q/(4*np.pi)))))
        
    return Iincoh_Q
    
    
def calc_JQ(Iincoh_Q, Ztot, fe_Q):
    """Function to calculate J(Q)
    
    arguments:
    Iincoh_Q: incoherent scattering intensity - array
    Ztot: total Z number - number
    fe_Q: effective electric form factor - array
    
    return:
    J_Q: J(Q) - array
    """
    
    # J_Q = Iincoh_Q/(Ztot**2 * fe_Q**2)
    J_Q = Iincoh_Q/(fe_Q**2)
    
    return J_Q


def calc_Kp(fe_Q, element, Q, calc_aff):
    """Function to calculate Kp

    """

    Kp = np.mean(calc_aff(element,Q)/fe_Q)

    return Kp


def calc_Sinf(elementList, fe_Q, Q, Ztot, calc_Kp, calc_aff):
    """Function to calculate Sinf
    
    arguments:
    Kp
    
    return:    
    """
    
    sum_Kp2 = 0

    for element, multiplicity in elementList.items():
        sum_Kp2 += multiplicity * calc_Kp(fe_Q, element, Q, calc_aff)**2
    
    Sinf = sum_Kp2 / Ztot**2

    return Sinf


def calc_alpha(J_Q, Sinf, Q, I_Q, fe_Q, Ztot, rho0):
    """Function to calculate alpha

    arguments:

    return:
    """
    
    Integral1 = simps((J_Q + Sinf) * Q**2, Q)
    Integral2 = simps((I_Q/fe_Q**2) * Q**2,Q)

    alpha = Ztot**2 * (((-2*np.pi**2*rho0) + Integral1) / Integral2)

    return alpha

#-------------------------------------------------------------
# From this point onwards I implemented functions as I see in the Igor code
# and the Gunnar email, but I am not understanding how they works and
# their meaning...
    
def calc_SQ(alpha, Isample, fe_Q, J_Q, Ztot, Sinf):
    """Function to calculate the S(Q)

    arguments:

    return:
    """

    S_Q = (alpha*Isample / fe_Q**2 - J_Q)/Ztot**2+Sinf
    
    return S_Q
    
    
def FFT_QiQ(QiQ):
    """ I don't understand very much this part...
    Function to calculate the Fast Fourier Transformation
    
    arguments:
    
    return:
    
    """
    
    val_fft = fftpack.fft(QiQ)
    imag_fft = val_fft.imag
    delta = QiQ[1] - QiQ[0]
    F_r = 2/np.pi * imag_fft * delta
    
    return F_r
    
    
def IFFT_Fr(Fr):
    """ I don't undeerstand very much this part...
    Function to calculate the Inverse Fast Fourier Transformation
    
    arguments:
    
    return:
    """

    val_fft = fftpack.fft(Fr)
    imag_fft = val_fft.imag
    delta = Fr[1] - Fr[0]
    QiQ = imag_fft * delta
    
    return QiQ
        
    
def ItterateDeltaF(Q, SQ_SS, SQCorr):
    """
    
    """
    
    Qmax = np.amax(Q)
    Damping = 0.5
    argQmax = np.aargmax(Q)
    
    QiQ_SS = exp(-(Damping/Qmax**2)*x^2)*Q*(SQ_SS-SQ_SS(argQmax+1))
    QiQCorr = Q*(SQCorr-SQCorr(argQmax+1))
    Del_iQ = (QiQ_SS-QiQCorr)/(Q*SQCorr)
    
    return (QiQ_SS, QiQCorr, Del_iQ)
    
def FitRemoveGofRPeaks():
    

def OptimizeScaleGofRCorr():

def OptimizeDensityGofRCorr():

def OptimizeWsampleGofRCorr():

def OptimizeWrefGofRCorr():
