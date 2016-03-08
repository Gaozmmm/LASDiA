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

"""Module containing useful functions used in LASDiA.

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

def read_file(path):
    """Read the file and return x and y as numpy vectors
    
    arguments:
    path: file path to read - string
    
    returns:
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
        # I move the factor scale from nm to A into the other functions
        # if dimension.lower() in header2.lower():
            # scale_factor = 10
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
#    xVect /= scale_factor
    yVect = np.array(y)
    
#    plot_data(xVect, yVect,  "Q", "I(Q)", "Data")
    
    return (xVect, yVect)
    
    
def write_file(path, Q, I_Q):
    """Function to write into file
    """
    
    file = open(path, "w")
    file_name = file.name
    
    output = np.column_stack((Q.flatten(),I_Q.flatten()))
    np.savetxt(file_name,output,delimiter='\t')
    file.close()
    
    
def calc_indices(Q, minQ, QmaxIntegrate, maxQ):
    """Function to calculate the ranges where S(Q) is constant
    
    arguments:
    Q: momentum transfer - array
    minQ: minimum Q value - number
    maxQ: maximum Q value - number
    QmaxIntegrate: maximum Q value for the intagrations - number
    
    returns:
    min_index: array index of element with Q<=minQ - array
    max_index: array index of element with Q>QmaxIntegrate & Q<=maxQ - array
    """
    
    min_index = np.where(Q<=minQ)
    max_index = np.where((Q>QmaxIntegrate) & (Q<=maxQ))
    
    return (min_index, max_index)
    
    
def calc_ranges(Q, minQ, QmaxIntegrate, maxQ):
    """Function to calculate the ranges used in the program
    
    arguments:
    Q: momentum transfer - array
    minQ: minimum Q value - number
    maxQ: maximum Q value - number
    QmaxIntegrate: maximum Q value for the intagrations - number
    
    returns:
    validation_index: range of valide Q
    integration_index: range of integration
    calculation_index: range where S(Q) is calculated
    """
    
    validation_index = np.where(Q<=maxQ)
    integration_index = np.where(Q<=QmaxIntegrate)
    calculation_index = np.where((Q>minQ) & (Q<=QmaxIntegrate))
    
    return(validation_index, integration_index, calculation_index)
    
    
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
    
    
def SQsmoothing(Q, S_Q, Sinf, smoothfactor, min_index, minQ, QmaxIntegrate, maxQ, NumPoints):
    """Function for smoothing S(Q).
    This function smooths S(Q) and resets the number of points for the variable Q
    
    arguments:
    Q: momentum transfer - array
    S_Q: structure factor - array
    Sinf: Sinf - number
    smoothfactor: smoothing factor - number
    min_index: array index of element with Q<=minQ - array
    minQ: minimum Q value - number
    maxQ: maximum Q value - number
    QmaxIntegrate: maximum Q value for the intagrations - number
    NumPoints: number of points in the smoothed S(Q) - number
    
    returns:
    newQ: new set of Q with NumPoints dimension - array
    S_Qsmoothed: smoothed S(Q) with NumPoints dimension - array
    """
    
    mask_smooth = np.where((Q>minQ) & (Q<=maxQ))
    smooth = interpolate.UnivariateSpline(Q[mask_smooth], S_Q[mask_smooth], k=3, s=smoothfactor)
    newQ = np.linspace(np.amin(Q), maxQ, NumPoints, endpoint=True)
    S_Qsmoothed = smooth(newQ)
    
    # mask_low = np.where(Q<=minQ)
    num_low = S_Qsmoothed[newQ<minQ].size
    smooth = interpolate.UnivariateSpline(Q[min_index], S_Q[min_index], k=3, s=smoothfactor)
    newQLow = np.linspace(np.amin(newQ), minQ, num_low, endpoint=True)
    S_QsmoothLow = smooth(newQLow)
    
    S_Qsmoothed[newQ<minQ] = S_QsmoothLow
    S_Qsmoothed[(newQ>QmaxIntegrate) & (newQ<=maxQ)] = Sinf
    
    return (newQ, S_Qsmoothed)
    
    
def fitline(X, f_X, index1, element1, index2, element2):
    """Function to flat the peak
    """
    
    xpoints = [element1, element2]
    ypoints = [f_X[index1], f_X[index2]]

    coefficients = np.polyfit(xpoints, ypoints, 1)
    polynomial = np.poly1d(coefficients)
    y_axis = polynomial(X)

    return y_axis
    
    
def find_nearest(array, value):
    """Function to find the nearest array element of a given value
    
    arguments:
    array: array of which it wants to find the nearest element - array
    value: comparing element - number
    
    returns:
    index: index of the nearest element - number
    element: nearest element - number
    """
    
    index = (np.abs(array-value)).argmin()
    element = array[index]
    
    return (index, element)
    
    
def removePeaks(Q, I_Q):
    """Function to remove Bragg's peaks
    
    arguments:
    Q: momentum transfer - array
    I_Q: measured scattering intensity - array
    
    returns:
    I_Q: measured scattering intensity without peaks - array
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