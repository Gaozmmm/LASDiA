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

For the functions arguments and the returns I followed this convetion for the notes:
arguments: description - type
returns: description - type
"""

import matplotlib.pyplot as plt

import sys
import os

import numpy as np
import scipy.constants as sc
from scipy import fftpack
from scipy.integrate import simps
import math
from scipy import interpolate
from scipy import signal

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
    """Function to write to file
    """
    
    file = open(path, "w")
    file_name = file.name
    
    output = np.column_stack((Q.flatten(),I_Q.flatten()))
    np.savetxt(file_name,output,delimiter='\t')
    file.close()
    
    
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
    
    
# def interpolation(X, f_X, rebinnedX):
    # """Function for the interpolation
    # """
    
    # interpolatedf_X = interpolate.interp1d(X, f_X)
    # newf_X = interpolatedf_X(rebinnedX)
    
    # return newf_X
    
    
def smoothing(X, f_X, smoothfactor):
    """Function for smoothing
    """
    
    sm = interpolate.UnivariateSpline(X, f_X, k=3, s=smoothfactor)
    smoothedf_X = sm(X)
    
    return smoothedf_X
    
    
def SQsmoothing(Q, S_Q, Sinf, smoothfactor, min_index, max_index, validation_index):
    """Function for smoothing S(Q)
    """
    
    sm = interpolate.UnivariateSpline(Q[validation_index], S_Q[validation_index], k=3, s=smoothfactor)
    S_Qsmooth = sm(Q[validation_index])
    
    S_Qsmooth[min_index].fill(0.0)
    S_Qsmooth[max_index].fill(Sinf)
    # S_Qmin = np.zeros(Q[min_index].size)
    # S_Qsmooth = np.concatenate([S_Qmin, S_Qsmooth])
    
    # S_Qmax = np.zeros(Q[max_index].size)
    # S_Qmax.fill(Sinf)
    # S_Qsmooth = np.concatenate([S_Qsmooth, S_Qmax])
    
    # S_Q[index] = S_Qsmooth
    
    return S_Qsmooth
    
    
def fitline(X, f_X, index1, element1, index2, element2):
    """Function to flat the peak
    """
    xpoints = [element1, element2]
    ypoints = [f_X[index1], f_X[index2]]

    coefficients = np.polyfit(xpoints, ypoints, 1)
    polynomial = np.poly1d(coefficients)
    y_axis = polynomial(X)

    return y_axis
    
    
def fitcurve(X, f_X, mask):
    """Function to flat the peak
    """
    xpoints = X[mask]
    ypoints = f_X[mask]

    coefficients = np.polyfit(xpoints, ypoints, 2)
    polynomial = np.poly1d(coefficients)
    y_axis = polynomial(xpoints)

    return y_axis
    
    
def find_nearest(array, value):
    """Function to find the nearest array element of a given value
    """
    
    index = (np.abs(array-value)).argmin()
    element = array[index]
    
    return (index, element)
    
    
def removePeaks(Q, I_Q):
    """Function to remove Bragg's peaks
    """
    # peakind = signal.find_peaks_cwt(I_Q, widths = np.arange(4,6))
    plt.figure('Remove Peaks')
    plt.plot(Q, I_Q)
    plt.grid()
    plt.xlabel('Q')
    plt.ylabel('I(Q)')
    
    points = np.array(plt.ginput(n=0, timeout=0, show_clicks=True, mouse_add=1, mouse_pop=3, mouse_stop=2))
    
    plt.close()
    
    idxelem = np.zeros(shape=(len(points),2))
    
    for i in range(0, len(points)):
        idxelem[i] = find_nearest(Q, points[i,0])
    
    zippedidx = np.array(list(zip(*[iter(idxelem[:,0])]*2)))
    zippedelem = np.array(list(zip(*[iter(idxelem[:,1])]*2)))
    
    
    for i in range(0, len(zippedelem)):
        mask = np.where((Q>=zippedelem[i,0]) & (Q<=zippedelem[i,1]))
        I_Q1 = fitline(Q[mask], I_Q, zippedidx[i,0], zippedelem[i,0], zippedidx[i,1], zippedelem[i,1])
        # I_Q1 = fitcurve(Q, I_Q, mask)
        I_Q[mask] = I_Q1
    
    return I_Q
    
    
def calc_damp(Q, QmaxIntegrate, damping_factor):
    """
    """
    
    # damping_factor = 0.5
    # damping_factor = np.log(10)
    exponent_factor = damping_factor / QmaxIntegrate**2
    damp_Q = np.exp(-exponent_factor * Q**2)
    
    return damp_Q
    
    
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



