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
    
    
def rebinning(X, maxQ):
    """
    """
    
    lenX = len(X)
    numX = 2**int(math.log(lenX,2))
    rebinnedX = np.linspace(np.amin(X), maxQ, numX, endpoint=True)
    
    return rebinnedX
    
def interpolation(X, f_X):
    """
    """
    
    s = UnivariateSpline(X, f_X, k=3, s=0.5)
    newf_X = s(X)
    
    return newf_X
    
    
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



