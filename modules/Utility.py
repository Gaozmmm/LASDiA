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

"""Module containing useful functions used in LASDiA.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by
an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""


import sys
import os

import numpy as np
import math
import re
import imp
import time
import datetime
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d.axes3d import Axes3D

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget


def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)


def calc_distMol(x1, y1, z1, x2, y2, z2):
    """Function to calculate the distance between 2 points.

    Parameters
    ----------
    x1, x2 : float
             x coordinates of the two points (nm)
    y1, y2 : float
             y coordinates of the two points (nm)
    z1, z2 : float
             z coordinates of the two points (nm)

    Returns
    -------
    d       : float
              distance between two points (nm)
    """

    d = np.sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

    return d


def molToelemList(molecule):
    """Function to convert the molecule name to dictionary.

    Parameters
    ----------
    molecule    : string
                  molecule name

    Returns
    -------
    elementList : dictionary("element": multiplicity)
                  chemical elements of the sample with their multiplicity
                  element      : string
                                 chemical element
                  multiplicity : int
                                 chemical element multiplicity
    """

    elemlist = re.findall('[A-Z][a-z]*|\d+', re.sub('[A-Z][a-z]*(?![\da-z])', '\g<0>1', molecule))
    elementList = dict(zip(elemlist[0::2], elemlist[1::2]))
    elementList = dict((k,int(v)) for k,v in elementList.items())

    return elementList


def path_xyz_file(molecule):
    """Function for the GUI to determinate the xyz file path from the molecule name.

    Parameters
    ----------
    molecule : string
               molecule name

    Returns
    -------
    path     : string
               path of the xyz file
    """

    path = "./xyzFiles/" + molecule + ".xyz"

    return path


def plot_chi2(chi2, scale_factor, scale_factor_idx, rho0, rho0_idx):
    """Function to plot the 2D graph of chi2 with the scale factor and the atomic density.

    chi2             : 2D numpy array
                       chi2 values
    scale_factor     : numpy array
                       scale factor
    scale_factor_idx : int
                       scale factor best value index
    rho0             : numpy array
                       average atomic density
    rho0_idx         : int
                       atomic density best value index
    """

    # plot 2d chi2
    plt.figure('chi2s')
    plt.plot(scale_factor,chi2[rho0_idx, : ])
    plt.xlabel('scale factor')
    plt.ylabel('chi2')
    plt.grid()
    plt.show

    plt.figure('chi2rho0')
    plt.plot(rho0,chi2[ : ,scale_factor_idx])
    plt.xlabel('rho0')
    plt.ylabel('chi2')
    plt.grid()
    plt.draw()


def plot_data(xVal, yVal, plotName, xName, yName, labName, overlapping):
    """Function to plot the data.

    Parameters
    ----------
    xVal        : numpy array
                  abscissa values
    yVal        : numpy array
                  ordinate values
    plotName    : string
                  canvas name
    xName       : string
                  abscissa name
    yName       : string
                  ordinate name
    labName     : string
                  label name
    overlapping : string
                  flag for the overlapping
    """

    plt.figure(plotName)
    if overlapping.lower() == "n":
        plt.clf()
    plt.plot(xVal, yVal, label=labName)
    plt.xlabel(xName)
    plt.ylabel(yName)
    plt.legend()
    plt.grid(True)
    plt.draw()


def plot_data2(xVal, yVal, plotName, xName, yName, labName, overlapping):
    """Function to plot the data.

    Parameters
    ----------
    xVal        : numpy array
                  abscissa values
    yVal        : numpy array
                  ordinate values
    plotName    : string
                  canvas name
    xName       : string
                  abscissa name
    yName       : string
                  ordinate name
    labName     : string
                  label name
    overlapping : string
                  flag for the overlapping
    """

    fig = Figure()
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    ax.set_xlabel(xName)
    ax.set_ylabel(yName)
    ax.grid(True)
    
    ax.hold(True)
    ax.plot(xVal, yVal, label=labName)
    ax.legend()
    
    ui = QWidget()
    
    toolbar = NavigationToolbar(canvas, ui, coordinates=True)
    
    vbox = QVBoxLayout()
    vbox.addWidget(canvas)
    vbox.addWidget(toolbar)
    ui.setLayout(vbox)
    ui.show()


def plot_data_2scale(plotName, xVal1, yVal1, xName1, yName1, labName1,
    xVal2, yVal2, yName2, labName2):
    """Function to plot the data.

    Parameters
    ----------
    xVal        : numpy array
                  abscissa values
    yVal        : numpy array
                  ordinate values
    plotName    : string
                  canvas name
    xName       : string
                  abscissa name
    yName       : string
                  ordinate name
    labName     : string
                  label name
    overlapping : string
                  flag for the overlapping
    """
    
    fig, ax1 = plt.subplots()
    ax1.set_title(plotName)
    
    ax2 = ax1.twinx()
    # ax2.set_title(plotName)
    
    ax1.plot(xVal1, yVal1, 'b-', label=labName1)
    ax2.plot(xVal2, yVal2, 'g-', label=labName2)

    ax1.set_xlabel(xName1)
    ax1.set_ylabel(yName1, color='b')
    ax2.set_ylabel(yName2, color='g')
    
    # ax1.legend()
    # ax2.legend(loc=0)
    ax1.grid(True)
    # ax2.grid(True)
    
    align_yaxis(ax1, 0, ax2, 0)
    
    plt.draw()


def plot3d(xVal, yVal, zVal, plotName, xName, yName, zName):
    """Function to plot the 3D graph of chi2 and its profile.

    Parameters
    ----------
    xVal     : numpy array
               abscissa values
    yVal     : numpy array
               ordinate values
    zVal     : numpy array
               applicate values
    plotName : string
               canvas name
    xName    : string
               abscissa name
    yName    : string
               ordinate name
    zName    : string
               applicate name
    """

    # plot the 3d and its profile
    x, y = np.meshgrid(xVal, yVal)

    # # # fig = plt.figure('3D')
    # # # ax = Axes3D(fig)
    # # # ax.set_xlabel(xName)
    # # # ax.set_ylabel(yName)
    # # # ax.set_zlabel(zName)
    # # # ax.plot_surface(x, y, zVal, rstride=1, cstride=1, cmap='rainbow')

    # nLevel = int(chi2Max - chi2Min)
    # minIndxRho0, minIndxS = np.unravel_index(z_val.argmin(), z_val.shape)
    # maxIndxRho0, maxIndxS = np.unravel_index(z_val.argmax(), z_val.shape)
    # chi2Min = z_val[minIndxRho0][minIndxS]
    # chi2Max = z_val[maxIndxRho0][maxIndxS]
    # levels = np.linspace(0, chi2Max, int(chi2Max - chi2Min))

    plt.figure(plotName)
    plt.contour(xVal, yVal, zVal, 150)
    plt.colorbar()
    plt.xlabel(xName)
    plt.ylabel(yName)
    plt.draw()


def read_file(path):
    """Function to read the data file.
    This function can read a .chi and a .xy file.
    In our case the function returns are Q and I(Q).

    Parameters
    ----------
    path : string
           path of the data file

    Returns
    -------
    xVal : numpy array
           abscissa values
    yVal : numpy array
           ordinate values
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
    xVal = np.array(x)
    yVal = np.array(y)

    return (xVal, yVal)


def read_inputFile(path):
    """Function to read variables from the inputFile.txt.

    Parameters
    ----------
    path      : string
                path of the input file

    Returns
    -------
    variables : module
                input variables setted by the user
    """

    file = open(path)
    variables = imp.load_source("data", "", file)
    file.close()
    return variables


def read_MCC_file(path, type):
    """Function to read the MCC file with Soller Slits characteristics.

    Parameters
    ----------
    path : string
           path of the file
    type : string
           MCC model

    Returns
    -------
    ws1  : float
           width of the inner slit (cm)
    ws2  : float
           width of the outer slit (cm)
    r1   : float
           curvature radius of first slit (cm)
    r2   : float
           curvature radius of second slit (cm)
    d    : float
           slit thickness (cm)
    """

    file = open(path, "r")
    header1 = file.readline()
    lines = file.readlines()
    file.close()

    for line in lines:
        columns = line.split()
        if columns[0] == type:
            ws1 = float(columns[1])
            ws2 = float(columns[2])
            r1 = float(columns[3])
            r2 = float(columns[4])
            d = float(columns[5])
            break

    return (ws1, ws2, r1, r2, d)


def read_parameters(elementList, path):
    """Function to read the file containing the atomic form factor and incoherent parameters.

    Parameters
    ----------
    elementList       : dictionary("element": multiplicity)
                        chemical elements of the sample with their multiplicity
                        element      : string
                                       chemical element
                        multiplicity : int
                                       chemical element multiplicity
    path              : string
                        path of elements parameters

    Returns
    -------
    elementParameters : dictionary("element": parameters)
                        chemical elements of the sample with their parameters
                        element    : string
                                     chemical element
                        parameters : list
                                     list of the parameter
                                     (Z, a1, b1, a2, b2, a3, b3, a4, b4, c, M, K, L)
    """

    file = open(path, "r")
    header1 = file.readline()
    lines = file.readlines()
    file.close()

    elementParameters = {}

    for element, multiplicity in elementList.items():
        for line in lines:
            columns = line.split()
            if element not in elementParameters:
                if columns[0] == element:
                    Z = float(columns[1])
                    a1 = float(columns[2])
                    b1 = float(columns[3])
                    a2 = float(columns[4])
                    b2 = float(columns[5])
                    a3 = float(columns[6])
                    b3 = float(columns[7])
                    a4 = float(columns[8])
                    b4 = float(columns[9])
                    c = float(columns[10])
                    M = float(columns[11])
                    K = float(columns[12])
                    L = float(columns[13])
                    param = [Z, a1, b1, a2, b2, a3, b3, a4, b4, c, M, K, L]
                    elementParameters[element] = param
                    break

    return elementParameters


def read_xyz_file(path):
    """Function to read the xyz file.
    The atomic coordinates are in Angstrom, the function converts them in nanometer.

    Parameters
    ----------
    path             : string
                       path of the xyz file

    Returns
    -------
    numAtoms         : int
                       number of atoms in the molecule
    element          : string array
                       array of the elements in the molecule
    xPos, yPos, zPos : float array
                       atomic coordinate in the xyz_file (nm)
    """

    file = open(path, "r")
    file_name = file.name
    ext = os.path.splitext(file_name)[1]

    numAtoms = int(file.readline())
    comment = file.readline()

    lines = file.readlines()
    file.close()

    element = []
    x = []
    y = []
    z = []

    for line in lines:
        columns = line.split()
        element.append(columns[0])
        x.append(float(columns[1]))
        y.append(float(columns[2]))
        z.append(float(columns[3]))

    xPos = np.array(x)
    xPos /= 10.0
    yPos = np.array(y)
    yPos /= 10.0
    zPos = np.array(z)
    zPos /= 10.0
    # elementPosition = {}

    # for line in lines:
        # elem, x, y, z = line.split()
        # elementPosition[elem] = [x, y, z]

    return (numAtoms, element, xPos, yPos, zPos)


def setArray(minVal, maxVal, stepVal):
    """Function to generate the numpy array for the scale factor and the atomic density.

    Parameters
    ----------
    minVal   : float
               minimum value of the array
    maxVal   : float
               maximum value of the array
    stepVal  : float
               step of the array

    Returns
    -------
    arrayVal : numpy array
               array with setted values
    """

    if minVal == maxVal:
        arrayVal = np.arange(minVal, maxVal+stepVal, stepVal)
    else:
        arrayVal = np.arange(minVal, maxVal, stepVal)

    return arrayVal


def write_file(path, xVal, yVal, xName, yName):
    """Function to write on file.

    Parameters
    ----------
    path  : string
            path of the file
    xVal  : numpy array
            abscissa values
    yVal  : numpy array
            ordinate values
    xName : string
            abscissa name
    yName : string
            ordinate name
    """

    dir = os.path.dirname(path)

    if not os.path.exists(dir):
        os.makedirs(dir)

    file = open(path, "w")
    file_name = file.name

    output = np.column_stack((xVal.flatten(), yVal.flatten()))
    # file.write(xName + " \t " + yName + "\n")
    np.savetxt(file_name, output, delimiter='\t')
    file.close()


def write_results(path, molecule, scale_factor, rho0):
    """Function to write scale factor and atomic density on file.

    Parameters
    ----------
    path         : string
                   path of the file
    molecule     : string
                   molecule analyzed
    scale_factor : float
                   scale factor value
    rho0         : float
                   atomic density value
    """

    file = open(path, "a")
    file_name = file.name

    ts = time.time()
    timeStamp = datetime.datetime.fromtimestamp(ts).strftime("%Y%m%d_%H:%M:%S")

    file.write(timeStamp + " \t " + molecule + " \t " + str(scale_factor) + " \t " + str(rho0) + "\n")
    file.close()


def resize_zero(f_x, newDim):
    """Function to resize the function with zeros.
    
    Parameters
    ----------
    f_x    : numpy array
             function to resize
    
    newDim : int
             function final dimension
    
    Returns
    -------
    f_x    : numpy array
             resized function
    """

    for n in range(len(f_x), newDim):
        f_x = np.append(f_x, 0)
    
    return f_x