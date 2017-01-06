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

"""Set of modules used in LASDiA to calculate the chi2 minimization.

The nomenclature and the procedure follow the article:
Eggert et al. 2002 PRB, 65, 174105.

For the functions arguments and the returns I followed this convetion for the notes:
arguments: description - type
returns: description - type.

For the variables name I used this convention:
if the variable symbolizes a function, its argument is preceded by an underscore: f(x) -> f_x
otherwise it is just the name.
"""


import numpy as np
import matplotlib.pyplot as plt

from modules import KaplowMethod
from modules import UtilityAnalysis
import time


def chi2_fit(variableArray, chi2Array):
    """Function to fit the chi2 function to find the minimal value.
    
    Parameters
    ----------
    variableArray : numpy array
                    array with variable values
    chi2Array     : numpy array
                    array with chi2 values
    
    Returns
    -------
    x_fit         : numpy array
                    array with the x values to plot the fit results
    y_fit         : numpy array
                    array with the fit y values
    minValue      : float
                    minimum value of the fit curve
    """
    
    pol_fit = np.poly1d(np.polyfit(variableArray, chi2Array, 3))
    x_fit = np.linspace(variableArray[0], variableArray[-1], 1000)
    y_fit = pol_fit(x_fit)
    minValue = x_fit[np.argmin(y_fit)]
    
    return (x_fit, y_fit, minValue)


def chi2_minimization(scaleFactor, Q, I_Q, Ibkg_Q, J_Q, fe_Q, Iincoh_Q, Sinf, Ztot,
    density, Fintra_r, r, minQ, QmaxIntegrate, maxQ, smoothFactor, dampFactor, iteration, rmin):
    """Function to calculate the whole loop for the chi2 minimization.
    """
    
    scaleStep = 0.05
    densityStep = 0.025
    numSample = 23
    numLoopIteration = 0
    
    plt.ion()
    figure, ax = plt.subplots()

    while True:
        NoPeak = 0
        while True:
            ax.cla()
            ax.grid(True)
            scaleArray = UtilityAnalysis.make_array_loop(scaleFactor, scaleStep, numSample)

            chi2Array = np.zeros(numSample)

            plt.xlabel("Scale")
            ax.relim()
            ax.autoscale_view()
            for i in range(len(scaleArray)):
                chi2Array[i], SsmoothDamp_Q, F_r, Fopt_r = KaplowMethod.Kaplow_method(Q, I_Q,
                    Ibkg_Q, J_Q, fe_Q, Iincoh_Q, Sinf, Ztot, scaleArray[i], density, Fintra_r, r,
                    minQ, QmaxIntegrate, maxQ, smoothFactor, dampFactor, iteration, rmin)
                
                plt.scatter(scaleArray[i], chi2Array[i])
                figure.canvas.draw()
            
            if np.amax(chi2Array) > 10**8:
                chi2Array = chi2Array[0:np.amax(chi2Array)]
            
            scaleFactor = scaleArray[np.argmin(chi2Array)] - scaleStep*1.1
            
        ax.cla()
        ax.grid(True)

        density0 = density
        densityArray = UtilityAnalysis.make_array_loop(density, densityStep, numSample)
        chi2Array = np.zeros(numSample)
        
        plt.xlabel("Density")
        ax.relim()
        ax.autoscale_view()
        for i in range(len(densityArray)):
            chi2Array[i], SsmoothDamp_Q, F_r, Fopt_r = KaplowMethod.Kaplow_method(Q, I_Q,
                Ibkg_Q, J_Q, fe_Q, Iincoh_Q, Sinf, Ztot, scaleFactor, densityArray[i], Fintra_r, r,
                minQ, QmaxIntegrate, maxQ, smoothFactor, dampFactor, iteration, rmin)
            
            plt.scatter(densityArray[i], chi2Array[i])
            figure.canvas.draw()
            
        
        xfit, yfit, density = chi2_fit(densityArray, chi2Array)
        plt.plot(xfit, yfit)
        figure.canvas.draw()
        
        if np.abs(density-density0) > density0/25:
            scaleStep = 0.006
            densityStep = density0/10
        elif np.abs(density-density0) > density0/75:
            scaleStep = 0.0006
            densityStep = density0/100
        else:
            scaleStep = 0.00006
            densityStep = density0/1000

        numLoopIteration += 1
        if (np.abs(density-density0) < density0/2500 or numLoopIteration > 30):
            break
    
    plt.ioff()
    
    return (density, scaleFactor)
