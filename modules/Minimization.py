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


def calc_min_chi2(scale_factor, rho0, chi2):
    """Function to calculate the minimum of chi2 matrix
    
    Parameters
    ----------
    scale_factor                : numpy array
                                  scale factor
    rho0                        : numpy array
                                  average atomic density
    chi2                        : 2D numpy array
                                  chi2 values
    
    
    Returns
    -------
    scale_factor[minIndxS]      : float
                                  scale factor minimum value
    rho0[minIndxRho0]           : float
                                  atomic density minimum value
    """
    
    minIndxRho0, minIndxS = np.unravel_index(chi2.argmin(), chi2.shape)
    
    return (scale_factor[minIndxS], rho0[minIndxRho0])


def chi2_fit(density_array, chi2_array):
    """Function to fit the chi2 function to find the minimal value.
    
    Parameters
    ----------
    density_array : numpy array
                    array with density values
    chi2_array    : numpy array
                    array with chi2 values
    
    Returns
    -------
    x_fit          : numpy array
                     array with the x values to plot the fit results
    y_fit          : numpy array
                     array with the fit y values
    min_density    : float
                     minimum value of the fit curve
    """
    
    pol_fit = np.poly1d(np.polyfit(density_array, chi2_array, 2))
    x_fit = np.linspace(density_array[0], density_array[-1], 1000)
    y_fit = pol_fit(x_fit)
    min_density = x_fit[np.argmin(y_fit)]
    
    return (x_fit, y_fit, min_density)
    
    