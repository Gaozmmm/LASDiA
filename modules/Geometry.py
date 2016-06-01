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

"""Set of modules used in LASDiA to calculate the diamond absorption and Soller Slits correction.

The nomenclature and the procedure follow the articles:
Weck et al. Rev. Sci. Instrum. 84, 063901 (2013)
Yaoita et al. Rev. Sci. Instrum. 68, 2106 (1997)

For the functions arguments and the returns I followed this convetion for the notes:
argument: description - type
return: description - type.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""

import numpy as np
from modules.UtilityAnalysis import Qto2theta
# from modules.Utility import write_file
from scipy.integrate import simps
import matplotlib.pyplot as plt

def calc_diamond_absorption(wavelenght, abs_length, _2theta, dimension, Q, I_Q, diamond_angle):
    """Function to calculate the diamond absorption correction.
    Diamond characteristics can be found here:
    http://www.almax-easylab.com/DiamondSelectionPage.aspx
    
    Parameters
    ----------
    wavelenght    : float
                    XRay beam wavelenght, @ESRF ID27 0.03738nm (nm)
    abs_length    : float
                    absorption length, @33keV 12.08mm (mm)
    _2theta       : numpy array
                    diffraction angle (rad)
    dimension     : float
                    diamond dimension (mm)
    Q             : numpy array
                    momentum transfer (nm)
    I_Q           : numpy array
                    measured scattering intensity (nm)
    diamond_angle : float
                    diamond rotation angle respect the XRay beam (deg)
    
    
    Returns
    -------
    I_Qeff        : numpy array
                    corrected scattering intensity (nm)
    corr_factor   : numpy array
                    correction factor
    """
    
    # for now...
    # wavelenght = 0.03738 # nm
    # _2theta = Qto2theta(Q) # rad
    
    mu_l = 1/abs_length # mm^-1 at 33 keV
    diamond_angle = np.radians(diamond_angle)
    
    path_lenght = dimension / np.cos(_2theta - diamond_angle)
    
    corr_factor = np.exp(mu_l * path_lenght)
    
    I_Qeff = corr_factor * I_Q
    
    return (I_Qeff, corr_factor)


def calc_phi_angle(ws1, ws2, r1, r2, d, _2theta, xth):
    """Function to calculate the dispersion angle for the Soller Slits correction.
    
    Parameters
    ----------
    ws1     : float
              width of the inner slit (cm)
    ws2     : float
              width of the outer slit (cm)
    r1      : float
              curvature radius of first slit (cm)
    r2      : float
              curvature radius of second slit (cm)
    d       : float
              slit thickness (cm)
    _2theta : float
              diffraction angle (rad)
    xth     : float
              positin of i-th point on x-axis (cm)
    
    
    Returns
    -------
    phi     : float
              dispersion angle (rad)
    """
    
    gamma_2 = np.arcsin(ws1/(2*r1))
    
    alpha1 = np.arctan( r1 * np.sin(_2theta + gamma_2) / (r1*np.cos(_2theta + gamma_2) - xth ))
    alpha2 = np.arctan( r1 * np.sin(_2theta - gamma_2) / (r1*np.cos(_2theta - gamma_2) - xth ))
    
    beta1 = np.arctan( (r2+d) * np.sin(_2theta + gamma_2) / ((r2+d)*np.cos(_2theta + gamma_2) - xth ))
    beta2 = np.arctan( (r2+d) * np.sin(_2theta - gamma_2) / ((r2+d)*np.cos(_2theta - gamma_2) - xth ))
    
    psi = min(alpha1, beta1) - max(alpha2, beta2)
    phi = max(0, psi)
    
    return phi
    
    
def calc_phi_matrix(thickness, _2theta, ws1, ws2, r1, r2, d):
    """Function to calculate the dispersion angle matrix.
    half_thick in cm
    
    Parameters
    ----------
    thickness  : float
                 thickness of the sample (or sample+DAC)
    _2theta    : numpy array
                 diffraction angle (rad)
    ws1        : float
                 width of the inner slit (cm)
    ws2        : float
                 width of the outer slit (cm)
    r1         : float
                 curvature radius of first slit (cm)
    r2         : float
                 curvature radius of second slit (cm)
    d          : float
                 slit thickness (cm)
    
    
    Returns
    -------
    phi_matrix : 2D numpy array
                 dispersion angle matrix (rad)
    """
    
    sth = np.linspace(-thickness, thickness, num=100)
    phi_matrix = np.zeros((sth.size, _2theta.size))
    
    for i, val_sth in enumerate(sth):
        for j, val_theta2 in enumerate(_2theta):
            phi_matrix[i][j] = calc_phi_angle(ws1, ws2, r1, r2, d, _2theta[j], sth[i])
    
    return phi_matrix
    
    
def calc_T_MCC_sample(phi_matrix):
    """Function to calculate the MCC sample transfer function.
    
    Parameters
    ----------
    phi_matrix : 2D numpy array
                 dispersion angle matrix (rad)
    
    
    Returns
    -------
    T_MCC_sample : numpy array
                 MCC sample transfer function
    """
    
    T_MCC_sample = simps(phi_matrix, axis=0, even="first")
        
    return T_MCC_sample
    
    
def calc_T_MCC_DAC(phi_matrix, T_MCC_samp):
    """Function to calculate the MCC DAC transfer function.
    
    Parameters
    ----------
    phi_matrix : 2D numpy array
                 dispersion angle matrix (rad)
    T_MCC_samp : numpy array
                 MCC sample transfer function
    
    
    Returns
    -------
    T_MCC_DAC  : numpy array
                 MCC DAC transfer function
    """
    
    T_MCC_ALL = simps(phi_matrix, axis=0, even="first") 
    T_MCC_DAC = T_MCC_ALL - T_MCC_samp
        
    return (T_MCC_ALL, T_MCC_DAC)