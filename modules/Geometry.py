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

def calc_absorption_correction(abs_length, _2theta, thickness, I_Q, angle):
    """Function to calculate the absorption correction.
    This function can be used to calculate the absorption correction for the diamond 
    or for any other object between the sample and the detector.
    
    The characteristics for some diamonds can be found here:
    http://www.almax-easylab.com/DiamondSelectionPage.aspx
    
    Parameters
    ----------
    abs_length  : float
                  absorption length  (cm), @33keV 1.208cm
    _2theta     : numpy array
                  diffraction angle (rad)
    thickness   : float
                  object thickness (cm)
    I_Q         : numpy array
                  measured scattering intensity
    angle       : float
                  object rotation angle respect the XRay beam (deg)
    
    Returns
    -------
    I_Qeff      : numpy array
                  corrected scattering intensity (nm)
    corr_factor : numpy array
                  correction factor
    """
    
    # for now...
    # wavelenght  : float
    #               XRay beam wavelenght (nm), @ESRF ID27 0.03738nm
    # wavelenght = 0.03738 # nm
    # _2theta = Qto2theta(Q) # rad
    
    mu_l = 1/abs_length
    angle = np.radians(angle)
    
    path_lenght = thickness / np.cos(_2theta - angle)
    
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


def calc_phi_matrix(thickness, _2theta, ws1, ws2, r1, r2, d, num_point):
    """Function to calculate the dispersion angle matrix.
    half_thick in cm
    
    Parameters
    ----------
    thickness_value    : float
                         object thickness (sample or sample+DAC)
    _2theta            : numpy array
                         diffraction angle (rad)
    ws1                : float
                         width of the inner slit (cm)
    ws2                : float
                         width of the outer slit (cm)
    r1                 : float
                         curvature radius of first slit (cm)
    r2                 : float
                         curvature radius of second slit (cm)
    d                  : float
                         slit thickness (cm)
    num_point          : int
                         number of point for the thickness array
    
    Returns
    -------
    thickness_sampling : numpy array
                         array with the thickness values
    phi_matrix         : 2D numpy array
                         dispersion angle matrix (rad)
    """
    
    thickness_sampling = np.linspace(-thickness, thickness, num=num_point) # num=500)
    phi_matrix = np.zeros((thickness_sampling.size, _2theta.size))
    
    for i, val_sth in enumerate(thickness_sampling):
        for j, val_2theta in enumerate(_2theta):
            phi_matrix[i][j] = calc_phi_angle(ws1, ws2, r1, r2, d, _2theta[j], thickness_sampling[i])
    
    return (thickness_sampling, phi_matrix)


def calc_T_MCC(sample_thickness, all_thickness_sampling, all_phi_angle_matrix):
    """Function to calculate the MCC transfer function for the sample, the DAC and sample+DAC.
    
    Parameters
    ----------
    sample_thickness             : float
                                   sample thickness
    all_thickness_sampling       : numpy array
                                   array with the thickness values for sample+DAC
    all_phi_angle_matrix         : 2D numpy array
                                   dispersion angle matrix for sample+DAC (rad)
    
    Returns
    -------
    T_MCC_sample                 : numpy array
                                   MCC sample transfer function
    T_MCC_DAC                    : numpy array
                                   MCC DAC transfer function
    T_MCC_ALL                    : numpy array
                                   MCC sample+DAC transfer function
    all_thickness_sampling[mask] : numpy array
                                   array with the thickness values for sample
    all_phi_angle_matrix[mask]   : 2D numpy array
                                   dispersion angle matrix for sample (rad)
    """
    
    mask = (all_thickness_sampling>= -sample_thickness/2) & (all_thickness_sampling <=sample_thickness/2)
    
    T_MCC_ALL = simps(all_phi_angle_matrix, axis=0, even="first")
    T_MCC_sample = simps(all_phi_angle_matrix[mask], axis=0, even="first")
    T_MCC_DAC = T_MCC_ALL - T_MCC_sample
        
    return (T_MCC_sample, T_MCC_DAC, T_MCC_ALL, all_thickness_sampling[mask], all_phi_angle_matrix[mask])