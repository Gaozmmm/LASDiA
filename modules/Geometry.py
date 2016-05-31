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

"""Module containing the geometric correction: Diamond and Soller Slit.

The nomenclature and the procedure follow the article:
Weck et al. Rev. Sci. Instrum. 84, 063901 (2013)
Yaoita et al. Rev. Sci. Instrum. 68, 2106 (1997)

For the functions arguments and the returns I followed this convetion for the notes:
arguments: description - type
returns: description - type.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""

import numpy as np
from modules.UtilityAnalysis import Qto2theta
# from modules.Utility import write_file
from scipy.integrate import simps
import matplotlib.pyplot as plt

def diamond(dimension, Q, I_Q, diamond_angle):
    """Function to calculate the diamond correction.
    http://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html
    """
    
    # for now...
    wavelenght = 0.03738 # nm
    mu_l = 1/12.08 # mm^-1 at 33 keV
    _2theta = Qto2theta(Q) # rad
    diamond_angle = np.radians(diamond_angle)
    
    path_lenght = dimension / np.cos(_2theta - diamond_angle)
    
    corr_factor = np.exp(mu_l * path_lenght)
    
    I_Qeff = corr_factor * I_Q
    
    return (I_Qeff, corr_factor)


def calc_phi_angle(ws1, ws2, r1, r2, d, _2theta, sth):
    """Function to calculate dispersion angle for the Soller Slits correction.
    """
    
    gamma_2 = np.arcsin(ws1/(2*r1))
    
    # _2theta = np.radians(_2theta)
    
    alpha1 = np.arctan( r1 * np.sin(_2theta + gamma_2) / (r1*np.cos(_2theta + gamma_2) - sth ))
    alpha2 = np.arctan( r1 * np.sin(_2theta - gamma_2) / (r1*np.cos(_2theta - gamma_2) - sth ))
    
    beta1 = np.arctan( (r2+d) * np.sin(_2theta + gamma_2) / ((r2+d)*np.cos(_2theta + gamma_2) - sth ))
    beta2 = np.arctan( (r2+d) * np.sin(_2theta - gamma_2) / ((r2+d)*np.cos(_2theta - gamma_2) - sth ))
    
    psi = min(alpha1, beta1) - max(alpha2, beta2)
    phi = max(0, psi)
    
    # print(phi, psi)
    
    return (psi, phi)
    
    
def calc_T_MCC_samp(phi):
    """Function to calculate the MCC transmission for the sample.
    """
    
    print("phi ", phi.shape)
    DeltaPhi = np.diff(phi, axis=0)
    print("DeltaPhi ", DeltaPhi.shape)
    meanDeltaPhi = np.mean(DeltaPhi, axis=0)
    print("meanDeltaPhi ", meanDeltaPhi.shape)
    T_MCC_samp = np.sum(phi, axis=0) * meanDeltaPhi
    
    plt.figure("DeltaPhi")
    plt.plot(DeltaPhi[:,1999])
    plt.grid()
    plt.show()
        
    return T_MCC_samp
    
    
def calc_T_MCC_samp2(phi, sth):
    """Function to calculate the MCC transmission for the sample.
    """
    
    T_MCC_samp = simps(phi, sth, axis=0)
        
    return T_MCC_samp