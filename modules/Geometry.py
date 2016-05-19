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
Eggert et al. 2002 PRB, 65, 174105.

For the functions arguments and the returns I followed this convetion for the notes:
arguments: description - type
returns: description - type.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""

import numpy as np
from modules.UtilityAnalysis import Qto2theta

def diamond(path, type, Q, I_Q, diamond_angle):
    """Function to calculate the diamond correction.
    http://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html
    """
    
    file = open(path, "r")
    header1 = file.readline()
    header2 = file.readline()
    lines = file.readlines()
    file.close()
    
    for line in lines:
        columns = line.split()
        if columns[0] == type:
            dimension = float(columns[1])
            break
            
    # for now...
    # mu = 0.2562 # cm2/g
    wavelenght = 0.03738 # nm
    # diamond_density = 3.51 # g/cm3
    # mu_l = mu * diamond_density # cm
    mu_l = 1/12.08 # mm^-1 at 33 keV
    # theta_angle = np.arcsin(wavelenght*Q/(4*np.pi))
    _2theta_angle = Qto2theta(Q) # rad
    
    
    
    path_lenght = dimension / np.cos(_2theta_angle)
    
    corr_factor = np.exp(mu_l * path_lenght)
    
    I_Qeff = corr_factor * I_Q
    
    return (I_Qeff, corr_factor)
    
    
    
    
    
    
    
    
    
    
    