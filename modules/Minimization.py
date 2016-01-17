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

"""Set of modules used in LASDiA to calculate all the used functions.

The nomenclature and the procedure follow the article:
Eggert et al. 2002 PRB, 65, 174105.

For the functions arguments and the returns I followed this convetion for the notes:
arguments: description - type
returns: description - type.

For the variables name I used this convention:
if the variable symbolizes a function, its argument is preceded by an underscore: f(x) -> f_x
otherwise it is just the name.
"""

import matplotlib.pyplot as plt

import sys
import os

import numpy as np
import scipy.constants as sc
from scipy import fftpack
from scipy.integrate import simps
from scipy.stats import chisquare

from modules.LiquidStructure import *

def calc_deltaMinim(N, r, Q, rho0, s, Sinf, I_Q, I_Qbkg, Iincoh, J_Q, fe_Q, Ztot, Fintra):
    """Function to minimize the density
    
    """
    #numAtoms = sc.N_A
    Isample = I_Q - s * I_Qbkg
    alpha = calc_alpha(J_Q, Sinf, Q, I_Q, fe_Q, Ztot, rho0)
    Icoh = (N * alpha * Isample) - (N * Iincoh)
    S_Q = calc_SQ(Icoh, Ztot, fe_Q)
    i_Q = calc_iQ(S_Q, Sinf)
    F_r = calc_Fr(r, Q, i_Q)
    observed = F_r
    excepted = Fintra - (4*np.pi*rho0)
    
    chi2, p = chisquare(observed, f_exp=excepted)
    return chi2
    
    



