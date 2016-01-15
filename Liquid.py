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

"""Environment where to test the new modules
The nomenclature and the procedures follow the article: Eggert et al. 2002 PRB, 65, 174105
"""

from __future__ import (absolute_import, division, print_function, unicode_literals)
import six

import sys
import os

import scipy.constants as sc
from scipy import fftpack
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

from modules.Utility import *
from modules.LiquidStructure import *
from modules.InterpolateData import *
from modules.Optimization import *
from modules.Minimization import *

# import cmath
# from cmath import exp, pi

if __name__ == '__main__':
    N = 1 # sc.N_A
    
    Q, I_Q = read_file("./data/cea_files/HT2_034T++.chi")
    Qbkg, I_Qbkg = read_file("./data/cea_files/HT2_036T++.chi")
    
    minQ = 3
    maxQ = 109
    QmaxIntegrate = 90
    
    operative_index = np.where((Q>minQ) & (Q<=maxQ))
    operative_Q = Q[operative_index]
    operative_I_Q = I_Q[operative_index]
    operative_I_Qbkg = I_Qbkg[operative_index]
    
    min_index = np.where(Q<=minQ)
    max_index = np.where((Q>QmaxIntegrate) & (Q<=maxQ))
    validate_index = np.where(Q<=maxQ)
    
    integrate_index = np.where((Q>minQ) & (Q<=QmaxIntegrate))
    integrate_Q = Q[integrate_index]
    integrate_I_Q = I_Q[integrate_index]
    integrate_I_Qbkg = I_Qbkg[integrate_index]
    
    elementList = {"Ar":1}
    s =  0.5
    rho0 = 24.00
    
    # remember the electron unit in atomic form factor!!!
    fe_Q, Ztot = calc_eeff(elementList, Q)
    Iincoh_Q = calc_Iincoh(elementList, Q)
    J_Q = calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf = calc_Sinf(elementList, fe_Q, Q, Ztot)
    Isample_Q = I_Q - s * I_Qbkg
    
    alpha = calc_alpha(J_Q, Sinf, Q, Isample_Q, fe_Q, Ztot, rho0, integrate_index)
    Icoh_Q = calc_Icoh(N, alpha, Isample_Q, Iincoh_Q)
    
    S_Q = calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, min_index, max_index, integrate_index)
    
    plt.figure(1)
    plt.plot(Q[validate_index], S_Q)
    plt.grid()
    plt.show
    
    i_Q = calc_iQ(S_Q, Sinf)
    Qi_Q = Q[validate_index] * i_Q
    
    DeltaQ = np.diff(Q)
    meanDeltaQ = np.mean(DeltaQ)
    r = fftpack.fftfreq(Q[validate_index].size, meanDeltaQ)
    mask = np.where(r>0)
    
    F_r = calc_Fr(r[mask], Q[validate_index], i_Q)
    
    plt.figure(2)
    plt.plot(r[mask], F_r)
    plt.grid()
    plt.show
    
    iteration = 2
    r_cutoff = 0.25
    F_rInt = calc_optimize_Fr(iteration, F_r, rho0, i_Q, Q[validate_index], Sinf, J_Q[validate_index], r[mask], r_cutoff)
    
    plt.figure(2)
    plt.plot(r[mask], F_rInt)
    plt.show()