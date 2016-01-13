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

    right_index = np.where((Q>minQ) & (Q<QmaxIntegrate))
    # inf_index = np.where((Q>QmaxIntegrate) & (Q<maxQ))
    newQ = Q[right_index]
    newI_Q = I_Q[right_index]
    newI_Qbkg = I_Qbkg[right_index]
    
    elementList = {"Ar":1}

    # remember the electronic unit in atomic form factor!!!
    fe_Q, Ztot = calc_eeff(elementList, newQ)
    
    s =  0.5
    rho0 = 24.00
    
    Iincoh_Q = calc_Iincoh(elementList, newQ)
    J_Q = calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf = calc_Sinf(elementList, fe_Q, newQ, Ztot)
    
    Isample_Q = newI_Q - s * newI_Qbkg

    alpha = calc_alpha(J_Q, Sinf, newQ, Isample_Q, fe_Q, Ztot, rho0)
    
    Icoh_Q = calc_Icoh(N, alpha, Isample_Q, Iincoh_Q)

    S_Q = calc_SQ(N, Icoh_Q, Ztot, fe_Q)

    DeltaQ = np.diff(Q)
    meanDeltaQ = np.mean(DeltaQ)
    # Qinf = Q[inf_index]
    Qinf = np.arange(QmaxIntegrate, maxQ, meanDeltaQ)
    Qzero = np.arange(0.0, minQ, meanDeltaQ)
    newQinf = np.concatenate([Qzero, newQ])
    newQinf = np.concatenate([newQinf, Qinf])

    SQinf = np.zeros(Qinf.size)
    SQinf.fill(Sinf)
    SQzero = np.zeros(Qzero.size)
    newSQinf = np.concatenate([S_Q, SQinf])
    newSQinf = np.concatenate([SQzero, newSQinf])
    
    plt.figure(1)
    plt.plot(newQinf,newSQinf)
    plt.grid()
    plt.show

    i_Q = calc_iQ(newSQinf, Sinf)
    Qi_Q = newQinf * i_Q
    
    F_r = fftpack.fft(Qi_Q)
    # num = Qi_Q.size
    F_rI = F_r.imag
    F_rI *= (2/np.pi*meanDeltaQ)
    # f_s = newQinf.size * meanDeltaQ
    r = fftpack.fftfreq(newQinf.size, meanDeltaQ)
    mask = np.where(r>0)

    plt.figure(2)
    plt.plot(r[mask], F_rI[mask])
    plt.grid()
    plt.show

    r3 = np.linspace(0.0, 9.0, newQinf.size)
    F_r3 = calc_Fr(r[mask], newQinf, i_Q)
   
    plt.figure(3)
    plt.plot(r[mask], F_r3)
    plt.grid()
    plt.show()