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

import matplotlib.pyplot as plt
import numpy as np

from modules.Utility import *
from modules.LiquidStructure import *
from modules.InterpolateData import *
from modules.Optimization import *

if __name__ == '__main__':

    Q, I_Q = read_file("./data/cea_files/HT2_034.chi") # Mac
    Qbkg, I_Qbkg = read_file("./data/cea_files/HT2_036.chi") # Mac

    #-------------------------------------------------------------
    
    elemList = {"Ar":1}
      
    fe_Q, Ztot = calc_eeff(elemList, Q)
    
    #-------------------------------------------------------------
    
    Iinc = calc_Iincoh(elemList, Q)
    
    #-------------------------------------------------------------
    
    J_Q = calc_JQ(Iinc, fe_Q, Ztot)
    
    #-------------------------------------------------------------

    kp = calc_Kp(fe_Q, "Ar", Q)
    
    #-------------------------------------------------------------

    Sinf = calc_Sinf(elemList, fe_Q, Q, Ztot)
    # print(Sinf)
    
    #-------------------------------------------------------------    

    alpha = calc_alpha(J_Q, Sinf, Q, I_Q, fe_Q, Ztot)
    Isample = I_Q - I_Qbkg
    
    S_Q = calc_SQ(alpha, Isample, Ztot, fe_Q)
    # plt.figure(1)
    # plt.plot(Q, S_Q)
    # plt.grid()
    # plt.show

    #--------------------------------------------------------------
    
    i_Q = calc_iQ(S_Q, Sinf)

    r, F_r = calc_Fr(Q, i_Q)
    # plt.plot(r, F_r)
    # plt.grid()
    # plt.show()
    
    #--------------------------------------------------------------
    # calc of rho0
       
    rho0 = calc_rho0(Q, i_Q)
    # print(rho0)
    
    #--------------------------------------------------------------
    rho0 = 25.0584
    r_cutoff = 0.2
    
    iteration = 2
    opt_Fr = calc_optimize_Fr(iteration, F_r, Fintra_r, rho0, i_Q, Q, Sinf, J_Q, r, r_cutoff)
    
    plt.figure(1)
    plt.plot(r, opt_Fr)
    plt.grid()
    plt.show()

    
    # #-------------------------------------------------------------
    # Very important for the FFT!!!
    # r = fftpack.fftfreq(S_Q.size, Q[1] - Q[0])    
    # F_r = fftpack.fft(S_Q)
    # F_rI = F_r.imag

    # pidxs = np.where(r > 0)
    # r = r[pidxs]
    # F_rI = F_rI[pidxs]
    
    # print(r.size, F_rI.size)
    
    # plt.figure(2)
    # plt.plot(r, F_rI)
    # plt.grid()
    # plt.show()
    
    # #-------------------------------------------------------------    

    # QiQ = Q*(S_Q - Sinf)
    
    # r, F_r = FFT_QiQ(QiQ)