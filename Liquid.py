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

if __name__ == '__main__':

    Q, I_Q = read_file("./data/cea_files/HT2_034.chi") # Mac
    Qbkg, I_Qbkg = read_file("./data/cea_files/HT2_036.chi") # Mac

    #-------------------------------------------------------------
    
    # plt.plot(Q, I_Q, "b-", label="Data")
    # plt.show()

    #-------------------------------------------------------------
    
    fQ = calc_aff("Ar", Q)
    
    # plt.plot(Q, fQ)
    # plt.show()
    
    #-------------------------------------------------------------
    
    elemList = {"Ar":1}
      
    fe, Ztot = calc_eeff(elemList, Q, calc_aff)
    # print(fe)
    # print(Ztot)
    # plt.plot(Q, fe)
    # plt.show()
    
    #-------------------------------------------------------------
    
    Iinc = calc_Iincoh(elemList, Q, calc_aff)
    # plt.plot(Q, Iinc)
    # plt.show()
    
    #-------------------------------------------------------------
    
    JQ = calc_JQ(Iinc, fe, Ztot)
    # plt.plot(Q, JQ)
    # plt.show()

    #-------------------------------------------------------------

    kp = calc_Kp(fe, "Ar", Q, calc_aff)
    # print(kp)

    #-------------------------------------------------------------

    sinf = calc_Sinf(elemList, fe, Q, Ztot, calc_Kp, calc_aff)
    # print(sinf)

    #-------------------------------------------------------------    

    shit = interpolateSpectra(Q, I_Q)
    plt.plot(Q, shit)
    plt.show()
    
    
    
    # print("Ztot: ", Ztot)
    # print("sum_Kp2: ", sum_Kp2)
    # print("Sinf: ", Sinf)

    # alpha = calc_alpha(JQ, sinf, Q, I_Q, fe, Ztot, 5)
    # print(alpha)
