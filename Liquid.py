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

from modules.MainModules import *

if __name__ == '__main__':
#    Q, I_Q = read_file(".\data\cea_files\HT2_034.chi") # Win
    Q, I_Q = read_file("./data/cea_files/HT2_034.chi") # Mac
    # Q, I_Q = read_file(".\data\my_test\HT2_034_20151104.chi")
    #plot_data(1, Q, I_Q, "Q", "I(Q)", "b-", "Data", False)
    
#    Qbkg, I_Qbkg = read_file(".\data\cea_files\HT2_036.chi") # Win
    Qbkg, I_Qbkg = read_file("./data/cea_files/HT2_036.chi") # Mac
    # Qbkg, I_Qbkg = read_file(".\data\my_test\HT2_036_20151104.chi")
    #plot_data(1, Qbkg, I_Qbkg, "Q", "I(Q)", "k--", "Bkg", True)

    #-------------------------------------------------------------
    
    # plt.plot(Q, I_Q, "b-", label="Data")
    # plt.plot(Qbkg, I_Qbkg, "k--", label="Bkg")
    # plt.legend()
    # plt.xlabel("Q")
    # plt.ylabel("I(Q)")
    # plt.show()

    #-------------------------------------------------------------
    
    # fakeQ = np.arange(0.0,25.0,1.0)
    # Q /= 10 # conversion from inverse nm to inverse A
    # fQ = calc_aff("Ar", Q)
    
    # print(fQ.shape)
    
    # plt.plot(Q, fQ)
    # plt.show()
    
    #-------------------------------------------------------------
    
    elemList = {"Ar":1}
      
    fe, Ztot = calc_eeff(elemList, Q, calc_aff)
#    print(fe)
    # print(Ztot)
    # plt.plot(Q, fe)
    # plt.show()
    
    #-------------------------------------------------------------
    
    Iinc = calc_Iincoh(elemList, Q, calc_aff)
 #   print(Iinc)
    # plt.plot(Q, Iinc)
    # plt.show()
    
    #-------------------------------------------------------------
    
    JQ = calc_JQ(Iinc, fe, Ztot)
 #   print(JQ)
    # plt.plot(Q, JQ)
    # plt.show()

    #-------------------------------------------------------------

    kp = calc_Kp(fe, "Ar", Q, calc_aff)
#    print(kp)

    #-------------------------------------------------------------

    sinf = calc_Sinf(elemList, fe, Q, Ztot, calc_Kp, calc_aff)
    print(sinf)

    #-------------------------------------------------------------    

    # print("Ztot: ", Ztot)
    # print("sum_Kp2: ", sum_Kp2)
    # print("Sinf: ", Sinf)

    alpha = calc_alpha(JQ, sinf, Q, I_Q, fe, Ztot, 5)
    print(alpha)
