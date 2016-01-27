# The MIT License (MIT)

# Copyright (c) 2016 Francesco Devoto

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

"""Environment for testing the software
The nomenclature and the procedures follow the article: Eggert et al. 2002 PRB, 65, 174105
"""

from __future__ import (absolute_import, division, print_function, unicode_literals)
import six

import sys
import os

import scipy.constants as sc
from scipy import fftpack
from scipy import signal
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
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

    # s = UnivariateSpline(Q, I_Q, k=3, s=0.5)
    # I_Qs = s(Q)

    # plt.figure(3)
    # plt.plot(Q, I_Q)
    # plt.plot(Q, I_Qs)
    # plt.grid()
    # plt.show()
    
    minQ = 3
    maxQ = 109
    QmaxIntegrate = 90
    
    min_index = np.where(Q<=minQ)
    max_index = np.where((Q>QmaxIntegrate) & (Q<=maxQ))
    validation_index = np.where(Q<=maxQ)
    integration_index = np.where(Q<=QmaxIntegrate)
    
    calculation_index = np.where((Q>minQ) & (Q<=QmaxIntegrate))
    
    elementList = {"Ar":1}
    # test values
    # s = np.arange(0.2, 0.8, 0.1)
    # rho0 = np.arange(24, 28, 1)
    
    # real values
    # s = np.arange(0.2, 0.8, 0.01)
    # rho0 = np.arange(24.0, 28.0, 0.1)
    
    # best values
    s = np.array([0.57])
    rho0 = np.array([26.10])
    
    chi2 = np.zeros((rho0.size, s.size))
    
    # remember the electron unit in atomic form factor!!!
    fe_Q, Ztot = calc_eeff(elementList, Q)
    Iincoh_Q = calc_Iincoh(elementList, Q)
    J_Q = calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf = calc_Sinf(elementList, fe_Q, Q, Ztot)

    for i, val_rho0 in enumerate(rho0):
        for j, val_s in enumerate(s):
            # print(rho0[i], s[j])
            Isample_Q = calc_IsampleQ(I_Q, s[j], I_Qbkg)
            
            # plt.figure(1)
            # plt.plot(Q, Isample_Q)
            # plt.grid()
            # plt.show()
            
            alpha = calc_alpha(J_Q, Sinf, Q, Isample_Q, fe_Q, Ztot, rho0[i], integration_index)
            Icoh_Q = calc_Icoh(N, alpha, Isample_Q, Iincoh_Q)
            
            S_Q, S_Qs = calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, min_index, max_index, calculation_index)
            
            plt.figure(1)
            plt.plot(Q[validation_index], S_Q)
            # plt.plot(Q[validation_index], S_Qs)
            plt.grid()
            plt.show
            
            i_Q = calc_iQ(S_Q, Sinf)
        #    i_Q = calc_iQ(S_Qs, Sinf)
            Qi_Q = Q[validation_index] * i_Q
            
            DeltaQ = np.diff(Q)
            meanDeltaQ = np.mean(DeltaQ)
            r = fftpack.fftfreq(Q[validation_index].size, meanDeltaQ)
            mask = np.where(r>0)
            
            F_r = calc_Fr(r[mask], Q[integration_index], i_Q[integration_index])
            
            plt.figure(2)
            plt.plot(r[mask], F_r)
            plt.grid()
            plt.show
            
            iteration = 2
            rmin = 0.24
            F_rInt = calc_optimize_Fr(iteration, F_r, rho0[i], i_Q[integration_index], Q[integration_index], Sinf, J_Q[integration_index], r[mask], rmin)
            
            plt.figure(2)
            plt.plot(r[mask], F_rInt)
            plt.show()
            
            maskInt = np.where((r>0) & (r < rmin))
            rInt = r[maskInt]
            Fintra_r = calc_Fintra(rInt)
            deltaF_r = calc_deltaFr(F_rInt[maskInt], Fintra_r, rInt, rho0[i])
            chi2[i][j] = simps(deltaF_r**2, rInt)
        
    minIndxRho0, minIndxS = np.unravel_index(chi2.argmin(), chi2.shape)
    # print(chi2[minIndxRho0][minIndxS])
    # print(rho0[minIndxRho0], s[minIndxS])
    
    # x, y = np.meshgrid(s, rho0)
    # fig = plt.figure(3)
    # ax = Axes3D(fig)
    # ax.plot_surface(x, y, chi2, rstride=1, cstride=1, cmap='rainbow')
    
    # plt.figure(4)
    # plt.contour(s, rho0, chi2, 200)
    
    # plt.show()
    

    
