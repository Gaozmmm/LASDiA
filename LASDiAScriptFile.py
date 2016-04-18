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

"""LASDiA main script file.
This script is mainly used for testing the software, but it can be used to run LASDiA in text mode.
The nomenclature and the procedures follow the article: Eggert et al. 2002 PRB, 65, 174105
"""

from __future__ import (absolute_import, division, print_function, unicode_literals)
import six

import sys
import os

import matplotlib.pyplot as plt

from modules.MainFunctions import *
from modules.Utility import *
# from modules.UtilityAnalysis import *
# from modules.Optimization import *
# from modules.Minimization import *
# from modules.Formalism import *
# from modules.IgorFunctions import *


if __name__ == "__main__":
    N = 1 # sc.N_A

    variables = read_inputFile("inputFile.txt")

    elementList = molToelemList(variables.molecule)
    Q, I_Q = read_file(variables.data_file)
    Qbkg, I_Qbkg = read_file(variables.bkg_file)

    fe_Q, Ztot = calc_eeff(elementList, Q)
    Iincoh_Q = calc_Iincoh(elementList, Q)
    J_Q = calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf = calc_Sinf(elementList, fe_Q, Q, Ztot)

    min_index, max_index = calc_indices(Q, variables.minQ, variables.QmaxIntegrate, variables.maxQ)
    validation_index, integration_index, calculation_index = calc_ranges(Q, variables.minQ, variables.QmaxIntegrate, variables.maxQ)

    scale_factor = setArray(variables.s_min, variables.s_max, variables.s_step)
    rho0 = setArray(variables.rho0_min, variables.rho0_max, variables.rho0_step)

    chi2 = np.zeros((rho0.size, scale_factor.size))

    for i, val_rho0 in enumerate(rho0):
        for j, val_s in enumerate(scale_factor):
            Isample_Q = calc_IsampleQ(I_Q, scale_factor[j], I_Qbkg)
            alpha = calc_alpha(J_Q[integration_index], Sinf, Q[integration_index], Isample_Q[integration_index], fe_Q[integration_index], Ztot, rho0[i])
            Icoh_Q = calc_Icoh(N, alpha, Isample_Q, Iincoh_Q)

            S_Q = calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, max_index, integration_index)
            newQ, S_Qsmoothed = calc_SQsmoothing(Q[validation_index], S_Q[validation_index], Sinf, variables.smooth_factor, min_index, variables.minQ, variables.QmaxIntegrate, variables.maxQ, 550)
            S_QsmoothedDamp = calc_SQdamp(S_Qsmoothed, newQ, Sinf, variables.QmaxIntegrate, variables.damp_factor)

            plot_data(newQ, S_QsmoothedDamp, "Q($nm^{-1}$)", "S(Q)", "S_Q")
            plot_data(Q[validation_index], S_Q, "Q($nm^{-1}$)", "S(Q)", "S_Q")

            # plt.figure("S_Q")
            # plt.plot(Q[validation_index], S_Q, label="S_Q")
            # plt.plot(newQ, S_QsmoothedDamp, label="S_Q")
            # plt.xlabel("Q")
            # plt.ylabel("S(Q)")
            # plt.legend()
            # plt.grid()
            # plt.show()
