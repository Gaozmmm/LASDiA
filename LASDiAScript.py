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

"""LASDiA main script file.
This script is mainly used for testing the software, but it can be used to run LASDiA 
in text mode.

The nomenclature and the procedures follow the article: Eggert et al. 2002 PRB, 65, 174105.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by
an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""


from __future__ import (absolute_import, division, print_function, unicode_literals)
import six

import matplotlib.pyplot as plt
import numpy as np
from itertools import product
from timeit import default_timer as timer

from modules import Formalism
from modules import Geometry
from modules import IgorFunctions
from modules import KaplowMethod
from modules import MainFunctions
from modules import Minimization
from modules import Optimization
from modules import Utility
from modules import UtilityAnalysis


if __name__ == '__main__':
    variables = Utility.read_inputFile("./inputFile.txt")
    
    elementList = Utility.molToelemList(variables.molecule)
    elementParameters = Utility.read_parameters(elementList, variables.element_params)
    
    path = Utility.path_xyz_file(variables.molecule)
    numAtoms, element, x, y, z = Utility.read_xyz_file(path)
    
    Q, I_Q = Utility.read_file(variables.data_file)
    Qbkg, Ibkg_Q  = Utility.read_file(variables.bkg_file)
    
    Q, I_Q, Qbkg, Ibkg_Q = UtilityAnalysis.check_data_length(Q, I_Q, Qbkg, Ibkg_Q, \
        variables.minQ, variables.maxQ)
    
    I_Q, Ibkg_Q = Geometry.geometry_correction(Q, I_Q, Qbkg, Ibkg_Q, variables, "n")
    
    fe_Q, Ztot = MainFunctions.calc_eeff(elementList, Q, elementParameters)
    Iincoh_Q = MainFunctions.calc_Iincoh(elementList, Q, elementParameters)
    J_Q = MainFunctions.calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf = MainFunctions.calc_Sinf(elementList, fe_Q, Q, Ztot, elementParameters)
    
    r, Fintra_r = Optimization.calc_intraComponent(Q, fe_Q, Ztot, variables.QmaxIntegrate, \
        variables.maxQ, elementList, element, x, y, z, elementParameters, variables.damping_factor)
    
    scale_factor = variables.s_value
    density = variables.rho0_value
    percentage = 20
    
    # while True:
    sf_array = UtilityAnalysis.make_array(scale_factor, percentage)
    rho0_array = UtilityAnalysis.make_array(density, percentage)
    chi2 = np.zeros((rho0_array.size, sf_array.size))
    
    chi2_min = 1000000
    
    for (idx_rho0, rho0), (idx_sf, sf) in product(enumerate(rho0_array), enumerate(sf_array)):
        chi2[idx_rho0][idx_sf], S_Q, F_r = KaplowMethod.Kaplow_method(numAtoms, variables, \
            Q, I_Q, Ibkg_Q, J_Q, fe_Q, Iincoh_Q, Sinf, Ztot, sf, rho0, Fintra_r, r)
        if chi2[idx_rho0][idx_sf] < chi2_min:
            chi2_min = chi2[idx_rho0][idx_sf]
            best_sf = sf
            best_rho0 = rho0
            FOpt_r = F_r
            Sbest_Q = S_Q
    
    scale_factor, density = Minimization.calc_min_chi2(sf_array, rho0_array, chi2)
    # percentage /= 2
    # if conditions...
    
    print(scale_factor, best_sf)
    print(density, best_rho0)
    
    SOpt_Q = MainFunctions.calc_SQCorr(FOpt_r, r, Q, Sinf)
    
    Utility.plot_data(Q, Sbest_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S(Q)$", "y")
    Utility.plot_data(Q, SOpt_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S_{Opt}(Q)$", "y")
    
    plt.show()
    
    
    
    
    
    
    #-----------------------------------------------------
    
    # # test formalism
    
    # Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scale_factor, Ibkg_Q)
    # alpha = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, \
        # Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], \
        # fe_Q[Q<=variables.QmaxIntegrate], Ztot, density)
    # Icoh_Q = MainFunctions.calc_Icoh(numAtoms, alpha, Isample_Q, Iincoh_Q)
    
    # aff_sq_mean = Formalism.calc_aff_squared_mean(numAtoms, elementList, Q, elementParameters)
    # aff_mean_sq = Formalism.calc_aff_mean_squared(numAtoms, elementList, Q, elementParameters)
    
    # SFZ_Q = Formalism.calc_SFZ_Q(Q, Icoh_Q, aff_sq_mean, aff_mean_sq, variables.minQ, \
        # variables.QmaxIntegrate, variables.maxQ)
    # SAL_Q = Formalism.calc_SAL_Q(Q, Icoh_Q, aff_sq_mean, variables.minQ, variables.QmaxIntegrate, \
        # variables.maxQ)
    
    # S_Q = MainFunctions.calc_SQ(numAtoms, Icoh_Q, Ztot, fe_Q, Sinf, Q, variables.minQ, \
        # variables.QmaxIntegrate, variables.maxQ)
    
    # # Utility.plot_data(Q, SAL_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{AL}(Q)$", "y")
    # # Utility.plot_data(Q, SFZ_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{FZ}(Q)$", "y")
    # # Utility.plot_data(Q, S_Q, "S_Q", r"$Q(nm^{-1})$", r"$S(Q)$", r"$S^{E}(Q)$", "y")
    
    # # end formalism test
    
    #-----------------------------------------------------