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
This script is mainly used for testing the software, but it can be used to run 
LASDiA in text mode.

The nomenclature and the procedures follow the article:
    Eggert et al. 2002 PRB, 65, 174105.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by
an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""


from __future__ import (absolute_import, division, print_function, unicode_literals)

import matplotlib.pyplot as plt
import numpy as np


# from modules import Formalism
from modules import Geometry
from modules import IgorFunctions
# from modules import KaplowMethod
from modules import MainFunctions
from modules import Minimization
from modules import Optimization
from modules import Utility
from modules import UtilityAnalysis

#from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget


if __name__ == "__main__":
    
    #---------------------------Files reading----------------------------------
    
    variables = Utility.read_inputFile("./inputFile.txt")
    
    elementList = Utility.molToElemList(variables.molecule)
    elementParameters = Utility.read_parameters(elementList, variables.element_params_path)
    
    path = Utility.path_xyz_file(variables.molecule)
    numAtoms, element, x, y, z = Utility.read_xyz_file(path)
    
    Q, I_Q = Utility.read_file(variables.data_file)
    Qbkg, Ibkg_Q  = Utility.read_file(variables.bkg_file)

    #--------------------Preliminary calculation-------------------------------
    
    # Q, I_Q, Qbkg, Ibkg_Q = UtilityAnalysis.check_data_length(Q, I_Q, Qbkg, Ibkg_Q, \
        # variables.minQ, variables.maxQ)
    
    fe_Q, Ztot = MainFunctions.calc_eeff(elementList, Q, elementParameters)
    Iincoh_Q = MainFunctions.calc_Iincoh(elementList, Q, elementParameters)
    J_Q = MainFunctions.calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf = MainFunctions.calc_Sinf(elementList, fe_Q, Q, Ztot, elementParameters)

    dampingFunction = UtilityAnalysis.calc_dampingFunction(Q, variables.dampingFactor,
        variables.QmaxIntegrate, variables.typeFunction)
        
    #-------------------Intra-molecular components-----------------------------
    
    iintra_Q = Optimization.calc_iintra(Q, fe_Q, Ztot, variables.QmaxIntegrate,
        variables.maxQ, elementList, element, x, y, z, elementParameters)
    iintradamp_Q = UtilityAnalysis.calc_iintradamp(iintra_Q, dampingFunction)
    Qiintradamp_Q = Q*iintradamp_Q
    rintra, Fintra_r = MainFunctions.calc_Fr(Q[Q<=variables.QmaxIntegrate], 
        Qiintradamp_Q[Q<=variables.QmaxIntegrate])
    
    # ---------------------Geometrical correction------------------------------
    
    thickness_sampling, phi_matrix = Geometry.check_phi_matrix(Q, phi_matrix_flag,
        ws1, ws2, r1, r2, d)
    
    absCorrFactor = IgorFunctions.absorption(Q)
    # two_theta = UtilityAnalysis.Qto2theta(Q)
    # absCorrFactor = Geometry.calcAbsCorrection(variables.abs_length, two_theta, 
        # variables.dac_thickness, 0.0)
    I_Q = I_Q/absCorrFactor
    Ibkg_Q = Ibkg_Q/absCorrFactor
    
    # ------------------------Starting minimization----------------------------

    scaleFactor = variables.scaleFactor
    density0 = variables.density
    
    # ----------------------First scale minimization---------------------------
    
    scaleStep = 0.05
    
    scaleFactor = Minimization.OptimizeScale(Q, I_Q, Ibkg_Q, J_Q, Iincoh_Q,
        fe_Q, variables.maxQ, variables.minQ, variables.QmaxIntegrate, Ztot,
        density0, scaleFactor, Sinf, variables.smoothingFactor, variables.rmin, 
        dampingFunction, Fintra_r, variables.iterations, scaleStep,
        sth, s0th, thickness_sampling, phi_matrix)
    
    # ----------------------First density minimization-------------------------
    
    densityStep = density0/50
    densityStepEnd = density0/250
    
    density = Minimization.OptimizeDensity(Q, I_Q, Ibkg_Q, J_Q, Iincoh_Q,
        fe_Q, variables.maxQ, variables.minQ, variables.QmaxIntegrate, Ztot, 
        density0, scaleFactor, Sinf, variables.smoothingFactor, variables.rmin, 
        dampingFunction, Fintra_r, variables.iterations, densityStep,
        densityStepEnd, sth, s0th, thickness_sampling, phi_matrix)
    
    print("density0, density", density0, density)
    numLoopIteration = 0
    
    while 1:
        if np.abs(density-density0) > density/25:
            print("First")
            scaleStep = 0.006
            densityStep = density/10
            WSamplestep=0.0008
            WRefstep=0.0008
        elif np.abs(density-density0) > density/75:
            print("Second")
            scaleStep = 0.0006
            densityStep = density/100
            WSamplestep=0.0002
            WRefstep=0.0002
        else:
            print("Third")
            scaleStep = 0.00006
            densityStep = density/1000
            WSamplestep=0.0001
            WRefstep=0.0001
        
        scaleFactor = Minimization.OptimizeScale(Q, I_Q, Ibkg_Q, J_Q, Iincoh_Q,
            fe_Q, variables.maxQ, variables.minQ, variables.QmaxIntegrate, Ztot,
            density, scaleFactor, Sinf, variables.smoothingFactor, variables.rmin, 
            dampingFunction, Fintra_r, variables.iterations, scaleStep)

        density0=density

        density = Minimization.OptimizeDensity(Q, I_Q, Ibkg_Q, J_Q, Iincoh_Q,
            fe_Q, variables.maxQ, variables.minQ, variables.QmaxIntegrate, Ztot, 
            density0, scaleFactor, Sinf, variables.smoothingFactor, variables.rmin, 
            dampingFunction, Fintra_r, variables.iterations, densityStep, density/250)
        
        numLoopIteration += 1
        print("numLoopIteration", numLoopIteration, scaleFactor, density)
        if (np.abs(density-density0) > np.abs(density/2500)) and (numLoopIteration <= 30):
           continue
        else:
            break
       
    print("final scale", scaleFactor, "final density", density)
    
    
    #scaleFactor = 1.0603638618
    #density = 29.9874611902
    
    
    Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scaleFactor, Ibkg_Q)
    alpha = MainFunctions.calc_alpha(J_Q[Q<=variables.QmaxIntegrate], Sinf, 
        Q[Q<=variables.QmaxIntegrate], Isample_Q[Q<=variables.QmaxIntegrate], 
        fe_Q[Q<=variables.QmaxIntegrate], Ztot, density)
    Icoh_Q = MainFunctions.calc_Icoh(alpha, Isample_Q, Iincoh_Q)

    S_Q = MainFunctions.calc_SQ(Icoh_Q, Ztot, fe_Q, Sinf, Q, variables.minQ, 
        variables.QmaxIntegrate, variables.maxQ)

    Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, 
        variables.smoothingFactor, variables.minQ, variables.QmaxIntegrate, variables.maxQ)

    SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Sinf,
        dampingFunction)

    i_Q = MainFunctions.calc_iQ(SsmoothDamp_Q, Sinf)
    
    Qi_Q = Q*i_Q
    r, F_r = MainFunctions.calc_Fr(Q[Q<=variables.QmaxIntegrate], 
        Qi_Q[Q<=variables.QmaxIntegrate])
    Fopt_r, deltaFopt_r = Optimization.calc_optimize_Fr(variables.iterations, F_r,
                Fintra_r, density, i_Q[Q<=variables.QmaxIntegrate], Q[Q<=variables.QmaxIntegrate],
                Sinf, J_Q[Q<=variables.QmaxIntegrate], r, variables.rmin, "n")
    
    Scorr_Q = MainFunctions.calc_SQCorr(Fopt_r, r, Q, Sinf)
    
    Utility.plot_data(Q, SsmoothDamp_Q, "S_Q", "Q", "S_Q", "S(Q)", "y")
    Utility.plot_data(Q, Scorr_Q, "S_Q", "Q", "S_Q", "Scorr(Q)", "y")
    Utility.plot_data(r, F_r, "F_r", "r", "F_r", "F(r)", "y")
    Utility.plot_data(r, Fopt_r, "F_r", "r", "F_r", "Fopt(r)", "y")
    
    Utility.write_file("./S_Q.txt", Q, SsmoothDamp_Q)
    Utility.write_file("./Scorr_Q.txt", Q, Scorr_Q)
    Utility.write_file("./F_r.txt", r, F_r)
    Utility.write_file("./Fopt_r.txt", r, Fopt_r)
    
    plt.show()