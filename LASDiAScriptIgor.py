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

import matplotlib.pyplot as plt
import numpy as np
# from itertools import product
# from timeit import default_timer as timer
#from scipy.integrate import simps

#from modules import Formalism
#from modules import Geometry
from modules import IgorFunctions
#from modules import KaplowMethod
from modules import MainFunctions
from modules import Minimization
from modules import Optimization
from modules import Utility
from modules import UtilityAnalysis

#from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget


if __name__ == "__main__":
    
    #---------------------------Files reading----------------------------------
    
    variables = Utility.read_inputFile("./inputFile.txt")
    
    elementList = Utility.molToelemList(variables.molecule)
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
    J_Q = IgorFunctions.calc_JQ(Iincoh_Q, fe_Q)
    Sinf = MainFunctions.calc_Sinf(elementList, fe_Q, Q, Ztot, elementParameters)
    
    dampingFunction = UtilityAnalysis.calc_dampingFunction(Q, variables.dampingFactor,
        variables.QmaxIntegrate, variables.typeFunction)

    #-------------------Intra-molecular components-----------------------------

    iintra_Q = Optimization.calc_iintra(Q, fe_Q, Ztot, variables.QmaxIntegrate, 
        variables.maxQ, elementList, element, x, y, z, elementParameters)
    iintradamp_Q = UtilityAnalysis.calc_iintradamp(iintra_Q, Q, variables.QmaxIntegrate, 
        dampingFunction)
    rintra, Fintra_r = IgorFunctions.calc_FFT_QiQ(Q, Q*iintradamp_Q, variables.QmaxIntegrate)
    
    _, Fintra_r = UtilityAnalysis.rebinning(rintra, Fintra_r, np.amin(rintra), 
        np.amax(rintra), 8192)
    
    _, dampingFunction = UtilityAnalysis.rebinning(Q, dampingFunction, 0.0, 
        variables.maxQ, variables.NumPoints)
    
    # ---------------------Geometrical correction------------------------------
    
    absCorrFactor = IgorFunctions.absorption(Q)
    I_Q = I_Q/absCorrFactor
    Ibkg_Q = Ibkg_Q/absCorrFactor
    
    # ------------------------Starting minimization----------------------------
    # To keep in mind:
    # gi_1 = density0
    # gi = density
    # Init1 = density

    scaleFactor = variables.scaleFactor
    density0 = variables.density
    
    # ----------------------First scale minimization---------------------------
    
    scaleStep = 0.05
    # scaleStepEnd = 0.00006
    
    scaleFactor = IgorFunctions.OptimizeScaleGofRCorr(Q, I_Q, Ibkg_Q,
        J_Q, fe_Q, variables.maxQ, variables.minQ, variables.QmaxIntegrate, Ztot,
        density0, scaleFactor, Sinf, variables.smoothingFactor, variables.rmin, 
        variables.NumPoints, dampingFunction, Fintra_r, variables.iterations,
        scaleStep)
    
    # ----------------------First density minimization-------------------------
    
    densityStep = density0/50
    densityStepEnd = density0/250
    
    density = IgorFunctions.OptimizeDensityGofRCorr(Q, I_Q, Ibkg_Q,
        J_Q, fe_Q, variables.maxQ, variables.minQ,
        variables.QmaxIntegrate, Ztot, density0, scaleFactor, Sinf,
        variables.smoothingFactor, variables.rmin, variables.NumPoints,
        dampingFunction, Fintra_r, variables.iterations, densityStep, densityStepEnd)
    
    print("density0, density", density0, density)
    numLoopIteration = 0
    
    while (np.abs(density-density0) > np.abs(density/2500)) and (numLoopIteration <= 30):
        if np.abs(density-density0) > density/25:
            scaleStep = 0.006
            densityStep = density/10
        elif np.abs(density-density0) > density/75:
            scaleStep = 0.0006
            densityStep = density/100
        else:
            scaleStep = 0.00006
            densityStep = density/1000
        
        scaleFactor = IgorFunctions.OptimizeScaleGofRCorr(Q, I_Q, Ibkg_Q, 
            J_Q, fe_Q, variables.maxQ, variables.minQ,
            variables.QmaxIntegrate, Ztot, density, scaleFactor, Sinf, 
            variables.smoothingFactor, variables.rmin, variables.NumPoints,
            dampingFunction, Fintra_r, variables.iterations, scaleStep)
        
        density0=density
        density = IgorFunctions.OptimizeDensityGofRCorr(Q, I_Q, Ibkg_Q,
            J_Q, fe_Q, variables.maxQ, variables.minQ,
            variables.QmaxIntegrate, Ztot, density0, scaleFactor, Sinf,
            variables.smoothingFactor, variables.rmin, variables.NumPoints,
            dampingFunction, Fintra_r, variables.iterations, densityStep,
            density/250)
        
        numLoopIteration += 1
        print("numLoopIteration", numLoopIteration, scaleFactor, density)
       
    print("final scale", scaleFactor, "final density", density)
    plt.show()