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


def open_file(default_dir=None):
    """Function to load the input file with a pop-up windows for the script
    version.
    """

    from sys import argv
    from os import environ
    from PyQt5.QtWidgets import QApplication, QFileDialog
    if default_dir == None:
        default_dir = environ["HOME"]
    app = QApplication(argv)
    filename, _ = QFileDialog.getOpenFileName(caption="Load Input File",
        directory=default_dir, filter="*txt")
    app.exit()
    return filename

if __name__ == "__main__":
    
    #---------------------------Files reading----------------------------------

    # inputFile_path = open_file("./")

    inputVariables = Utility.read_inputFile("./inputFile.txt")

    elementList = Utility.molToElemList(inputVariables["molecule"])
    elementParameters = Utility.read_parameters(elementList, inputVariables["elementParamsPath"])
    elementPosition = Utility.read_xyz_file(inputVariables["xyzPath"])
    
    Q, I_Q = Utility.read_file(inputVariables["dataFile"])
    Qbkg, Ibkg_Q  = Utility.read_file(inputVariables["bkgFile"])

    # plt.plot(Q, I_Q)
    # plt.plot(Qbkg, Ibkg_Q)
    # plt.show

    #--------------------Preliminary calculation-------------------------------

    Q, I_Q = UtilityAnalysis.data_interpolation(Q, I_Q, inputVariables["minQ"],
        inputVariables["maxQ"], inputVariables["numPoints"])
    Qbkg, Ibkg_Q = UtilityAnalysis.data_interpolation(Qbkg, Ibkg_Q, inputVariables["minQ"],
        inputVariables["maxQ"], inputVariables["numPoints"])

    # plt.plot(Q, I_Q)
    # plt.plot(Qbkg, Ibkg_Q)
    # plt.show

    fe_Q, Ztot = MainFunctions.calc_eeff(elementList, Q, elementParameters)
    Iincoh_Q = MainFunctions.calc_Iincoh(elementList, Q, elementParameters)
    J_Q = MainFunctions.calc_JQ(Iincoh_Q, Ztot, fe_Q)
    Sinf = MainFunctions.calc_Sinf(elementList, fe_Q, Q, Ztot, elementParameters)
    dampingFunction = UtilityAnalysis.calc_dampingFunction(Q, inputVariables["dampingFactor"],
        inputVariables["QmaxIntegrate"], inputVariables["typeFunction"])
        
    #-------------------Intra-molecular components-----------------------------

    iintra_Q = Optimization.calc_iintra(Q, fe_Q, Ztot, inputVariables["QmaxIntegrate"],
        inputVariables["maxQ"], elementList, elementPosition["element"],
        elementPosition["x"], elementPosition["y"], elementPosition["z"], elementParameters)
    iintradamp_Q = UtilityAnalysis.calc_iintradamp(iintra_Q, dampingFunction)
    Qiintradamp_Q = Q*iintradamp_Q
    rintra, Fintra_r = MainFunctions.calc_Fr(Q[Q<=inputVariables["QmaxIntegrate"]], 
        Qiintradamp_Q[Q<=inputVariables["QmaxIntegrate"]])
    
    # ---------------------Geometrical correction------------------------------

    # global phi_matrix
    # global thickness_sampling

    if inputVariables["mccFlag"].lower() == "y":
        print("Start matrix calculation or reading")
        thickness_sampling, phi_matrix = Geometry.check_phi_matrix(Q, inputVariables["ws1"],
            inputVariables["ws2"], inputVariables["r1"], 
            inputVariables["r2"], inputVariables["d"], inputVariables["phiMatrixCalcFlag"],
            inputVariables["phiMatrixPath"])
        print("End matrix calculation or reading")
    else:
        phi_matrix = 0.0
        thickness_sampling = 0.0

    absCorrFactor = Geometry.calc_abs_correction(Q, inputVariables["absLength"], 
        inputVariables["dacThickness"], 0.0)
    I_Q = I_Q/absCorrFactor
    Ibkg_Q = Ibkg_Q/absCorrFactor

    # plt.plot(Q, I_Q)
    # plt.plot(Qbkg, Ibkg_Q)
    # plt.show()
    
    # ------------------------Starting minimization----------------------------

    scaleFactor = inputVariables["scaleFactor"]
    density0 = inputVariables["density"]
    
    # ----------------------First scale minimization---------------------------
    
    scaleStep = 0.05

    print("Start first scale minimization")
    scaleFactor = Minimization.OptimizeScale(Q, I_Q, Ibkg_Q, J_Q, Iincoh_Q,
        fe_Q, inputVariables["minQ"], inputVariables["QmaxIntegrate"], inputVariables["maxQ"],
        Ztot, density0, scaleFactor, Sinf, inputVariables["smoothingFactor"],
        inputVariables["rmin"], dampingFunction, Fintra_r, inputVariables["iterations"],
        scaleStep, inputVariables["sth"], inputVariables["s0th"], inputVariables["mccFlag"],
        thickness_sampling, phi_matrix)
    print("End first scale minimization")
    
    # ----------------------First density minimization-------------------------
    
    densityStep = density0/50
    densityStepEnd = density0/250
    
    print("Start first density minimization")
    density = Minimization.OptimizeDensity(Q, I_Q, Ibkg_Q, J_Q, Iincoh_Q,
        fe_Q, inputVariables["minQ"], inputVariables["QmaxIntegrate"], inputVariables["maxQ"],
        Ztot, density0, scaleFactor, Sinf, inputVariables["smoothingFactor"],
        inputVariables["rmin"], dampingFunction, Fintra_r, inputVariables["iterations"],
        densityStep, inputVariables["sth"], inputVariables["s0th"], inputVariables["mccFlag"],
        thickness_sampling, phi_matrix)
    print("End first density minimization")

    # --------------------Free parameters minimization-------------------------

    # print("density0, density", density0, density)
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
        
        print("Start scale minimization")
        scaleFactor = Minimization.OptimizeScale(Q, I_Q, Ibkg_Q, J_Q, Iincoh_Q,
            fe_Q, inputVariables["minQ"], inputVariables["QmaxIntegrate"], inputVariables["maxQ"],
            Ztot, density, scaleFactor, Sinf, inputVariables["smoothingFactor"],
            inputVariables["rmin"], dampingFunction, Fintra_r, inputVariables["iterations"],
            scaleStep, inputVariables["sth"], inputVariables["s0th"], inputVariables["mccFlag"],
            thickness_sampling, phi_matrix)
        print("End scale minimization")

        density0=density

        print("Start density minimization")
        density = Minimization.OptimizeDensity(Q, I_Q, Ibkg_Q, J_Q, Iincoh_Q,
            fe_Q, inputVariables["minQ"], inputVariables["QmaxIntegrate"], inputVariables["maxQ"],
            Ztot, density0, scaleFactor, Sinf, inputVariables["smoothingFactor"],
            inputVariables["rmin"], dampingFunction, Fintra_r, inputVariables["iterations"],
            densityStep, inputVariables["sth"], inputVariables["s0th"], inputVariables["mccFlag"],
            thickness_sampling, phi_matrix)
        print("End density minimization")
        
        numLoopIteration += 1
        print("numLoopIteration", numLoopIteration, scaleFactor, density)
        if (np.abs(density-density0) > np.abs(density/2500)) and (numLoopIteration <= 30):
           continue
        else:
            break
       
    print("final scale", scaleFactor, "final density", density)
    
    Isample_Q = MainFunctions.calc_IsampleQ(I_Q, scaleFactor, Ibkg_Q)
    alpha = MainFunctions.calc_alpha(J_Q[Q<=inputVariables["QmaxIntegrate"]], Sinf, 
        Q[Q<=inputVariables["QmaxIntegrate"]], Isample_Q[Q<=inputVariables["QmaxIntegrate"]], 
        fe_Q[Q<=inputVariables["QmaxIntegrate"]], Ztot, density)
    Icoh_Q = MainFunctions.calc_Icoh(alpha, Isample_Q, Iincoh_Q)

    S_Q = MainFunctions.calc_SQ(Icoh_Q, Ztot, fe_Q, Sinf, Q, inputVariables["minQ"], 
        inputVariables["QmaxIntegrate"], inputVariables["maxQ"])

    Ssmooth_Q = UtilityAnalysis.calc_SQsmoothing(Q, S_Q, Sinf, 
        inputVariables["smoothingFactor"], inputVariables["minQ"], inputVariables["QmaxIntegrate"], inputVariables["maxQ"])

    SsmoothDamp_Q = UtilityAnalysis.calc_SQdamp(Ssmooth_Q, Sinf,
        dampingFunction)

    i_Q = MainFunctions.calc_iQ(SsmoothDamp_Q, Sinf)
    
    Qi_Q = Q*i_Q
    r, F_r = MainFunctions.calc_Fr(Q[Q<=inputVariables["QmaxIntegrate"]], 
        Qi_Q[Q<=inputVariables["QmaxIntegrate"]])
    Fopt_r, deltaFopt_r = Optimization.calc_optimize_Fr(inputVariables["iterations"], F_r,
                Fintra_r, density, i_Q[Q<=inputVariables["QmaxIntegrate"]], Q[Q<=inputVariables["QmaxIntegrate"]],
                Sinf, J_Q[Q<=inputVariables["QmaxIntegrate"]], r, inputVariables["rmin"], "n")
    
    Scorr_Q = MainFunctions.calc_SQCorr(Fopt_r, r, Q, Sinf)
    
    Utility.plot_data(Q, SsmoothDamp_Q, "S_Q", "Q", "S_Q", "S(Q)", "y")
    Utility.plot_data(Q, Scorr_Q, "S_Q", "Q", "S_Q", "Scorr(Q)", "y")
    Utility.plot_data(r, F_r, "F_r", "r", "F_r", "F(r)", "y")
    Utility.plot_data(r, Fopt_r, "F_r", "r", "F_r", "Fopt(r)", "y")
    
    # Utility.write_file("./S_Q_Ar.txt", Q, SsmoothDamp_Q)
    # Utility.write_file("./Scorr_Q_Ar.txt", Q, Scorr_Q)
    # Utility.write_file("./F_r_Ar.txt", r, F_r)
    # Utility.write_file("./Fopt_r_Ar.txt", r, Fopt_r)
    
    plt.show()