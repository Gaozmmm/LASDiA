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

"""Set of modules used in the Liquid program

The nomenclature and the procedures follow the article:
Eggert et al. 2002 PRB, 65, 174105

For the functions arguments and the returns I followed this convetion for the notes:
arguments: description - type
returns: description - type

For the variables name I used this convention:
if the variable symbolizes a function, its argument is preceded by an underscore: f(x) -> f_x
"""

import matplotlib.pyplot as plt

import sys
import os

import numpy as np
import scipy.constants as sc
from scipy import fftpack
from scipy.integrate import simps


def IterateDeltaF(Q, SQ_SS, SQCorr):
    """
    
    """
    
    Qmax = np.amax(Q)
    Damping = 0.5
    argQmax = np.aargmax(Q)
    
    QiQ_SS = exp(-(Damping/Qmax**2)*x^2)*Q*(SQ_SS-SQ_SS(argQmax+1))
    QiQCorr = Q*(SQCorr-SQCorr(argQmax+1))
    Del_iQ = (QiQ_SS-QiQCorr)/(Q*SQCorr)
    
    return (QiQ_SS, QiQCorr, Del_iQ)
    
def FitRemoveGofRPeaks():
    """
    arguments:
    
    returns:    
    """

    IFFT_GofR(DeltaG,"QiQCorrection",1)	
    QiQ1[1,x2pnt(QiQ1,QMaxIntegrate)-1]=QiQ1-(QiQ1/(x*(SQinf+JofQ/Ztot^2))+1)*QiQCorr
    FFT_QIofQ("QiofQ1", "GofRCorrection")
    DelG=0
    DelG[,x2pnt(DelG,Rnn)]=GR1-GRIdealSmallR
    Wavestats/Q/R=(0.95*Rnn,) GR1
    Rnn=0.99*V_minloc

    

# def OptimizeScaleGofRCorr():

# def OptimizeDensityGofRCorr():

# def OptimizeWsampleGofRCorr():

# def OptimizeWrefGofRCorr():
