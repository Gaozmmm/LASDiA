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

"""Set of modules used in LASDiA to calculate the relations among the different set of Structure Factors.

The formalisms used are Ashcroft-Langreth and Faber-Ziman.

The nomenclature and the procedure follow the book:
Waseda - The Structure of Non-Crystalline Materials.

For the functions arguments and the returns I followed this convetion for the notes:
arguments: description - type
returns: description - type.

For the variables name I used this convention:
if the variable symbolizes a function, its argument is preceded by an underscore: f(x) -> f_x
otherwise it is just the name.
"""

import matplotlib.pyplot as plt

import sys
import os

import numpy as np
import scipy.constants as sc
from scipy import fftpack
from scipy.integrate import simps
from scipy.stats import chisquare

from modules.LiquidStructure import *

import math

def ALtoFZ(N1, N2, S11, S12, S22):
    """Function to transform the AL to FZ (eq 1.2.17)
    
    arguments:
    N1:
    N2:
    S11:
    S12:
    S22:
    
    returns:
    a11:
    a12:
    a22:
    """
    
    c1 = N1/(N1+N2)
    c2 = N2/(N1+N2)
    
    a11_Q = (S11 - c2)/c1
    a12_Q = S12/math.sqrt(c1*c2) + 1
    a22_Q = (S22 - c1)/c2
    
    return (a11_Q, a12_Q, a22_Q)
    
    
def FZtoAL(N1, N2, a11, a12, a22):
    """Function to transform the AL to FZ (eq 1.2.17)
    
    arguments:
    N1:
    N2:
    a11:
    a12:
    a22:
    
    returns:
    S11:
    S12:
    S22:
    """
    
    c1 = N1/(N1+N2)
    c2 = N2/(N1+N2)
    
    S11_Q = 1 + c1*(a11_Q - 1)
    S12_Q = math.sqrt(c1*c2) * (a12_Q - 1)
    S11_Q = 1 + c2*(a22_Q - 1)
    
    return (S11_Q, S12_Q, S22_Q)
    
    
def FaberZiman():
    """
    
    """