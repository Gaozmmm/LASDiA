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

"""Set of modules used in LASDiA to calculate the diamond absorption and Soller
Slits correction.

The nomenclature and the procedure follow the articles:
Weck et al. Rev. Sci. Instrum. 84, 063901 (2013)
Yaoita et al. Rev. Sci. Instrum. 68, 2106 (1997)
Eggert et al. 2002 PRB, 65, 174105.

For the functions arguments and the returns I followed this convetion for the notes:
argument: description - type
return: description - type.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by
an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""


import numpy as np
from scipy.integrate import simps

from modules import Utility
from modules import UtilityAnalysis


def calcAbsCorrection(abs_length, two_theta, thickness, angle):
    """Function to calculate the absorption correction.
    This function can be used to calculate the absorption correction for the diamond
    or for any other object between the sample and the detector.

    The characteristics for some diamonds can be found here:
    http://www.almax-easylab.com/DiamondSelectionPage.aspx

    Parameters
    ----------
    abs_length  : float
                  absorption length (cm), @33keV 1.208cm
    two_theta   : numpy array
                  diffraction angle (rad)
    thickness   : float
                  object thickness (cm)
    angle       : float
                  object rotation angle respect the XRay beam (deg)

    Returns
    -------
    corr_factor : numpy array
                  correction factor
    """

    # for now...
    # wavelenght  : float
    #               XRay beam wavelenght (nm), @ESRF ID27 0.03738nm
    # wavelenght = 0.03738 # nm
    # two_theta = Qto2theta(Q) # rad

    mu_l = 1/abs_length
    angle = np.radians(angle)

    path_lenght = thickness / np.cos(two_theta - angle)

    corr_factor = np.exp(-mu_l * path_lenght)

    # I_Qeff = I_Q / corr_factor

    return corr_factor


def calc_phi_angle(ws1, ws2, r1, r2, d, two_theta, xth):
    """Function to calculate the dispersion angle for the Soller Slits correction.

    Parameters
    ----------
    ws1     : float
              width of the inner slit (cm)
    ws2     : float
              width of the outer slit (cm)
    r1      : float
              curvature radius of first slit (cm)
    r2      : float
              curvature radius of second slit (cm)
    d       : float
              slit thickness (cm)
    two_theta : float
              diffraction angle (rad)
    xth     : float
              i-th point position on x-axis (cm)

    Returns
    -------
    phi     : float
              dispersion angle (rad)
    """

    gamma_2 = np.arcsin(ws1/(2*r1))

    alpha1 = np.arctan( r1 * np.sin(two_theta + gamma_2) / (r1*np.cos(two_theta + gamma_2) - xth ))
    alpha2 = np.arctan( r1 * np.sin(two_theta - gamma_2) / (r1*np.cos(two_theta - gamma_2) - xth ))

    beta1 = np.arctan( (r2+d) * np.sin(two_theta + gamma_2) / ((r2+d)*np.cos(two_theta + gamma_2) - xth ))
    beta2 = np.arctan( (r2+d) * np.sin(two_theta - gamma_2) / ((r2+d)*np.cos(two_theta - gamma_2) - xth ))

    psi = min(alpha1, beta1) - max(alpha2, beta2)
    phi = max(0, psi)

    return phi


def calc_phi_matrix(thickness, two_theta, ws1, ws2, r1, r2, d, num_point):
    """Function to calculate the dispersion angle matrix.
    half_thick in cm

    Parameters
    ----------
    thickness          : float
                         object thickness (sample or sample+DAC) (cm)
    two_theta          : numpy array
                         diffraction angle (rad)
    ws1                : float
                         width of the inner slit (cm)
    ws2                : float
                         width of the outer slit (cm)
    r1                 : float
                         curvature radius of first slit (cm)
    r2                 : float
                         curvature radius of second slit (cm)
    d                  : float
                         slit thickness (cm)
    num_point          : int
                         number of point for the thickness array

    Returns
    -------
    thickness_sampling : numpy array
                         array with the thickness values (cm)
    phi_matrix         : 2D numpy array
                         dispersion angle matrix (rad)
    """

    # thickness_sampling = np.linspace(-thickness/2, thickness/2, num=num_point) # num=500)
    thickness_sampling = np.linspace(0, thickness, num=num_point) # num=500)
    phi_matrix = np.zeros((thickness_sampling.size, two_theta.size))

    for i, val_sth in enumerate(thickness_sampling):
        for j, val_2theta in enumerate(two_theta):
            phi_matrix[i][j] = calc_phi_angle(ws1, ws2, r1, r2, d, two_theta[j], \
            thickness_sampling[i])

    return (thickness_sampling, phi_matrix)


def calc_T_MCC(sample_thickness, thickness_sampling, phi_matrix, norm):
    """Function to calculate the MCC transfer function for the sample, the DAC 
        and sample+DAC (W. eq. 10, 11).

    Parameters
    ----------
    sample_thickness   : float
                         sample thickness
    thickness_sampling : numpy array
                         array with the thickness values for sample+DAC
    phi_matrix         : 2D numpy array
                         dispersion angle matrix for sample+DAC (rad)
    norm               : string
                         flag to normalize the MCC transfer function to start to 1

    Returns
    -------
    T_MCC_sample       : numpy array
                         MCC sample transfer function
    T_MCC_DAC          : numpy array
                         MCC DAC transfer function
    T_MCC_ALL          : numpy array
                         MCC sample+DAC transfer function
    """

    # mask = (thickness_sampling >= -sample_thickness/2) & (thickness_sampling <= sample_thickness/2)
    mask = (thickness_sampling >= 0) & (thickness_sampling <= sample_thickness/2)

    T_MCC_ALL = simps(phi_matrix, axis=0, even="first")
    T_MCC_sample = simps(phi_matrix[mask], axis=0, even="first")
    T_MCC_DAC = T_MCC_ALL - T_MCC_sample

    if norm.lower() == "y":
        T_MCC_sample /= T_MCC_sample[0]
        T_MCC_DAC /= T_MCC_DAC[0]
        T_MCC_ALL /= T_MCC_ALL[0]

    return (T_MCC_sample, T_MCC_DAC, T_MCC_ALL)


def calc_T_DAC_MCC_bkg_corr(T_DAC_MCC_sth, T_DAC_MCC_s0th):
    """Function to calculate the background correction factor for the MCC (W. eq. 12).

    Parameters
    ----------
    T_DAC_MCC_sth  : numpy array
                     MCC DAC transfer function
    T_DAC_MCC_s0th : numpy array
                     MCC DAC transfer function with sample thickness for the 
                     reference spectra
    
    Returns
    -------
    corr_factor    : numpy array
                     the background correction factor
    """

    corr_factor = T_DAC_MCC_sth / T_DAC_MCC_s0th

    return corr_factor


def calc_empty_cell_bkg(Q, Ibkg_Q, diamond_abs_corr_factor, Iincoh_Q, Ztot, fe_Q, mu2):
    """Function to calculate the empty-cell background from the solid-sample 
       background (E. eq. 32).
    
    Parameters
    ----------
    Q                       : numpy array
                              momentum transfer (nm^-1)
    Ibkg_Q                  : numpy array
                              background scattering intensity
    diamond_abs_corr_factor : numpy array
                              diamond correction factor
    Iincoh_Q                : numpy array
                              incoherent scattering intensity
    Ztot                    : int
                              total Z number
    fe_Q                    : numpy array
                              effective electric form factor
    mu2                     : float ???
                              mean-square component of the displacement of the
                              molecular units in the scattering direction
    
    Returns
    -------
    Ibkg_Q                  : numpy array
                              empty-cell background scattering intensity
    """
    
    FermiWidth=0
    FermiCutoff=70
    TDS=(1-1/(1+np.exp(FermiWidth*(Q-FermiCutoff)))) * Ztot**2*fe_Q**2 * (1-np.exp(-mu2*Q^2))
    
    # DiffuseIntensity = N/alpha'
    # Bkgd=BkgdSolid-DiffuseIntensity*AbsRefCalc*(Iincoh+TDS)
    
    Ibkg_Q = Ibkg_Q - diamond_abs_corr_factor*(Iincoh_Q - TDS)
    
    return Ibkg_Q


def geometry_correction(Q, I_Q, Qbkg, Ibkg_Q, variables, phi_matrix_flag):
    """Function to calcultate all intensity geometrical corrections.
    
    Parameters
    ----------
    Q               : numpy array
                      momentum transfer (nm^-1)
    I_Q             : numpy array
                      measured scattering intensity
    Qbkg            : numpy array
                      background momentum transfer (nm^-1)
    Ibkg_Q          : numpy array
                      background scattering intensity
    variables       : module
                      input variables setted by the user
    phi_matrix_flag : string
                      flag for the phi matrix calculation:
                      "y": calculate phi matrix and save on file
                      "n": read the phi matrix from file
    
    Returns
    -------
    I_Q             : numpy array
                      corrected measured sample intensity
    Ibkg_Q          : numpy array
                      corrected background intensity
    """
    
    two_theta = UtilityAnalysis.Qto2theta(Q)
    
    abs_corr_factor = calcAbsCorrection(variables.abs_length, \
        two_theta, variables.dac_thickness, 0)
    
    num_point = 1000
    
    phi_matrix_path = "./phi_matrix_" + variables.molecule
    
    if phi_matrix_flag.lower() == "y": 
        ws1, ws2, r1, r2, d = Utility.read_MCC_file(variables.MCC_path, variables.MCC_type)
        thickness_sampling, phi_matrix = calc_phi_matrix(variables.phi_matrix_thickness, \
            two_theta, ws1, ws2, r1, r2, d, num_point)
        np.save(phi_matrix_path, phi_matrix)
    else:
        thickness_sampling = np.linspace(0, variables.phi_matrix_thickness, num=num_point)
        phi_matrix = np.load(phi_matrix_path + ".npy")
    
    T_MCC_sample3, T_MCC_DAC3, T_MCC_ALL3 = calc_T_MCC(0.003, thickness_sampling, \
        phi_matrix, "y")
    T_MCC_sample4, T_MCC_DAC4, T_MCC_ALL4 = calc_T_MCC(0.004, thickness_sampling, \
        phi_matrix, "y")
    T_MCC_corr_factor_bkg = calc_T_DAC_MCC_bkg_corr(T_MCC_DAC4, T_MCC_DAC3)
    
    I_Q = I_Q /(abs_corr_factor * T_MCC_sample4)
    Ibkg_Q  = Ibkg_Q * T_MCC_corr_factor_bkg / (abs_corr_factor * T_MCC_sample4)
    
    return (I_Q, Ibkg_Q)


def MCCCorrection(Q, I_Q, Qbkg, Ibkg_Q, variables, phi_matrix_flag):
    """Function to calcultate all intensity geometrical corrections.
    
    Parameters
    ----------
    Q               : numpy array
                      momentum transfer (nm^-1)
    I_Q             : numpy array
                      measured scattering intensity
    Qbkg            : numpy array
                      background momentum transfer (nm^-1)
    Ibkg_Q          : numpy array
                      background scattering intensity
    variables       : module
                      input variables setted by the user
    phi_matrix_flag : string
                      flag for the phi matrix calculation:
                      "y": calculate phi matrix and save on file
                      "n": read the phi matrix from file
    
    Returns
    -------
    I_Q             : numpy array
                      corrected measured sample intensity
    Ibkg_Q          : numpy array
                      corrected background intensity
    """
    
    two_theta = UtilityAnalysis.Qto2theta(Q)
    
    num_point = 1000
    
    phi_matrix_path = "./phi_matrix_" + variables.molecule
    
    if phi_matrix_flag.lower() == "y": 
        ws1, ws2, r1, r2, d = Utility.read_MCC_file(variables.MCC_path, variables.MCC_type)
        thickness_sampling, phi_matrix = calc_phi_matrix(variables.phi_matrix_thickness, \
            two_theta, ws1, ws2, r1, r2, d, num_point)
        np.save(phi_matrix_path, phi_matrix)
    else:
        thickness_sampling = np.linspace(0, variables.phi_matrix_thickness, num=num_point)
        phi_matrix = np.load(phi_matrix_path + ".npy")
    
    T_MCC_sample3, T_MCC_DAC3, T_MCC_ALL3 = calc_T_MCC(0.007, thickness_sampling, \
        phi_matrix, "y") # s0th
    T_MCC_sample4, T_MCC_DAC4, T_MCC_ALL4 = calc_T_MCC(0.01, thickness_sampling, \
        phi_matrix, "y") #sth
    T_MCC_corr_factor_bkg = calc_T_DAC_MCC_bkg_corr(T_MCC_DAC4, T_MCC_DAC3)
    
    I_Q = I_Q /T_MCC_sample4
    Ibkg_Q  = Ibkg_Q * T_MCC_corr_factor_bkg / (T_MCC_sample4)
    
    return (I_Q, Ibkg_Q)