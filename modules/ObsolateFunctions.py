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

"""Module containing the obsolate functions used in LASDiA.
These functions were used to test some part of the code or they were substituted by a new version.

The nomenclature and the procedure follow the article:
Eggert et al. 2002 PRB, 65, 174105.

For the functions arguments and the returns I followed this convetion for the notes:
arguments: description - type
returns: description - type.

For the variables name I used this convention:
if the variable symbolizes a mathematical function, its argument is preceded by an underscore: f(x) -> f_x
otherwise it is symbolized with just its name.
"""


def calc_aff(element, Q, aff_path):
    """Function to calculate the Atomic Form Factor.
    The atomic form factor is calculated with the formula from:
    http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php

    and the parameters from the article:
    Hajdu Acta Cryst. (1972). A28, 250

    Parameters
    ----------
    element  : string
               chemical element
    Q        : numpy array
               momentum transfer (nm^-1)
    aff_path : string
               path of atomic scattering form factor parameters
    
    Returns
    -------
    f_Q      : numpy array
               atomic form factor
    """

    # open, read and close the parameters file
    # the first line contain the header and it is useless for the calculation
    # file = open("./affParamCEA.txt", "r")
    file = open(aff_path, "r")
    header1 = file.readline()
    lines = file.readlines()
    file.close()

    # scan the lines and when it find the right element, save the parameters in variables
    for line in lines:
        columns = line.split()
        if columns[0] == element:
            a1 = float(columns[1])
            b1 = float(columns[2])
            a2 = float(columns[3])
            b2 = float(columns[4])
            a3 = float(columns[5])
            b3 = float(columns[6])
            a4 = float(columns[7])
            b4 = float(columns[8])
            c = float(columns[9])
            break

    # Calculate the atomic form factor as:
    # f(Q) = f1(Q) + f2(Q) + f3(Q) + f4(Q) + c
    # fi(Q) = ai * exp(-bi * (Q/4pi)^2)

    f1_Q = a1 * np.exp(-b1 * (Q/(4*10*np.pi))**2)
    f2_Q = a2 * np.exp(-b2 * (Q/(4*10*np.pi))**2)
    f3_Q = a3 * np.exp(-b3 * (Q/(4*10*np.pi))**2)
    f4_Q = a4 * np.exp(-b4 * (Q/(4*10*np.pi))**2)

    f_Q = f1_Q + f2_Q + f3_Q + f4_Q + c

    return f_Q


def calc_eeff(elementList, Q, incoh_path, aff_path):
    """Function to calculate the effective electron Form Factor, fe (eq. 10).

    Parameters
    ----------
    elementList  : dictionary("element": multiplicity)
                   chemical elements of the sample with their multiplicity
                   element      : string
                                  chemical element
                   multiplicity : int
                                  chemical element multiplicity
    Q            : numpy array
                   momentum transfer (nm^-1)
    incoh_path   : string
                   path of incoherent scattered intensities parameters
    aff_path     : string
                   path of atomic scattering form factor parameters
    
    Returns
    -------
    fe_Q         : numpy array
                   effective electric form factor
    Ztot         : int
                   total Z number
    """

    fe_Q = 0
    Ztot = 0

    # file = open("./incohParamCEA.txt", "r")
    file = open(incoh_path, "r")
    header1 = file.readline()
    lines = file.readlines()
    file.close()

    # scan the lines and when it find the right element take the Z number
    for element, multiplicity in elementList.items():
        # print (element, multiplicity)
        for line in lines:
            columns = line.split()
            if columns[0] == element:
                Ztot += multiplicity * float(columns[1])
                break

    # print (Ztot)

    for element, multiplicity in elementList.items():
        fe_Q += multiplicity * calc_aff(element, Q, aff_path)

    fe_Q /= Ztot

    return (fe_Q, Ztot)


def calc_Kp(fe_Q, element, Q, aff_path):
    """Function to calculate the average of effective atomic number Kp (eq. 11, 14).

    Parameters
    ----------
    fe_Q    : numpy array
              effective electric form factor
    element : string
              chemical element of the sample
    Q       : numpy array
              momentum transfer (nm^-1)
    
    Returns
    -------
    Kp      : float
              average of effective atomic number
    """

    # effective atomic number
    Kp_Q = calc_aff(element, Q, aff_path)/fe_Q

    # average effective atomic number
    Kp = np.mean(Kp_Q)

    return Kp


def calc_Sinf(elementList, fe_Q, Q, Ztot, aff_path):
    """Function to calculate Sinf (eq. 19).

    Parameters
    ----------
    elementList  : dictionary("element": multiplicity)
                   chemical elements of the sample with their multiplicity
                   element      : string
                                  chemical element
                   multiplicity : int
                                  chemical element multiplicity
    fe_Q         : numpy array
                   effective electric form factor
    Q            : numpy array
                   momentum transfer (nm^-1)
    Ztot         : int
                   total Z number
    
    Returns
    -------
    Sinf         : float
                   Sinf
    """
    
    sum_Kp2 = 0
    
    for element, multiplicity in elementList.items():
        sum_Kp2 += multiplicity * calc_Kp(fe_Q, element, Q, aff_path)**2

    Sinf = sum_Kp2 / Ztot**2

    return Sinf


def calc_Iincoh(elementList, Q, incoh_path, aff_path):
    """Function to calculate the incoherent scattering intensity Iincoh(Q).
    The incoherent scattering intensity is calculated with the formula from the article:
    Hajdu Acta Cryst. (1972). A28, 250

    Parameters
    ----------
    elementList  : dictionary("element": multiplicity)
                   chemical elements of the sample with their multiplicity
                   element      : string
                                  chemical element
                   multiplicity : int
                                  chemical element multiplicity
    Q            : numpy array
                   momentum transfer (nm^-1)
    
    Returns
    -------
    Iincoh_Q     : numpy array
                   incoherent scattering intensity
    """
    
    # file = open("./incohParamCEA.txt", "r")
    file = open(incoh_path, "r")
    header1 = file.readline()
    lines = file.readlines()
    file.close()
    
    Iincoh_Q = 0
    
    # scan the lines and when it find the right element take the Z number
    for element, multiplicity in elementList.items():
        aff = calc_aff(element, Q, aff_path)
        for line in lines:
            columns = line.split()
            if columns[0] == element:
                Z = float(columns[1])
                M = float(columns[2])
                K = float(columns[3])
                L = float(columns[4])
                break
        
        Iincoh_Q += multiplicity * ((Z - aff**2/Z ) * (1 - M * (np.exp(-K*Q/(4*10*np.pi)) - np.exp(-L*Q/(4*10*np.pi)))))

    return Iincoh_Q


def calc_alpha(J_Q, Sinf, Q, Isample_Q, fe_Q, Ztot, rho0, index):
    """Function to calculate the normalization factor alpha (eq. 34)

    arguments:
    J_Q: J(Q) - array
    Sinf: Sinf - number
    Q: momentum transfer - array
    Isample_Q: sample scattering intensity - array
    fe_Q: effective electric form factor - array
    Ztot: total Z number - number
    rho0: average atomic density - number
    index: array index of element in the calculation range - array

    returns:
    alpha: normalization factor - number
    """

    Integral1 = simps((J_Q[index] + Sinf) * Q[index]**2, Q[index])
    Integral2 = simps((Isample_Q[index]/fe_Q[index]**2) * Q[index]**2,Q[index])
    alpha = Ztot**2 * (((-2*np.pi**2*rho0) + Integral1) / Integral2)

    # DeltaQ = np.diff(Q)
    # meanDeltaQ = np.mean(DeltaQ)
    # Int1 = np.sum((J_Q[index] + Sinf) * Q[index]**2) * meanDeltaQ
    # Int2 = np.sum( (Isample_Q[index]/fe_Q[index]**2) * Q[index]**2  ) * meanDeltaQ
    # alpha = Ztot**2 * (((-2*np.pi**2*rho0) + Int1) / Int2)

    return alpha


def calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, min_index, max_index, calculation_index):
    """Function to calculate the structure factor S(Q) (eq. 18)

    arguments:
    N: number of atoms - number
    Icoh: cohrent scattering intensity - array
    Ztot: total Z number - number
    fe_Q: effective electric form factor - array
    Sinf: Sinf - number
    Q: momentum transfer - array
    min_index: array index of element with Q<minQ - array
    max_index: array index of element with Q>QmaxIntegrate & Q<=maxQ - array
    calculation_index: array index of element in the calculation range Q>minQ & Q<=QmaxIntegrate - array

    returns:
    S_Q: structure factor - array
    """

    S_Q = Icoh_Q[calculation_index] / (N * Ztot**2 * fe_Q[calculation_index]**2)

    S_Qmin = np.zeros(Q[min_index].size)
    S_Q = np.concatenate([S_Qmin, S_Q])

    S_Qmax = np.zeros(Q[max_index].size)
    S_Qmax.fill(Sinf)
    S_Q = np.concatenate([S_Q, S_Qmax])

    return S_Q


def calc_SQ(N, Icoh_Q, Ztot, fe_Q, Sinf, Q, max_index, integration_index):
    """Function to calculate the structure factor S(Q) (eq. 18) with Igor range.
    This function doesn't set the value 0 for Q<minQ!!!
    
    Parameters
    ----------
    N                 : int
                        number of atoms
    Icoh_Q            : numpy array
                        cohrent scattering intensity
    Ztot              : int
                        total Z number
    fe_Q              : numpy array 
                        effective electric form factor
    Sinf              : float
                        Sinf
    Q                 : numpy array 
                        momentum transfer (nm^-1)
    max_index         : numpy array 
                        index of element with Q>QmaxIntegrate & Q<=maxQ
    integration_index : numpy array 
                        index of element in the integration range Q<=QmaxIntegrate
    
    Returns
    -------
    S_Q               : numpy array
                        structure factor
    """
    
    S_Q = Icoh_Q[integration_index] / (N * Ztot**2 * fe_Q[integration_index]**2)
    
    S_Qmax = np.zeros(Q[max_index].size)
    S_Qmax.fill(Sinf)
    S_Q = np.concatenate([S_Q, S_Qmax])
    
    return S_Q


def smoothing(X, f_X, smoothfactor):
    """Function for smoothing
    """

    smooth = interpolate.UnivariateSpline(X, f_X, k=3, s=smoothfactor)
    smoothedf_X = smooth(X)

    return smoothedf_X


def SQsmoothing(Q, S_Q, Sinf, smoothfactor, min_index, max_index, validation_index):
    """Function for smoothing S(Q)
    """

    S_Qsmooth = smoothing(Q[validation_index], S_Q[validation_index], smoothfactor)

    S_Qsmooth[min_index] = 0.0
    S_Qmax = np.zeros(Q[max_index].size)
    S_Qmax.fill(Sinf)
    S_Qsmooth = np.concatenate([S_Qsmooth, S_Qmax])
    # S_Qsmooth[max_index] = Sinf

    return S_Qsmooth


def SQsmoothing(Q, S_Q, Sinf, smoothfactor, min_index, max_index, calculation_index):
    """Function for smoothing S(Q)
    """

    S_Qsmooth = smoothing(Q[calculation_index], S_Q[calculation_index], smoothfactor)

    S_Qmin = np.zeros(Q[min_index].size)
    S_Qsmooth = np.concatenate([S_Qmin, S_Qsmooth])

    S_Qmax = np.zeros(Q[max_index].size)
    S_Qmax.fill(Sinf)
    S_Qsmooth = np.concatenate([S_Qsmooth, S_Qmax])

    return S_Qsmooth


def calc_SQsmoothing(Q, S_Q, Sinf, smooth_factor, min_index, minQ, QmaxIntegrate, maxQ, NumPoints):
    """Function for smoothing S(Q).
    This function smooths S(Q) and resets the number of points for the variable Q.
    
    Parameters
    ----------
    Q             : numpy array 
                    momentum transfer (nm^-1)
    S_Q           : numpy array 
                    structure factor
    Sinf          : float
                    Sinf
    smooth_factor : float
                    smoothing factor
    min_index     : numpy array 
                    indices of elements with Q<=minQ
    minQ          : float
                    minimum Q value
    maxQ          : float
                    maximum Q value
    QmaxIntegrate : float
                    maximum Q value for the intagrations
    NumPoints     : int
                    number of points in the smoothed S(Q)
    
    Returns
    -------
    newQ          : numpy array 
                    new set of Q with NumPoints dimension 
    S_Qsmoothed   : numpy array
                    smoothed S(Q) with NumPoints dimension
    """
    
    mask_smooth = np.where((Q>minQ) & (Q<=maxQ))
    smooth = interpolate.UnivariateSpline(Q[mask_smooth], S_Q[mask_smooth], k=3, s=smooth_factor)
    newQ = np.linspace(np.amin(Q), maxQ, NumPoints, endpoint=True)
    S_Qsmoothed = smooth(newQ)
    
    # mask_low = np.where(Q<=minQ)
    num_low = S_Qsmoothed[newQ<minQ].size
    smooth = interpolate.UnivariateSpline(Q[min_index], S_Q[min_index], k=3, s=smooth_factor)
    newQLow = np.linspace(np.amin(newQ), minQ, num_low, endpoint=True)
    S_QsmoothLow = smooth(newQLow)
    
    S_Qsmoothed[newQ<minQ] = S_QsmoothLow
    S_Qsmoothed[(newQ>QmaxIntegrate) & (newQ<=maxQ)] = Sinf
    
    return (newQ, S_Qsmoothed)


def calc_SQsmoothing(Q, newQ, S_Q, Sinf, smooth_factor, minQ, QmaxIntegrate, maxQ):
    """Function for smoothing S(Q).
    This function smooths S(Q) and resets the number of points for the variable Q.
    
    Parameters
    ----------
    Q             : numpy array 
                    momentum transfer (nm^-1)
    S_Q           : numpy array 
                    structure factor
    Sinf          : float
                    Sinf
    smooth_factor : float
                    smoothing factor
    minQ          : float
                    minimum Q value
    QmaxIntegrate : float
                    maximum Q value for the intagrations
    maxQ          : float
                    maximum Q value
    NumPoints     : int
                    number of points in the smoothed S(Q)
    
    Returns
    -------
    S_Qsmoothed   : numpy array
                    smoothed S(Q) with NumPoints dimension
    """
    
    smooth = interpolate.UnivariateSpline(Q[(Q>minQ) & (Q<=QmaxIntegrate)], \
        S_Q[(Q>minQ) & (Q<=QmaxIntegrate)], k=3, s=smooth_factor)
    # newQ = np.linspace(np.amin(Q), maxQ, NumPoints, endpoint=True)
    S_Qsmoothed = smooth(newQ)
    
    S_Qsmoothed[newQ<minQ] = 0
    S_Qsmoothed[(newQ>QmaxIntegrate) & (newQ<=maxQ)] = Sinf
    
    return S_Qsmoothed


def calc_SQdamp(S_Q, Q, Sinf, QmaxIntegrate, damping_factor):
    """
    """

    # damping_factor = 0.5
    # damping_factor = np.log(10)
    exponent_factor = damping_factor / QmaxIntegrate**2
    damp_Q = np.exp(-exponent_factor * Q**2)

    S_Qdamp = (damp_Q * (S_Q - Sinf)) + Sinf

    return S_Qdamp


def calc_damp(Q, QmaxIntegrate, damping_factor):
    """Function to calculate the damping function
    """

    # damping_factor = 0.5 # np.log(10)
    exponent_factor = damping_factor / QmaxIntegrate**2
    damp_Q = np.exp(-exponent_factor * Q**2)

    return damp_Q


def calc_r(Q):
    """Function to calculate the r value and range used into F(r) calculation

    """

    DeltaQ = np.diff(Q)
    meanDeltaQ = np.mean(DeltaQ)
    r = fftpack.fftfreq(Q.size, meanDeltaQ)
    mask = np.where(r>=0)

    return (r[mask], r, mask)


def calc_Fr(r, Q, i_Q):
    """Function to calculate F(r) (eq. 20)

    arguments:
    r: radius - array
    Q: momentum transfer - array
    i_Q: i(Q) - array

    returns:
    F_r: F(r) - array
    """

    F_r = (2.0 / np.pi) * simps(Q * i_Q * np.array(np.sin(np.mat(Q).T * np.mat(r))).T, Q)

    DeltaQ = np.diff(Q)
    meanDeltaQ = np.mean(DeltaQ)
    rQ = np.outer(r,Q)
    sinrQ = np.sin(rQ)
    F_r2 = (2.0 / np.pi) * np.sum(Q * i_Q * sinrQ, axis=1) * meanDeltaQ

    return (F_r, F_r2)


def calc_Fr(Q, Qi_Q, QmaxIntegrate):
    """Function to calculate F(r) (eq. 20) with the FFT.

    Parameters
    ----------
    Q             : numpy array
                    momentum transfer (nm^-1)
    Qi_Q          : numpy array 
                    Qi(Q)
    QmaxIntegrate : float
                    maximum Q value for the intagration
    
    Returns
    -------
    F_r  : numpy array
           F(r)
    """

    DeltaQ = np.diff(Q)

    # mask = np.where(Q<=QmaxIntegrate)
    # print(len(mask[0]))
    # print(len(Qi_Q))
    # zeroPad = np.zeros(2**275)
    # Qi_Q = np.concatenate([Qi_Q, zeroPad])

    # F_r2 = np.fft.fft(Qi_Q)
    F_r2 = fftpack.fft(Qi_Q, 2**12)
    F_r2 = np.imag(F_r2)
    F_r2 *= DeltaQ[0] *2/np.pi

    F_r3 = fftpack.dst(Qi_Q, type=2, n=2**11)
    F_r3 *= DeltaQ[0] * 2/np.pi

    # for i in range(len(F_r2)):
        # print(F_r2[i], "----", F_r2imag[i], "----", F_r3[i])

    return (F_r2, F_r3)


def calc_NewDimFFT(newQ, maxQ, QmaxIntegrate, Qi_Q):
    """Function to redefine the array dimension to use the FFT.
    """

    idx, elem = find_nearest(newQ, QmaxIntegrate)
    newDim = 2*2*2**math.ceil(math.log(5*(idx+1))/math.log(2))
    Qi_Q2 = np.resize(Qi_Q, newDim)
    Qi_Q2[idx:] = 0.0

    print(len(Qi_Q2))
    print(newDim)

    newQ2 = np.linspace(np.amin(newQ), maxQ, newDim, endpoint=True)
    print(len(newQ2))
    # DeltaQ = np.diff(newQ)
    # deltaR = 2*np.pi/(DeltaQ[0]*newDim)

    return (newQ2, Qi_Q2)


def calc_gr(r, F_r, rho0):
    """Function to calculate g(r)

    """

    g_r = 1 + F_r / (4 * np.pi * r * rho0)

    return g_r


def calc_iQi(i_Q, Q, Sinf, J_Q, deltaF_r, r, rmin):
    """Function to calculate the i-th iteration of i(Q) (eq. 46, 49)

    arguments:
    i_Q: i(Q) - array
    Q: momentum transfer - array
    Sinf: Sinf - number
    J_Q: J(Q) - array
    deltaF_r: deltaF(r) - array
    rmin: value of r cutoff - number

    returns:
    i_Qi: i-th iteration of i(Q) - array
    """

    mask = np.where(r < rmin)
    rInt = r[mask]
    deltaF_rInt = deltaF_r[mask]

    integral = simps(deltaF_rInt * (np.array(np.sin(np.mat(rInt).T *  np.mat(Q)))).T, rInt)

    i_Qi = i_Q - ( 1/Q * ( i_Q / (Sinf + J_Q) + 1)) * integral

    return i_Qi


def calc_iintra(Q, max_index):
    """Function to calculate the intramolecular contribution of i(Q) (eq. 41)

    To implemente!!! -> For now just for CO2!!!
    """

    # Fintra_r = np.zeros(r.size)

    # dCO = 0.1165 # nm
    dCO = 0.1514076 # nm
    dOO = 2 * dCO

    elementList = {"C":1,"O":2}
    fe_Q, Ztot = calc_eeff(elementList, Q)
    KC = calc_Kp(fe_Q, "C", Q)
    KO = calc_Kp(fe_Q, "O", Q)

    constCO = 4/Ztot**2
    constOO = 2/Ztot**2

    sinCO = np.zeros(Q.size)
    sinOO = np.zeros(Q.size)

    for i in range(Q.size):
        if Q[i] == 0.0:
            sinCO[i] = 1
            sinOO[i] = 1
        else:
            sinCO[i] = np.sin(dCO*Q[i])/(dCO*Q[i])
            sinOO[i] = np.sin(dOO*Q[i])/(dOO*Q[i])

    iintra_Q_CO = constCO * KC * KO * sinCO
    iintra_Q_OO = constOO * KO * KO * sinOO

    iintra_Q = iintra_Q_CO + iintra_Q_OO

    iintra_Q[max_index] = 0.0

    return iintra_Q


def calc_iintra(Q, max_index, elementList, element, x, y, z, incoh_path, aff_path):
    """Function to calculate the intramolecular contribution of i(Q) (eq. 41).
    
    Parameters
    ----------
    Q, max_index, elementList, element, x, y, z, incoh_path, aff_path
    
    Q           : numpy array
                  momentum transfer (nm^-1)
    max_index   : numpy array 
                  index of element with Q>QmaxIntegrate & Q<=maxQ
    elementList : dictionary("element": multiplicity)
                  chemical elements of the sample with their multiplicity
                  element      : string
                                 chemical element
                  multiplicity : int
                                 chemical element multiplicity
    element     : string array
                  array with the elements in the xyz_file
    x, y, z     : float
                  atomic coordinate in the xyz_file (nm)
    incoh_path  : string
                  path of incoherent scattered intensities parameters
    aff_path    : string
                  path of atomic scattering form factor parameters
    
    Returns
    -------
    iintra_Q    : numpy array
                  intramolecular contribution of i(Q)
    fe_Q        : numpy array
                  effective electric form factor
    """

    fe_Q, Ztot = calc_eeff2(elementList, Q, incoh_path, aff_path)

    # numAtoms, element, x, y, z = read_xyz_file(path)
    iintra_Q = np.zeros(Q.size)
    sinpq = np.zeros(Q.size)

    for ielem in range(len(element)):
        for jelem in range(len(element)):
            if ielem != jelem:
                KK = calc_Kp(fe_Q, element[ielem], Q, aff_path) * calc_Kp2(fe_Q, element[jelem], Q, aff_path)
                d = calc_distMol(x[ielem], y[ielem], z[ielem], x[jelem], y[jelem], z[jelem])
                if d != 0.0:
                    iintra_Q += KK * np.sin(d*Q) / (d*Q)
                    iintra_Q[Q==0.0] = KK

    iintra_Q[max_index] = 0.0
    iintra_Q /= Ztot**2

    return (iintra_Q, fe_Q)


def calc_Fintra(r, Q, QmaxIntegrate, aff_path):
    """Function to calculate the intramolecular contribution of F(r) (eq. 42).
    To implemente!!! -> For now just for CO2!!!
    For now I calculate Fintra from iintra.
    
    Parameters
    ----------
    r             : numpy array
                    atomic distance (nm)
    Q             : numpy array
                    momentum transfer (nm^-1)
    QmaxIntegrate : float
                    maximum Q value for the intagration
    aff_path      : string
                    path of atomic scattering form factor parameters
    
    
    Returns
    -------
    Fintra_r      : numpy array
                    intramolecular contribution of F(r)
    """

    # Fintra_r = np.zeros(r.size)

    dCO = 0.1165 # nm
    dOO = 2 * dCO

    elementList = {"C":1,"O":2}
    fe_Q, Ztot = calc_eeff(elementList, Q)
    KC = calc_Kp(fe_Q, "C", Q, aff_path)
    KO = calc_Kp(fe_Q, "O", Q, aff_path)

    constCO = 4/(np.pi * Ztot**2 * dCO)
    constOO = 2/(np.pi * Ztot**2 * dOO)

    Fintra_r_CO = constCO * KC * KO * \
        ((np.sin((r - dCO)*QmaxIntegrate)) / (r - dCO) - (np.sin((r + dCO)*QmaxIntegrate)) / (r + dCO))

    Fintra_r_OO = constOO * KO * KO * \
        ((np.sin((r - dOO)*QmaxIntegrate)) / (r - dOO) - (np.sin((r + dOO)*QmaxIntegrate)) / (r + dOO))

    Fintra_r = Fintra_r_CO + Fintra_r_OO

    return Fintra_r


def calc_deltaMinim(N, r, Q, rho0, s, Sinf, I_Q, I_Qbkg, Iincoh, J_Q, fe_Q, Ztot, Fintra):
    """Function to minimize the density

    """
    #numAtoms = sc.N_A
    Isample = I_Q - s * I_Qbkg
    alpha = calc_alpha(J_Q, Sinf, Q, I_Q, fe_Q, Ztot, rho0)
    Icoh = (N * alpha * Isample) - (N * Iincoh)
    S_Q = calc_SQ(Icoh, Ztot, fe_Q)
    i_Q = calc_iQ(S_Q, Sinf)
    F_r = calc_Fr(r, Q, i_Q)
    observed = F_r
    excepted = Fintra - (4*np.pi*rho0)

    chi2, p = chisquare(observed, f_exp=excepted)
    return chi2


def calc_alphaFZ(numAtoms, Q, Isample_Q, Iincoh_Q, rho0, elementParameters):
    """Function to calcultate alpha for the FZ formalism.
    
    """
    
    f2 = MainFunctions.calc_aff('C', Q, elementParameters)**2 + 2*MainFunctions.calc_aff('O',Q, elementParameters)**2
    f2 /= numAtoms
    f = MainFunctions.calc_aff('C', Q, elementParameters)**2 + \
        4*MainFunctions.calc_aff('C', Q, elementParameters)*MainFunctions.calc_aff('O',Q, elementParameters) + \
        4*MainFunctions.calc_aff('O',Q, elementParameters)**2
    f /= numAtoms
    
    Integral1 = simps((Iincoh_Q + (f2/f)) * Q**2, Q)
    Integral2 = simps((Isample_Q/f) * Q**2,Q)
    alpha = ((-2*np.pi**2*rho0) + Integral1) / Integral2
    
    return alpha


def calc_S_QFZ(numAtoms, Icoh_Q, Ztot, Q, elementParameters):
    """Function to calculate S(Q) with Faber-Ziman formalism.
    
    Parameters
    ----------
    Icoh_Q            : numpy array (nm)
                        coherent scattering intensity
    Q                 : numpy array
                        momentum transfer (nm)
    Sinf              : float
                        Sinf
    min_index         : numpy array
                        array index of element with Q<=minQ
    max_index         : numpy array
                        array index of element with Q>QmaxIntegrate & Q<=maxQ
    calculation_index : numpy array
                        range where S(Q) is calculated
    
    
    Returns
    -------
    S_Q               : numpy array
                        structure factor in FZ formalism
    """
    
    f2 = MainFunctions.calc_aff('C', Q, elementParameters)**2 + 2*MainFunctions.calc_aff('O',Q, elementParameters)**2
    f2 /= numAtoms
    f = MainFunctions.calc_aff('C', Q, elementParameters)**2 + \
        4*MainFunctions.calc_aff('C', Q, elementParameters)*MainFunctions.calc_aff('O',Q, elementParameters) + \
        4*MainFunctions.calc_aff('O',Q, elementParameters)**2
    f /= numAtoms
   
    S_Q = (Icoh_Q - (f2 - f)) / (f)
    # S_Q /= (numAtoms*Ztot**2)

    # S_Qmin = np.zeros(Q[min_index].size)
    # S_Q = np.concatenate([S_Qmin, S_Q])

    # S_Qmax = np.ones(Q[max_index].size)
    # # S_Qmax.fill(Sinf)
    # S_Q = np.concatenate([S_Q, S_Qmax])

    return (S_Q, f2, f)


def diamond(path, type, Q, I_Q, diamond_angle):
    """Function to calculate the diamond correction.
    http://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z06.html
    """
    
    file = open(path, "r")
    header1 = file.readline()
    header2 = file.readline()
    lines = file.readlines()
    file.close()
    
    for line in lines:
        columns = line.split()
        if columns[0] == type:
            dimension = float(columns[1])
            break
            
    # for now...
    # mu = 0.2562 # cm2/g
    wavelenght = 0.03738 # nm
    # diamond_density = 3.51 # g/cm3
    # mu_l = mu * diamond_density # cm
    mu_l = 1/12.08 # mm^-1 at 33 keV
    # theta_angle = np.arcsin(wavelenght*Q/(4*np.pi))
    _2theta_angle = Qto2theta(Q) # rad
    
    
    
    path_lenght = dimension / np.cos(_2theta_angle)
    
    corr_factor = np.exp(mu_l * path_lenght)
    
    I_Qeff = corr_factor * I_Q
    
    return (I_Qeff, corr_factor)
    
    
def calc_T_MCC_samp(phi, sth):
    """Function to calculate the MCC transmission for the sample.
    """
    
    DeltaSth = np.diff(sth)
    meanDeltaSth = np.mean(DeltaSth)
    T_MCC_samp = np.sum(phi, axis=0) * meanDeltaSth
    
    return T_MCC_samp


def calc_T_MCC_sample(phi_matrix):
    """Function to calculate the MCC sample transfer function.
    
    Parameters
    ----------
    phi_matrix : 2D numpy array
                 dispersion angle matrix (rad)
    
    Returns
    -------
    T_MCC_sample : numpy array
                 MCC sample transfer function
    """
    
    T_MCC_sample = simps(phi_matrix, axis=0, even="first")
        
    return T_MCC_sample


def calc_T_MCC_DAC(phi_matrix, T_MCC_sample):
    """Function to calculate the MCC DAC transfer function.
    
    Parameters
    ----------
    phi_matrix   : 2D numpy array
                   dispersion angle matrix (rad)
    T_MCC_sample : numpy array
                   MCC sample transfer function
    
    Returns
    -------
    T_MCC_DAC    : numpy array
                   MCC DAC transfer function
    """
    
    T_MCC_ALL = simps(phi_matrix, axis=0, even="first") 
    T_MCC_DAC = T_MCC_ALL - T_MCC_sample
        
    return (T_MCC_ALL, T_MCC_DAC)


def calc_T_DAC_MCC_bkg_corr(I_Qbkg, T_DAC_MCC_sth, T_DAC_MCC_s0th):
    """Function to calculate the bkg correction for the MCC (W. eq. 12)
    
    Parameters
    ----------
    I_Qbkg         : numpy array
                     measured scattering intensity for the bkg
    T_DAC_MCC_sth  : numpy array
                     MCC DAC transfer function
    T_DAC_MCC_s0th : numpy array
                     MCC DAC transfer function with sample thickness for the reference spectra
    
    Returns
    -------
    I_Qbkg_corr    : numpy array
                     corrected scattering intensity for the bkg
    """
    
    
    corr_factor = np.zeros(T_DAC_MCC_sth.size)
    
    for i in range(len(T_DAC_MCC_sth)):
        if (T_DAC_MCC_sth[i] == 0.0 or T_DAC_MCC_s0th[i] == 0.0):
            corr_factor[i] = 1
        else:
            corr_factor[i] = T_DAC_MCC_sth[i] / T_DAC_MCC_s0th[i]
    
    I_Qbkg_corr = corr_factor * I_Qbkg
    
    return (corr_factor, I_Qbkg_corr)


def calc_T_DAC_MCC_bkg_corr(I_Qbkg, T_DAC_MCC_sth, T_DAC_MCC_s0th):
    """Function to calculate the bkg correction for the MCC (W. eq. 12)
    
    Parameters
    ----------
    I_Qbkg         : numpy array
                     measured scattering intensity for the bkg
    T_DAC_MCC_sth  : numpy array
                     MCC DAC transfer function
    T_DAC_MCC_s0th : numpy array
                     MCC DAC transfer function with sample thickness for the reference spectra
    
    Returns
    -------
    corr_factor    : numpy array
                     the correction factor
    I_Qbkg_corr    : numpy array
                     corrected scattering intensity for the bkg
    """
    
    np.seterr(divide="ignore")
    corr_factor = T_DAC_MCC_sth / T_DAC_MCC_s0th
    corr_factor[T_DAC_MCC_sth == 0.0] = 1
    corr_factor[T_DAC_MCC_s0th == 0.0] = 1
    
    I_Qbkg_corr = corr_factor * I_Qbkg
    
    return (corr_factor, I_Qbkg_corr)


def calc_indices(Q, minQ, QmaxIntegrate, maxQ):
    """Function to calculate the Q ranges where S(Q) is constant.
    
    Parameters
    ----------
    Q             : numpy array 
                    momentum transfer (nm^-1)
    minQ          : float
                    minimum Q value
    QmaxIntegrate : float
                    maximum Q value for the intagrations
    maxQ          : float
                    maximum Q value
    
    Returns
    -------
    min_index     : numpy array
                    indices of elements with Q<=minQ
    max_index     : numpy array
                    indices of elements with Q>QmaxIntegrate & Q<=maxQ
    """

    min_index = np.where(Q<=minQ)
    max_index = np.where((Q>QmaxIntegrate) & (Q<=maxQ))

    return (min_index, max_index)


def calc_ranges(Q, minQ, QmaxIntegrate, maxQ):
    """Function to calculate the Q ranges used in the program.
    
    Parameters
    ----------
    Q                 : numpy array 
                        momentum transfer (nm^-1)
    minQ              : float
                        minimum Q value
    QmaxIntegrate     : float
                        maximum Q value for the intagrations
    maxQ              : float
                        maximum Q value
    
    Returns
    -------
    validation_index  : numpy array
                        range of valide Q (Q<=maxQ)
    integration_index : numpy array
                        range where the integration is calculated (Q<=QmaxIntegrate)
    calculation_index : numpy array
                        range where S(Q) is calculated (minQ<Q<=QmaxIntegrate)
    """

    validation_index = np.where(Q<=maxQ)
    integration_index = np.where(Q<=QmaxIntegrate)
    calculation_index = np.where((Q>minQ) & (Q<=QmaxIntegrate))

    return (validation_index, integration_index, calculation_index)


def interpolation(X, f_X, rebinnedX):
    """Function for the interpolation
    """

    interpolatedf_X = interpolate.interp1d(X, f_X)
    newf_X = interpolatedf_X(rebinnedX)

    return newf_X


def fitcurve(X, f_X, mask):
    """Function to flat the peak
    """

    xpoints = X[mask]
    ypoints = f_X[mask]

    coefficients = np.polyfit(xpoints, ypoints, 2)
    polynomial = np.poly1d(coefficients)
    y_axis = polynomial(xpoints)

    return y_axis


def calc_NewDimFFT():
    """Function to redefine the array dimension to use the FFT
    """

    idx, elem = find_nearest(newQ, QmaxIntegrate)
    newDim = 2*2*2**math.ceil(math.log(5*(idx+1))/math.log(2))
    Qi_Q2 = np.resize(Qi_Q, newDim)
    Qi_Q2[idx:] = 0.0

    newQ2 = np.linspace(np.amin(Q), maxQ, newDim, endpoint=True)

    DeltaQ = np.diff(newQ)
    deltaR = 2*np.pi/(DeltaQ[0]*newDim)


def rebinning(X, f_X, BinNum, Num, maxQ, minQ):
    """Function for the rebinning
    """

    newf_X = interpolate.interp1d(X, f_X)
    ShitX = np.linspace(np.amin(X), maxQ, BinNum*Num, endpoint=True)
    ShitY = newf_X(ShitX)

    min = (BinNum - 1)/2 * maxQ /(BinNum * Num - 1)
    max = maxQ - (BinNum - 1)/2 * maxQ / (BinNum*Num - 1)
    BinX = np.linspace(min, max, Num, endpoint=True)
    BinY = np.zeros(Num)

    for i in range(BinNum):
        for j in range(0, Num):
            BinY[j] += ShitY[j*BinNum+i]

    BinY /= BinNum

    mask = np.where(X<=minQ)
    BinY[mask] = 0.0

    # lenX = len(X)
    # numX = 2**int(math.log(lenX,2))
    # rebinnedX = np.linspace(np.amin(X), maxQ, numX, endpoint=True)
    # if min < np.amin(X):
        # min = np.amin(X)

    return (BinX, BinY)


def plot_raw_data(xVal, yVal, plotName, xName, yName, labName):
    """Function to plot the raw data.

    Parameters
    ----------
    xVal     : numpy array
               abscissa values
    yVal     : numpy array
               ordinate values
    plotName : string
               canvas name
    xName    : string
               abscissa name
    yName    : string
               ordinate name
    labName  : string
               label name
    """

    plt.figure(plotName)
    plt.plot(xVal, yVal, label=labName)
    plt.xlabel(xName)
    plt.ylabel(yName)
    plt.legend()
    plt.grid(True)
    plt.draw()