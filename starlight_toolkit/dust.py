import numpy as np


def calc_extinction(wl, EBV, Rv=3.1):
    """
    Calculates Galactic extinction at a given wavelength.

    Parameters:
        wl (float or array-like): Wavelength(s) in Angstroms.
        EBV (float): Color excess or the ratio of total to selective extinction (E(B-V)).
        Rv (float, optional): Ratio of total to selective extinction. Default is 3.1.

    Returns:
        float or array-like: Extinction at the specified wavelength(s).

    Notes:
        This function will not work if the input is a list.

    Reference:
        Cardelli, J. A., Clayton, G. C., & Mathis, J. S. (1989). 
        The relationship between infrared, optical, and ultraviolet extinction. 
        The Astrophysical Journal, 345, 245–256.
    """
 
    AV = Rv * EBV
    A_lambda = AV * CCM(wl)

    return A_lambda


def extinction_corr(flux, wl, EBV):
    """
    Corrects spectra for the effects of Galactic extinction.

    Parameters:
        flux (array-like): The flux values of the spectra.
        wl (array-like): The wavelengths corresponding to the flux values.
        EBV (float): Ratio of total to selective extinction (E(B-V)).

    Returns:
        array-like: Fluxes corrected for the effects of Milky Way dust.

    Notes:
        This function corrects input spectra for the effects of Galactic extinction due to interstellar dust.
        
    References:
        Cardelli, J. A., Clayton, G. C., & Mathis, J. S. (1989). 
        The relationship between infrared, optical, and ultraviolet extinction. 
        The Astrophysical Journal, 345, 245–256.
    """
 
    A_lambdas = calc_extinction(wl, EBV)

    corr_spectra = flux * 10 ** (0.4 * A_lambdas)

    return corr_spectra


def CCM(wl, R_V=3.1):
    """
    Calculates the Cardelli, Clayton & Mathis (CCM) extinction curve.

    Parameters:
        wl (array-like): Wavelength in Angstroms.
        R_V (float, optional): Ratio of total to selective extinction. Default is 3.1.

    Returns:
        array-like: Extinction A_lambda / A_V. Array of the same length as wl.

    References:
        Cardelli, J. A., Clayton, G. C., & Mathis, J. S. (1989).
        The relationship between infrared, optical, and ultraviolet extinction.
        The Astrophysical Journal, 345, 245–256.
    """

    a = np.zeros(np.shape(wl))
    b = np.zeros(np.shape(wl))
    F_a = np.zeros(np.shape(wl))
    F_b = np.zeros(np.shape(wl))
    x = np.zeros(np.shape(wl))
    y = np.zeros(np.shape(wl))
    q = np.zeros(np.shape(wl))

    x = 10000. / wl
    y = 10000. / wl - 1.82

    # Far-Ultraviolet: 8 <= x <= 10 ; 1000 -> 1250 Angs
    i = np.bitwise_and(x >= 8, x <= 10)

    a[i] = -1.073 - 0.628 * (x[i] - 8.) + 0.137 * (x[i] - 8.)**2 - 0.070 * (x[i] - 8.)**3
    b[i] = 13.670 + 4.257 * (x[i] - 8.) - 0.420 * (x[i] - 8.)**2 + 0.374 * (x[i] - 8.)**3

    # Ultraviolet: 3.3 <= x <= 8 ; 1250 -> 3030 Angs
    i =  np.bitwise_and(x >= 5.9, x < 8)
    F_a[i] = -0.04473 * (x[i] - 5.9)**2 - 0.009779 * (x[i] - 5.9)**3
    F_b[i] =  0.2130 * (x[i] - 5.9)**2 + 0.1207 * (x[i] - 5.9)**3

    i =  np.bitwise_and(x >= 3.3, x < 8)

    a[i] =  1.752 - 0.316 * x[i] - 0.104 / ((x[i] - 4.67)**2 + 0.341) + F_a[i]
    b[i] = -3.090 + 1.825 * x[i] + 1.206 / ((x[i] - 4.62)**2 + 0.263) + F_b[i]

    # Optical/NIR: 1.1 <= x <= 3.3 ; 3030 -> 9091 Angs ;
    i = np.bitwise_and(x >= 1.1, x < 3.3)

    a[i] = 1.+ 0.17699 * y[i] - 0.50447 * y[i]**2 - 0.02427 * y[i]**3 + \
        0.72085 * y[i]**4 + 0.01979 * y[i]**5 - 0.77530 * y[i]**6 + 0.32999 * y[i]**7
    b[i] = 1.41338 * y[i] + 2.28305 * y[i]**2 + 1.07233 * y[i]**3 - \
        5.38434 * y[i]**4 - 0.62251 * y[i]**5 + 5.30260 * y[i]**6 - 2.09002 * y[i]**7

    # Infrared: 0.3 <= x <= 1.1 ; 9091 -> 33333 Angs ;
    i = np.bitwise_and(x >= 0.3, x < 1.1)

    a[i] =  0.574 * x[i]**1.61
    b[i] = -0.527 * x[i]**1.61

    q = a + b / R_V

    return q


def CAL(wl, R_V=4.05):
    """
    Calculates the extinction curve according to the Calzetti law.

    Parameters:
        wl (array-like): Wavelength in Angstroms.
        R_V (float, optional): Ratio of total to selective extinction. Default is 4.05.

    Returns:
        array-like: Extinction A_lambda / A_V. Array of the same length as wl.

    Notes:
        This function computes the extinction curve described by the Calzetti law.
        
        Wavelength ranges:
        - UV to Optical: 1200 -> 6300 Angstroms
        - Red to Infrared: 6300 -> 25000 Angstroms

        This function issues a warning if the wavelength falls outside the valid range (1200 to 25000 Angstroms).
 
    References:
        Calzetti, D., Armus, L., Bohlin, R. C., et al. (2000). 
        The Dust Content and Opacity of Actively Star-forming Galaxies. 
        The Astrophysical Journal, 533(2), 682–695.
 
    """
 
    if not isinstance(wl, np.ma.MaskedArray):
        wl = np.asarray(wl, 'float64')

    x = 1e4 / wl
    q = np.zeros_like(wl)

    # UV -> Optical: 1200 -> 6300 Angs
    i = (wl >= 1200.) & (wl <= 6300.)
    q[i] = (2.659 / R_V) * (-2.156 + 1.509 * x[i] - 0.198 * x[i]**2 + 0.011 * x[i]**3) + 1.

    # Red -> Infrared
    i = (wl >= 6300.) & (wl <= 25000.)
    q[i] = (2.659 / R_V) * (-1.857 + 1.040 * x[i]) + 1.

    # Issue a warning if lambda falls outside 1200->22000 Angs range
    if ((wl < 1200.) | (wl > 25000.)).any():
        print('WARNING! some lambda outside valid range (1200, 22000.)')
    return q



def CSB(wl, R_V=3.1, B=0.5):
    """
    Calculates the extinction curve according to the Conroy, Schiminovich, and Blanton (CSB) law.

    Parameters:
        wl (array-like): Wavelength in Angstroms.
        R_V (float, optional): Ratio of total to selective extinction. Default is 3.1.
        B (float, optional): Parameter controlling the strength of the 2175 Å bump. Default is 0.5.

    Returns:
        array-like: Extinction A_lambda / A_V. Array of the same length as wl.

    Notes: 
        This function adds a bump to the Calzetti curve.

    References:
        Conroy, Charlie ; Schiminovich, David ; Blanton, Michael R. (2010).
        Dust Attenuation in Disk-dominated Galaxies: Evidence for the 2175 Å Dust Feature
        The Astrophysical Journal, 718(1), 184–202.
    """
    a   = np.zeros(np.shape(wl))
    b   = np.zeros(np.shape(wl))
    F_a = np.zeros(np.shape(wl))
    F_b = np.zeros(np.shape(wl))
    x   = np.zeros(np.shape(wl))
    y   = np.zeros(np.shape(wl))
    q   = np.zeros(np.shape(wl))

    x = 10000. / wl
    y = x - 1.82

    # Far-Ultraviolet: 8 <= x <= 10 ; 1000 -> 1250 Angs
    i = np.bitwise_and(x >= 8, x <= 11)

    a[i] = -1.073 - 0.628 * (x[i] - 8.) + 0.137 * (x[i] - 8.)**2 - 0.070 * (x[i] - 8.)**3
    b[i] = 13.670 + 4.257 * (x[i] - 8.) - 0.420 * (x[i] - 8.)**2 + 0.374 * (x[i] - 8.)**3

    # Ultraviolet: 3.3 <= x <= 8
    i = np.bitwise_and(x >= 5.9, x < 8)

    F_a[i] = -0.0447*(x[i]-5.9)**2 - 0.00978*(x[i]-5.9)**3
    F_b[i] = 0.213*(x[i]-5.9)**2  + 0.121*(x[i]-5.9)**3

    a[i] = F_a[i] + 1.752 - 0.316*x[i] -  0.104*B / ( (x[i]-4.67)**2 +0.341)
    b[i] = F_b[i] - 3.09 + 1.825*x[i] + 1.206*B / ( (x[i]-4.62)**2 +0.263)

    i = np.bitwise_and(x >= 3.3, x <= 5.9)

    F_a[i] = (3.3/x[i])**6 * (-0.0370 + 0.0469*B - 0.601*(B/R_V) + 0.542/R_V)

    a[i] = F_a[i] + 1.752 - 0.316*x[i] -  0.104*B / ( (x[i]-4.67)**2 +0.341)
    b[i] = -3.09 + 1.825*x[i] + 1.206*B / ( (x[i]-4.62)**2 +0.263)

    # Optical/NIR: 1.1 <= x <= 3.3 ; 3030 -> 9091 Angs ;
    i = np.bitwise_and(x >= 1.1, x < 3.3)

    a[i] = 1.+ 0.177 * y[i] - 0.504 * y[i]**2 - 0.0243 * y[i]**3 + \
        0.721 * y[i]**4 + 0.0198 * y[i]**5 - 0.775 * y[i]**6 + 0.330 * y[i]**7
    b[i] = 1.413 * y[i] + 2.283 * y[i]**2 + 1.072 * y[i]**3 - \
        5.384 * y[i]**4 - 0.622 * y[i]**5 + 5.303 * y[i]**6 - 2.090 * y[i]**7

    # Infrared: 0.3 <= x <= 1.1 ; 9091 -> 33333 Angs ;
    i = np.bitwise_and(x >= 0.3, x < 1.1)

    a[i] =  0.574 * x[i]**1.61
    b[i] = -0.527 * x[i]**1.61

    q = a + b / R_V

    return q


def CCC(wl, R_V=4.05):
    """
    Modified Calzetti curve for far UV wavelengths.

    Parameters:
        wl (array-like): Wavelength in Angstroms.
        R_V (float, optional): Ratio of total to selective extinction. Default is 4.05.

    Returns:
        array-like: Extinction A_lambda / A_V. Array of the same length as wl.

    """

    q = np.zeros(np.shape(wl))
    x = 10000. / wl

    i_CCC = np.bitwise_and(wl>=912, wl<1846)
    q[i_CCC] = (5.472 + 0.671 * x[i_CCC] - 9.218e-3 * x[i_CCC]**2 + 2.620e-3 * x[i_CCC]**3)/R_V

    i_CAL = np.where(wl>=1846)
    q[i_CAL] = CAL(wl[i_CAL], R_V)

    return q


def Salim2018(wl, delta=0, B=0, gamma=350., wl_bump=2175.):
    """
    Generates the modified Calzetti law from Salim et al. (2018).

    This modified law is based on parameterizations from Noll et al. (2009) and Kriek & Conroy (2011).

    Parameters:
        wl (array-like): Wavelength in Angstroms.
        delta (float, optional): Power-law index parameter. Default is 0.
        B (float, optional): Amplitude of the Drude profile. Default is 0.
        gamma (float, optional): Width of the Drude profile. Default is 350.
        wl_bump (float, optional): Central wavelength of the 2175 Å bump. Default is 2175.

    Returns:
        array-like: Extinction A_lambda / A_V. Array of the same length as wl.

    References:
        Salim, S., et al. (2018).
        Dust attenuation curves in the local universe: demographics and new laws for star-forming galaxies and high-redshift analyses.
        The Astrophysical Journal, 859(1), 11.
    """

    q = np.zeros(np.shape(wl))

    #Generate Drude profile:
    D = B * (wl**2) * (gamma**2) / ( (wl**2-wl_bump**2)**2 + (wl**2 * gamma**2))

    #Calculate q_lambda
    q = CCC(wl) * (wl/5500)**delta + D

    return q
