import numpy as np
from starlight_toolkit.synphot import resampler


def convert_x_lambda(popx, wl, wl_0, base_wl, base_f):
    """
    Convert the stellar population vector a different normalization wavelength.

    Parameters:
        popx (ndarray): Stellar population vector.
        wl (float): Target wavelength.
        wl_0 (float): Reference wavelength.
        base_wl (ndarray): Base wavelength array.
        base_f (ndarray): Base flux array.

    Returns:
        ndarray: Population vector normalized at target wavelength.
    """

    base_wl_res = np.arange(np.round(base_wl[0]), np.round(base_wl[-1]), 1)
    base_f_res  = np.array([resampler(base_wl, base_f[i], base_wl_res) for i in range(len(base_f))])
    
    base_f_norm = [base_f_res[i][base_wl_res==wl]/base_f_res[i][base_wl_res==wl_0]  for i in range(len(base_f))]
    base_f_norm = np.array([base_f_norm[i][0]  for i in range(len(base_f))])

    popx_wl = np.array([popx[i]*base_f_norm/np.sum(popx[i]*base_f_norm) for i in range(len(popx))])
    
    return popx_wl


def calc_sfh(age_base, popmu):
    """
    Calculates the star formation history (SFH) and cumulative SFH.

    Parameters:
        age_base (ndarray): Array of stellar population ages.
        popmu (ndarray): Population mass fractions.

    Returns:
       tuple: A tuple containing three elements: agevec (ndarray) - unique ages, sfh (ndarray) - SFH,
        and csfh (ndarray) - cumulative SFH.
    """
 
    agevec = np.unique(age_base)
 
    sfh = [np.sum(popmu[age_base==agevec[i]]) for i in range(len(agevec))]
    sfh /= popmu.sum() 
    sfh *= 100
 
    csfh = np.cumsum(sfh[::-1])
 
    return agevec, sfh, csfh[::-1]


def calc_sfh_x(age_base, popx):
    """
    Calculates the star formation history (SFH) and cumulative SFH for stellar population fractions.

    Parameters:
        age_base (ndarray): Array of stellar population ages.
        popx (ndarray): Stellar population light fractions.

    Returns:
        tuple: A tuple containing three elements: agevec (ndarray) - unique ages, sfh (ndarray) - SFH,
        and csfh (ndarray) - cumulative SFH.
    """
 
    agevec = np.unique(age_base)

    sfh = [np.sum(popx[age_base==agevec[i]]) for i in range(len(agevec))]
    sfh /= popx.sum() 
    sfh *= 100

    csfh = np.cumsum(sfh[::-1])

    return agevec, sfh, csfh[::-1]


def calc_atflux(age_base, popx, age_base_upp=None):
    """
    Calculates the luminosity-weighted age based on population fractions.

    Parameters:
        age_base (ndarray): Array of stellar population ages.
        popx (ndarray): Stellar population light fractions.
        age_base_upp (ndarray, optional): Upper age limits for each age bin. Defaults to None.

    Notes:
        If `age_base_upp` is given, the average is calculated based on the central age of the bins.

    Returns:
        float: Luminosity-weighted age of the object.
    """

    if age_base_upp is not None:
        log_t1 = np.log10(age_base)
        log_t2 = np.log10(age_base_upp)
        log_t = (log_t1 + log_t2) / 2.0
    else:
        log_t = np.log10(age_base)
    return np.sum(log_t * popx) / popx.sum()


def calc_atmass(age_base, popmu, age_base_upp=None):
    """
    Calculates the mass-weighted age based on population fractions.

    Parameters:
        age_base (ndarray): Array of stellar population ages.
        popx (ndarray): Stellar population mass fractions.
        age_base_upp (ndarray, optional): Upper age limits for each age bin. Defaults to None.

    Notes:
        If `age_base_upp` is given, the average is calculated based on the central age of the bins.

    Returns:
        float: Mass-weighted age of the object.
    """

    if age_base_upp is not None:
        log_t1 = np.log10(age_base)
        log_t2 = np.log10(age_base_upp)
        log_t  = (log_t1 + log_t2) / 2.0
    else:
        log_t  = np.log10(age_base)        
    return np.sum(log_t * popmu) / popmu.sum()

   
def calc_aZflux(Z_base, popx, Z_sun): 
    """
    Calculates the luminosity-weighted metallicity based on population fractions.

    Parameters:
        Z_base (ndarray): Array of metallicities for the stellar populations.
        popx (ndarray): Stellar population fractions.
        Z_sun (float): Solar metallicity value.

    Returns:
        float: Luminosity-weighted metallicity.
    """

    return (popx * np.log10(Z_base/Z_sun)).sum()/popx.sum()


def calc_aZmass(Z_base, popmu, Z_sun): 
    """
    Calculates the mass-weighted metallicity based on population fractions.

    Parameters:
        Z_base (ndarray): Array of metallicities for the stellar populations.
        popx (ndarray): Stellar population fractions.
        Z_sun (float): Solar metallicity value.

    Returns:
        float: Mass-weighted metallicity.
    """

    return (popmu * np.log10(Z_base/Z_sun)).sum()/popmu.sum()




