import numpy as np

def apply_redshift(wl, z, dest='rest'):
    """
    Apply redshift to wavelengths.

    Parameters:
        wl (ndarray): Array of input wavelengths.
        z (float or ndarray): Redshift value(s).
        dest (str, optional): Destination frame, either 'rest' (default) or 'observed'.

    Returns:
        ndarray: Redshifted wavelengths.
    """
   
    if dest == 'rest':
        op = lambda x, y: x / y
    elif dest == 'observed':
        op = lambda x, y: x * y

    if np.isscalar(z):
        return op(wl, 1. + z)
    else:
        return op(wl[:, np.newaxis], 1. + z[np.newaxis, :])


def spectra2restframe(l_obs, f_obs, z):
    """
    Shifts observed spectra to rest frame.

    Parameters:
        l_obs (ndarray): Array of observed wavelengths.
        f_obs (ndarray): Array of observed fluxes.
        z (float): Redshift value.

    Returns:
        ndarray: Rest-frame wavelengths.
        ndarray: Rest-frame fluxes.
    """
 
    l_rest = apply_redshift(l_obs, z)
    f_rest = f_obs * (1.0 + z)
    return l_rest, f_rest
