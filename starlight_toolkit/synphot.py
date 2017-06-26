"""

ariel@ufsc
May 2017

Provides tools for synthetic photometry, may be useful on 
Starlight fits including photometry

"""


import numpy as np
from scipy.interpolate import interp1d


def resampler(x_old, y_old, x_new):
    interp = interp1d(x_old, y_old, bounds_error = False
    , fill_value = (0.,0.))

    y_new = interp(x_new)
    
    return y_new
    

def synflux(wl, spectrum, filter_file):
    
    #Reading filter:
    wl_filter, T_filter = np.genfromtxt(filter_file).transpose()
    
    #Resampling filter and spectrum to 1\AA intervals:
    wl_filter_new = np.arange(np.round(wl_filter[0])-5, np.round(wl_filter[-1])+5)
    
    T_filter_new = resampler(wl_filter, T_filter, wl_filter_new)
    spectrum_new = resampler(wl, spectrum, wl_filter_new)
    
    #Convolution time:
    synflux = np.trapz(spectrum_new * T_filter_new * wl_filter_new, dx=1)
    synflux /= np.trapz(wl_filter_new * T_filter_new, dx=1)

    return synflux


def pivot_wavelength(filter_file):

    #Reading and resampling filter:
    wl_filter, T_filter = np.genfromtxt(filter_file).transpose()    
    
    wl_filter_new = np.arange(np.round(wl_filter[0])-5, np.round(wl_filter[-1])+5, 1)
    
    T_filter_new = resampler(wl_filter, T_filter, wl_filter_new)
    
    #Calculating pivot_wavelength
    
    pivot_wl = np.trapz(T_filter_new * wl_filter_new, dx=1) / np.trapz(T_filter_new * (wl_filter_new**-1), dx=1)
    pivot_wl = np.sqrt(pivot_wl)
    
    return pivot_wl


def effective_wavelength(wl, spectrum, filter_file):
    '''
    This is defined as the mean wavelength of the filter weighted by transmission of the filter
    and spectrum of the source
    '''
    
    #Reading filter:
    wl_filter, T_filter = np.genfromtxt(filter_file).transpose()
    
    #Resampling filter and spectrum to 1\AA intervals:
    wl_filter_new = np.arange(np.round(wl_filter[0])-5, np.round(wl_filter[-1])+5,1)
    
    T_filter_new = resampler(wl_filter, T_filter, wl_filter_new)
    spectrum_new = resampler(wl, spectrum, wl_filter_new)
    
    #Convolution time:
    effective_wl = np.trapz(wl_filter_new**2 * spectrum_new * T_filter_new, dx=1)
    effective_wl /= np.trapz(spectrum_new * T_filter_new * wl_filter_new, dx=1)
    
    return effective_wl
    
    
    
    
    
    
