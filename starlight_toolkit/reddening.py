'''
Created on 08/24/2016

@author: ariel

Provides functions to correct spectra for galactic extinction
'''

import numpy as np
from astropy.coordinates import SkyCoord as sc


def CCM(lamb, Rv= 3.1):
    '''

    Calculates the Cardelli, Clayton & Mathis (CCM) extinction curve in the
    optical, wavelengths should be in the 3030-9090 Angstrons range.

    Input:   Wavelengths, Rv (Optional, default is 3.1)
    Returns: A_lambda/Av

    Reference: http://adsabs.harvard.edu/abs/1989ApJ...345..245C

    '''

    #Turn lambda from angstrons to microns:
    lamb = lamb / 10000.

    #Calculate A_lambda/Av:
    x =  1./lamb
    y =  (x-1.82)

    a =  1.+(0.17699 * y) - ( 0.50447 * (y **2) ) - ( 0.02427 * (y **3) )
    a += ( 0.72085 * (y **4) ) + (0.01979 * (y **5) ) - ( 0.77530 * (y **6) )
    a += (0.32999 * (y **7) )

    b =  (1.41338 * y) + ( 2.28305 * (y **2) ) + ( 1.07233 * (y **3) )
    b += -( 5.38434 * (y **4) ) - (0.62251 * (y **5) ) + ( 5.30260 * (y **6) )
    b += -(2.09002 * (y **7) )

    A_lambda_Av = a + (b/Rv)

    return A_lambda_Av


def calc_extinction(lamb
    , EBV, Rv=3.1):


    #Calculate Av and A_lambda:
    Av = Rv * EBV
    A_lambda = Av * CCM(lamb)

    return A_lambda
    
    
def get_EBV(l,b,EBVmap):
    #Get the corresponting HEALPix index and the E(B-V) value:
    index = hp.ang2pix(nside=2048,theta=(np.pi/2) - b,phi=l)
    EBV = EBVmap[index]
    
    return EBV



def extinction_corr(spectra,lambdas,EBV):
    '''

    Corrects spectra for the effects of galactic extinction.

    Input: Wavelenghts, Fluxes, RA, Dec
    Returns: Fluxes corrected for the effects of Milky Way dust.

    '''

    #Get the extinction in each wavelenght:
    A_lambdas = calc_extinction(lambdas,EBV)

    #Calculate corrected spectra:
    corr_spectra = spectra * np.exp(A_lambdas/np.log(10))

    return corr_spectra


def extinction_decorr(spectra,lambdas,EBV):
    '''

    Input: Wavelenghts, Fluxes, RA, Dec
    Returns: Reddened fluxes.

    '''

    #Get the extinction in each wavelenght:
    A_lambdas = calc_extinction(lambdas,EBV)

    #Calculate corrected spectra:
    corr_spectra = spectra * np.exp(-A_lambdas/np.log(10))

    return corr_spectra
