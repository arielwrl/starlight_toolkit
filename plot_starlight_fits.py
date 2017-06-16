'''

Tools for plotting Starlight outputs.

ariel@ufsc
Created on 05/30/2017

'''

import numpy as np
import matplotlib.pyplot as plt
from starlight_tables import read_output_table

def plot_starlight_output_spec(out_file, plot_obs=True, plot_error=False, plot_labels=False, obs_color='k', syn_color='b', w0_color='y', syn_label='Fitted Spectrum'):
    '''
    TODO: Add marks for clipped wavelengths.
    
    Note that there is no plt.show() in this function, the intention is that you 
    can use this funcion multiple times to overplot different Starlight outputs.
    '''
    
    out = read_output_table(out_file)
    
    l_obs, f_obs, f_syn, f_wei = out['spectra']['l_obs'], out['spectra']['f_obs'], \
    out['spectra']['f_syn'], out['spectra']['f_wei']
    
    w0 = out['spectra']['f_wei'] <= 0
    
    error = 1/(f_wei[~w0]**2)

    if plot_obs==True:
        
        f_obs_masked = np.ma.masked_array(data = f_obs, mask = w0)
        f_w0         = np.ma.masked_array(data = f_obs, mask = ~w0)
        
        plt.plot(l_obs, f_obs, color=obs_color, lw=0.5, label='Observed Spectrum')
        plt.plot(l_obs, f_w0, color=w0_color, lw=0.5, label=r'$w_\lambda=0$')

    if plot_error==True:
        plt.plot(l_obs[~w0], error, '--r')

    if plot_labels==True:
        plt.ylabel(r'$F_\lambda/F_{\lambda0}$', fontsize=15)
        plt.xlabel(r'$\lambda\mathrm{[\AA]}$', fontsize=15)

    plt.plot(l_obs, f_syn, color=syn_color, lw=0.5, label=syn_label)
    
    plt.ylim(0, 1.4 * np.max(f_obs))
    plt.xlim(3200,12000)


