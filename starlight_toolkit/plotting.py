'''

Tools for plotting Starlight output.

ariel@ufsc
Created on 05/30/2017

'''

#TODO: reading output files in each function may be a bad idea, I should
#create a quick_plot function that reads the output, and other functions
#that take an output dictionary as input.


import numpy as np
import matplotlib.pyplot as plt
from starlight_toolkit.tables import read_output_table


def plot_spec(out_file, ax=None, plot_obs=True, plot_error=False
, plot_labels=True, obs_color='k', syn_color='b', w0_color='y'
, clip_color='r', syn_label='Fitted Spectrum'):
    '''
    TODO: Add marks for clipped wavelengths, documentation.

    Parameters
    ----------
        out_file : file, string
                Name of the output file to be plotted;

        ax : matplotlib axis
                Axis to draw plot;

        plot_obs : boolean
                If True, observed spectrum will be plotted;

    '''

    if ax==None:
        ax = plt.gca()

    out = read_output_table(out_file)

    l_obs, f_obs, f_syn, f_wei = out['spectra']['l_obs'], \
    out['spectra']['f_obs'], out['spectra']['f_syn'], out['spectra']['f_wei']

    w0 = out['spectra']['f_wei'] <= 0

    clipped = out['spectra']['f_wei'] == -1.0

    error = 1/(f_wei[~w0]**2)

    if plot_obs==True:

        f_obs_masked = np.ma.masked_array(data=f_obs, mask=w0)
        f_w0         = np.ma.masked_array(data=f_obs, mask=~w0)

        ax.plot(l_obs, f_obs, color=obs_color, lw=0.5, label='Observed Spectrum')
        ax.plot(l_obs, f_w0, color=w0_color, lw=0.5, label=r'$w_\lambda=0$')

        if clipped.sum() > 0:
            ax.plot(l_obs[clipped], f_obs[clipped], color=clip_color
            , marker='x', label=r'Clipped')

    if plot_error==True:
        ax.plot(l_obs[~w0], error, '--r', label=r'Error')

    if plot_labels==True:
        ax.set_ylabel(r'$F_\lambda/F_{\lambda0}$', fontsize=15)
        ax.set_xlabel(r'$\lambda\mathrm{[\AA]}$', fontsize=15)

    ax.plot(l_obs, f_syn, color=syn_color, lw=0.5, label=syn_label)

    ax.set_ylim(0, 1.3 * np.max(f_syn))
    ax.set_xlim(out['keywords']['l_ini'],out['keywords']['l_fin'])


def plot_filter(filter_file, ax=None, filter_color='k'
, filter_ls='dashed', redshift=0, scale_factor=1):
    '''
    Plots filter file.

    To plot filters in the restframe, you should provide the galaxy's redshift,
    which is set to zero as default.

    The variable scale_factor will be multiplied by the filter's transmission.

    '''

    if ax==None:
        ax = plt.gca()

    #Reading filter:
    wl_filter, T_filter = np.genfromtxt(filter_file).transpose()

    #Plotting:
    ax.plot(wl_filter/(1+redshift), T_filter*scale_factor, color=filter_color
    , linestyle=filter_ls)


def plot_residual_spec(out_file, ax=None, residual_color='k'
, plot_labels=True):

    if ax==None:
        ax = plt.gca()

    out = read_output_table(out_file)

    l_obs, f_obs, f_syn, f_wei = out['spectra']['l_obs'], \
    out['spectra']['f_obs'], out['spectra']['f_syn'], out['spectra']['f_wei']

    w0 = out['spectra']['f_wei'] <= -1

    residual = np.ma.masked_array(data=(f_obs-f_syn)/f_syn, mask=w0)

    ax.plot(l_obs, residual, lw=0.25, color=residual_color)

    if plot_labels==True:
        ax.set_xlabel(r'$\lambda[\mathrm{\AA}]$')
        ax.set_ylabel(r'Residual Spectrum')
