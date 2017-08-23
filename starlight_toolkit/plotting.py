'''

Tools for plotting Starlight output.

ariel@ufsc
Created on 05/30/2017

'''

import numpy as np
import matplotlib.pyplot as plt
from starlight_toolkit.output import read_output_file
import starlight_toolkit.post_processing as pp
import matplotlib.gridspec as gridspec


def plot_spec(out, ax=None, plot_obs=True, plot_error=True
, plot_labels=True, obs_color='k', syn_color='b', syn_lw=0.5, w0_color='y'
, clip_color='m', syn_label=r'$M_\lambda$'):
    '''

    Quick plots for Starlight output files.

    Parameters
    ----------
        out : Starlight output dictionary;

        ax  : matplotlib axis
                Axis to draw plot;

        plot_obs : boolean
                If True, observed spectrum will be plotted;

    '''

    if ax==None:
        ax = plt.gca()


    l_obs, f_obs, f_syn, f_wei = out['spectra']['l_obs'], \
    out['spectra']['f_obs'], out['spectra']['f_syn'], out['spectra']['f_wei']
    
    w0 = out['spectra']['f_wei'] <= 0

    clipped = out['spectra']['f_wei'] == -1.0

    error = np.ma.masked_array(1/f_wei, mask=w0)

    if plot_obs==True:

        f_obs_masked = np.ma.masked_array(data=f_obs, mask=w0)
        f_w0         = np.ma.masked_array(data=f_obs, mask=~w0)

        ax.plot(l_obs, f_obs, color=obs_color, lw=0.5, label='$O_\lambda$')
        ax.plot(l_obs, f_w0, color=w0_color, lw=0.5, label=r'$w_\lambda=0$')

        if clipped.sum() > 0:
            ax.scatter(l_obs, np.ma.masked_array(f_obs, mask=~clipped), color=clip_color
            , marker='.', label=r'Clipped', zorder=5)

    if plot_error==True:
        ax.plot(l_obs, error, '--r', label=r'$\epsilon_\lambda$')

    if plot_labels==True:
        ax.set_ylabel(r'$F_\lambda/F_{\lambda0}$', fontsize=11)
        ax.set_xlabel(r'$\lambda\mathrm{[\AA]}$', fontsize=11)
    
    ax.plot(l_obs, f_syn, color=syn_color, lw=syn_lw, label=syn_label)

    if out['keywords']['IsPHOcOn']==1:
        flux_mod = out['PHO_mod']['fY_mod']/out['keywords']['fobs_norm']
        flux_obs = out['PHO_mod']['fY_obs']/out['keywords']['fobs_norm']

        #Reading redshift from output file:
        z = out['keywords']['PHO_Redshift']

        #Getting rest-frame fluxes:
        flux_mod *= (1+z)
        flux_obs *= (1+z)

        #Reading wavelengths and shifting them to the rest-frame
        wl_pho = out['PHO_mod']['MeanLamb']/(1+z)
        
        #Plotting observed photometry:
        ax.plot(wl_pho, flux_obs, 'ok', markersize=3
        , label=r'$O_{\mathrm{PHO}}$')

        #Plotting modeled photometry:
        ax.plot(wl_pho, flux_mod, 'oc', markersize=3
        , label=r'$M_{\mathrm{PHO}}$')


    ax.set_ylim(0, 1.3 * np.max(f_syn))
    ax.set_xlim(out['keywords']['l_ini'],out['keywords']['l_fin'])
    
    



def plot_spec_from_file(out_file):
    '''

    Quick plots for Starlight output files.

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

    #Reading output file (a bit dirty because we have to handle exceptions):
    try:
        out = read_output_file(out_file)
    except (ValueError, IndexError, Exception):
        print "Check if the output file is ok."

    #Plotting spectra:
    plot_spec(out)


def plot_filter(filter_file, ax=None, filter_color='k'
, filter_ls='dashed', redshift=0, scale_factor=1):
    '''
    Plots filter transmission curve.

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


def plot_residual_spec(out, ax=None, residual_color='g'
, plot_labels=True):

    if ax==None:
        ax = plt.gca()
   
    l_obs, f_obs, f_syn, f_wei = out['spectra']['l_obs'], \
    out['spectra']['f_obs'], out['spectra']['f_syn'], out['spectra']['f_wei']

    w0 = out['spectra']['f_wei'] <= -1

    residual = np.ma.masked_array(data=(f_obs-f_syn)/f_syn, mask=w0)

    ax.plot(l_obs, residual, lw=0.5, color=residual_color)

    x = np.linspace(l_obs[0], l_obs[-1])
    y = np.zeros_like(x)

    plt.plot(x, y, '--k')

    if plot_labels==True:
        ax.set_xlabel(r'$\lambda[\mathrm{\AA}]$', fontsize=10)
        ax.set_ylabel(r'Residual', fontsize=10)
    

def plot_fit_complete(out, title=None):


    fig = plt.figure(figsize=(7.5,6.25))

    #Plot spectrum:
    
    gs1 = gridspec.GridSpec(3, 3)
    gs1.update(bottom=0.47, top=0.96, hspace=0.0, right=0.96)
    
    p1 = plt.subplot(gs1[0:2,:])
    
    plot_spec(out, ax=p1, plot_labels=False)
    p1.set_ylabel(r'$F_\lambda/F_{\lambda0}$', fontsize=11)
    plt.tick_params(axis='x', bottom='off'
    , labelbottom='off')

    p1.set_ylim(0, 2)


    if title == None:
        p1.set_title(out['keywords']['arq_synt'])
    else: 
        p1.set_title(title)
    
    #Create legend:
    p1.legend(frameon=False, fontsize=9)

    #Plot residual spectrum:
    p2 = plt.subplot(gs1[2, :], sharex=p1)
    plot_residual_spec(out, ax=p2)
    
    p2.set_ylim(-0.15,0.15)
    
    #Calculating and plotting SFH:
    age_base  = out['population']['popage_base']
    Z_base    = out['population']['popZ_base']
    popx      = out['population']['popx']
    popmu_ini = out['population']['popmu_ini']
    popmu_cor = out['population']['popmu_cor']

    agevec, sfh_x, csfh_x = pp.calc_sfh_x(age_base, popx)    
    agevec, sfh_m, csfh_mu = pp.calc_sfh(age_base, popmu_ini)
    
    gs2 = gridspec.GridSpec(1, 3)
    gs2.update(top=0.38, bottom=0.08, hspace=0.05, right=0.96, wspace=0.02)
    
    p3 = plt.subplot(gs2[0, 0])

    p3.plot(np.log10(agevec), csfh_x, color='b', label=r'$x$')
    p3.plot(np.log10(agevec), csfh_x, '.b')

    p3.plot(np.log10(agevec), csfh_mu, color='r', label=r'$\mu$')
    p3.plot(np.log10(agevec), csfh_mu, '.r')

    p3.grid()

    p3.set_xlabel(r'$\log t_*$', fontsize=10)
    p3.set_ylabel('Fraction', fontsize=10)
    
    p3.legend(frameon=False, fontsize=9)
    
    p3.set_ylim(0, 110)
    
    #Annotations:
    p4 = plt.subplot(gs2[0, 1:3])
    
    #Removing ticks:
    plt.tick_params(axis='both', bottom='off'
    , top='off', labelbottom='off', right='off', left='off'
    , labelleft='off')
    
    #Calculating mass and light weighted ages:
    atflux = pp.calc_atflux(age_base, popx)
    atmass = pp.calc_atmass(age_base, popmu_cor)
    Zflux  = pp.calc_meanZflux(Z_base, popx, 0.02)
    Zmass  = pp.calc_meanZmass(Z_base, popmu_cor, 0.02)
    
    #Annotating stuff:
    
    annotation_size = 9

    p4.annotate(r'$\log M_*=$ %0.2f'%np.log10(out['keywords']['Mcor_tot']),\
    (0.02, 0.025), textcoords='axes fraction', size=annotation_size)
    p4.annotate(r'$\sigma_*$=%0.2f'%out['keywords']['vd'], \
    (0.02, 0.15), textcoords='axes fraction', size=annotation_size)
    p4.annotate(r'$N_{\mathrm{Clipped}}$=%i'%out['keywords']['Ntot_clipped'], \
    (0.02, 0.3), textcoords='axes fraction', size=annotation_size)
    p4.annotate(r'$\chi^2/N_{eff}$=%0.2f'%out['keywords']['chi2/N_eff'], \
    (0.02, 0.45), textcoords='axes fraction', size=annotation_size)
    p4.annotate(r'$\chi^2_{\mathrm{OPT}}$=%0.2f'%out['keywords']['chi2_OPT'], \
    (0.02, 0.6), textcoords='axes fraction', size=annotation_size)
    p4.annotate(r'$\chi^2_{\mathrm{TOT}}$=%0.2f'%out['keywords']['chi2'], \
    (0.02, 0.75), textcoords='axes fraction', size=annotation_size)
    p4.annotate(r'$\overline{\Delta}$=%0.2f%%'%out['keywords']['adev'], \
    (0.02, 0.9), textcoords='axes fraction', size=annotation_size)
    

    p4.annotate(r'$\delta A_V$ = %0.2f'%out['keywords']['exAv'], \
    (0.3, 0.025), textcoords='axes fraction', size=annotation_size)
    p4.annotate(r'$A_V$ = %0.2f'%out['keywords']['Av'], \
    (0.3, 0.15), textcoords='axes fraction', size=annotation_size)
    p4.annotate(r'$S/N=$ %0.2f'%out['keywords']['SN_snwin'], \
    (0.3, 0.3), textcoords='axes fraction', size=annotation_size)
    p4.annotate(r'$\langle \log Z_* \rangle_M=$%0.2f'%Zmass, \
    (0.3, 0.45), textcoords='axes fraction', size=annotation_size)
    p4.annotate(r'$\langle \log Z_* \rangle_L=$%0.2f'%Zflux, \
    (0.3, 0.6), textcoords='axes fraction', size=annotation_size)
    p4.annotate(r'$\langle \log t_* \rangle_M=$%0.2f'%atmass, \
    (0.3, 0.75), textcoords='axes fraction', size=annotation_size)
    p4.annotate(r'$\langle \log t_* \rangle_L=$%0.2f'%atflux, \
    (0.3, 0.9), textcoords='axes fraction', size=annotation_size)
    
    
    
    
    if out['keywords']['IsPHOcOn']==1:
        p4.annotate(r'$\chi^2_{PHO}=%0.2f, k_{PHO}=%0.2f$'%(out['keywords']['chi2_PHO'],out['keywords']['k_PHO']), \
        (0.6, 0.6), textcoords='axes fraction', size=annotation_size)
    if out['keywords']['IsPHOcOn']==-1:
        p4.annotate('Predicting PHO', (0.6, 0.6), textcoords='axes fraction'
        , size=annotation_size)
    if out['keywords']['IsPHOcOn']==0:
        p4.annotate('PHO off', (0.6, 0.6), textcoords='axes fraction', size=annotation_size)
        

    if out['keywords']['IsQHRcOn']==1:
        p4.annotate(r'$\chi^2_{QHR}=%0.2f, k^2_{QHR}=%0.2f$'%(out['keywords']['chi2_QHR'],out['keywords']['k_QHR']), \
        (0.6, 0.75), textcoords='axes fraction', size=annotation_size)
    if out['keywords']['IsQHRcOn']==-1:
        p4.annotate('Predicting QHR', (0.6, 0.9), textcoords='axes fraction'
        , size=annotation_size)
    if out['keywords']['IsQHRcOn']==0:
        p4.annotate('QHRc off', (0.6, 0.75), textcoords='axes fraction', size=annotation_size)


    if out['keywords']['IsFIRcOn']==1:
        p4.annotate(r'$\chi^2_{FIR}=%0.2f, k^2_{FIR}=%0.2f$'%(out['keywords']['chi2_FIR'],out['keywords']['k_FIR']), \
        (0.6, 0.9), textcoords='axes fraction', size=annotation_size)
    if out['keywords']['IsFIRcOn']==-1:
        p4.annotate('Predicting FIR', (0.6, 0.9), textcoords='axes fraction'
        , size=annotation_size)
    if out['keywords']['IsFIRcOn']==0:
        p4.annotate('FIRc off', (0.6, 0.9), textcoords='axes fraction', size=annotation_size)


#additional notes: lambda0, ESM, ELR => obs and mod
    
    

def plot_fit_complete_from_file(out_file):
    try:
        out = read_output_file(out_file)
    except (ValueError, IndexError, Exception):
        print "Check if the output file is ok."

    #Plotting spectra:
    plot_fit_complete(out)

    
