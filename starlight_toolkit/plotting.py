"""

Tools for plotting Starlight output.

ariel@ufsc
Created on 05/30/2017

"""

import numpy as np
import matplotlib.pyplot as plt
import starlight_toolkit.output as stout
import starlight_toolkit.post_processing as pp
import matplotlib.gridspec as gridspec


def plot_spec(out, ax=None, plot_obs=True, plot_syn=True, plot_error=True
              , plot_labels=True, obs_color='k', syn_color='b', syn_lw=0.6, obs_lw=0.5, w0_color='y'
              , clip_color='m', flag_color='g', syn_label=r'$M_\lambda$'
              , plot_PHO=True, PHO_color='cyan', PHO_label=r'$M_l$', PHO_obs_label=r'$O_l$'):
    """
    Quick plots for Starlight output files.
    Parameters
    ----------
        out : Starlight output dictionary or path to output file;
        ax  : matplotlib axis
                Axis to draw plot;
        plot_obs : boolean
                If True, observed spectrum will be plotted;
    """

    if ax is None:
        ax = plt.gca()

    if type(out) is str:
        try:
            out = stout.read_output_file(out)
        except (ValueError, IndexError, Exception):
            print("Check if the output file is ok.")

    l_obs, f_obs, f_syn, f_wei = out['spectra']['l_obs'], \
                                 out['spectra']['f_obs'], out['spectra']['f_syn'], out['spectra']['f_wei']

    w0 = out['spectra']['f_wei'] <= 0

    clipped = out['spectra']['f_wei'] == -1.0

    flagged = out['spectra']['f_wei'] < -1.0

    error = np.ma.masked_array(1 / f_wei, mask=w0)

    if plot_obs is True:

        f_obs_masked = np.ma.masked_array(data=f_obs, mask=w0)
        f_w0 = np.ma.masked_array(data=f_obs, mask=~w0)
        f_obs_flag = np.ma.masked_array(data=f_obs, mask=~flagged)

        ax.plot(l_obs, f_obs_masked, color=obs_color, lw=obs_lw, label=r'$O_\lambda$')
        ax.plot(l_obs, f_w0, color=w0_color, lw=obs_lw, label=r'$w_\lambda=0$')
        ax.plot(l_obs, f_obs_flag, color=flag_color, lw=obs_lw, label='$Flag$', zorder=10)

        if clipped.sum() > 0:
            ax.scatter(l_obs, np.ma.masked_array(f_obs, mask=~clipped), color=clip_color,
                       marker='x', label=r'Clipped', zorder=10)

    if plot_error is True:
        ax.plot(l_obs, error, '--r', label=r'$\epsilon_\lambda$', lw=0.5)

    if plot_labels is True:
        ax.set_ylabel(r'$F_\lambda/F_{\lambda%i}$' % out['keywords']['l_norm'], fontsize=11)
        ax.set_xlabel(r'$\lambda\mathrm{[\AA]}$', fontsize=11)
    
    if plot_syn is True:
        ax.plot(l_obs, f_syn, color=syn_color, lw=syn_lw, label=syn_label, zorder=5)

    if (out['keywords']['IsPHOcOn'] != 0) & (plot_PHO == True):
        # Reading fluxes and calculating errors
        flux_mod = out['PHO']['fY_mod'] / out['keywords']['fobs_norm']
        flux_obs = out['PHO']['fY_obs'] / out['keywords']['fobs_norm']

        flux_err = flux_obs * out['PHO']['magYErr'] / np.log10(np.e)

        # Reading redshift from output file:
        z = out['keywords']['PHO_Redshift']

        # Getting rest-frame fluxes:
        flux_mod *= (1 + z)
        flux_obs *= (1 + z)

        # Reading wavelengths and shifting them to the rest-frame
        wl_pho = out['PHO']['MeanLamb'] / (1 + z)

        if out['keywords']['IsPHOcOn'] == 1:
            # Plotting observed photometry:
            ax.errorbar(wl_pho, flux_obs, flux_err, fmt='ok', ecolor='k', markersize=4,
                        label=PHO_obs_label)

        # Plotting modeled photometry:
        ax.plot(wl_pho, flux_mod, 'o', color=PHO_color, markersize=4
                , label=PHO_label, zorder=15)

    ax.set_ylim(0, 1.5 * np.max(f_syn))
    ax.set_xlim(out['keywords']['l_ini'], out['keywords']['l_fin'])


def plot_spec_simple(out, ax=None, plot_obs=True, plot_syn=True, plot_error=True,
                     plot_labels=True, obs_color='k', syn_color='b', syn_lw=0.6, obs_lw=0.5, w0_color='y',
                     clip_color='m', syn_label=r'$M_\lambda$',
                     plot_PHO=True, PHO_color='cyan', PHO_edgecolor=None,
                     PHO_label=r'$M_l$', PHO_obs_label=r'$O_l$', PHO_markersize=5):
    """
    Quick plots for Starlight output files.
    Parameters
    ----------
        ax  : matplotlib axis
        out : Starlight output dictionary or path to output file;
        plot_obs : boolean
                If True, observed spectrum will be plotted;
    """

    if ax is None:
        ax = plt.gca()

    if PHO_edgecolor is None:
        PHO_edgecolor = PHO_color

    if type(out) is str:
        try:
            out = stout.read_output_file(out)
        except (ValueError, IndexError, Exception):
            print("Check if the output file is ok.")

    l_obs, f_obs, f_syn, f_wei = out['spectra']['l_obs'], \
                                 out['spectra']['f_obs'], out['spectra']['f_syn'], out['spectra']['f_wei']

    w0 = out['spectra']['f_wei'] <= 0

    error = np.ma.masked_array(1 / f_wei, mask=w0)

    if plot_obs is True:
        f_obs_masked = np.ma.masked_array(data=f_obs, mask=w0)
        f_w0 = np.ma.masked_array(data=f_obs, mask=~w0)

        ax.plot(l_obs, f_obs_masked, color=obs_color, lw=obs_lw, label=r'$O_\lambda$')
        ax.plot(l_obs, f_w0, color=w0_color, lw=obs_lw, label=r'$w_\lambda=0$')

    if plot_error is True:
        ax.plot(l_obs, error, '--r', label=r'$\epsilon_\lambda$', lw=0.5)

    if plot_labels is True:
        ax.set_ylabel(r'$F_\lambda/F_{\lambda%i}$' % out['keywords']['l_norm'], fontsize=11)
        ax.set_xlabel(r'$\lambda\mathrm{[\AA]}$', fontsize=11)

    if plot_syn is True:
        ax.plot(l_obs, f_syn, color=syn_color, lw=syn_lw, label=syn_label, zorder=5)

    if (out['keywords']['IsPHOcOn'] != 0) & (plot_PHO == True):
        flux_mod = out['PHO']['fY_mod'] / out['keywords']['fobs_norm']
        flux_obs = out['PHO']['fY_obs'] / out['keywords']['fobs_norm']

        flux_err = flux_obs * out['PHO']['magYErr'] / np.log10(np.e)

        # Reading redshift from output file:
        z = out['keywords']['PHO_Redshift']

        # Getting rest-frame fluxes:
        flux_mod *= (1 + z)
        flux_obs *= (1 + z)

        # Reading wavelengths and shifting them to the rest-frame
        wl_pho = out['PHO']['MeanLamb'] / (1 + z)

        if out['keywords']['IsPHOcOn'] == 1:
            # Plotting observed photometry:
            ax.errorbar(wl_pho, flux_obs, flux_err, fmt='ok', ecolor='k', markersize=PHO_markersize
                        , label=PHO_obs_label, zorder=14)

        # Plotting modeled photometry:
        ax.plot(wl_pho, flux_mod, 'o', color=PHO_color, markersize=PHO_markersize
                , label=PHO_label, zorder=15, markeredgecolor=PHO_edgecolor, markeredgewidth=1)

    ax.set_ylim(0, 1.5 * np.max(f_syn))
    ax.set_xlim(out['keywords']['l_ini'], out['keywords']['l_fin'])



def plot_filter(filter_file, ax=None, filter_color='k'
                , filter_ls='dashed', redshift=0, scale_factor=1):
    """
    Plots filter transmission curve.

    To plot filters in the restframe, you should provide the galaxy's redshift,
    which is set to zero as default.

    The variable scale_factor will be multiplied by the filter's transmission.

    """

    if ax is None:
        ax = plt.gca()

    # Reading filter:
    wl_filter, T_filter = np.genfromtxt(filter_file).transpose()

    # Plotting:
    ax.plot(wl_filter / (1 + redshift), T_filter * scale_factor, color=filter_color
            , linestyle=filter_ls)


def plot_residual_spec(out, ax=None, residual_color='g'
                       , plot_labels=True):
    if ax is None:
        ax = plt.gca()

    l_obs, f_obs, f_syn, f_wei = out['spectra']['l_obs'], \
                                 out['spectra']['f_obs'], out['spectra']['f_syn'], out['spectra']['f_wei']

    w0 = out['spectra']['f_wei'] <= -1

    residual = np.ma.masked_array(data=(f_obs - f_syn) / f_syn, mask=w0)

    ax.plot(l_obs, residual, lw=0.3, color=residual_color)

    x = np.linspace(l_obs[0], l_obs[-1])
    y = np.zeros_like(x)

    ax.plot(x, y, '--k')

    if plot_labels is True:
        ax.set_xlabel(r'$\lambda[\mathrm{\AA}]$', fontsize=10)
        ax.set_ylabel(r'Residual', fontsize=10)


def plot_sfh(out, ax=None, plot_axlabels=True):
    """

    Plot's STARLIGHT's star-formation histories

    Parameters
    ----------

    out: dictionary
        STARLIGHT output dictionary

    ax: matplotlib axis
        axis on which to plot

    plot_axlabels: boolean
        wether to plot labels or not

    """

    if ax is None:
        ax = ax.gca()

    age_base = out['population']['popage_base']
    popx = out['population']['popx']
    popmu_ini = out['population']['popmu_ini']

    agevec, sfh_x, csfh_x = pp.calc_sfh_x(age_base, popx)
    agevec, sfh_m, csfh_mu = pp.calc_sfh(age_base, popmu_ini)

    ax.plot(np.log10(agevec), csfh_x, color='b', label=r'$x$')
    ax.plot(np.log10(agevec), csfh_x, '.b')

    ax.plot(np.log10(agevec), csfh_mu, color='r', label=r'$\mu$')
    ax.plot(np.log10(agevec), csfh_mu, '.r')

    if out['keywords']['IsQHRcOn'] == 1:
        popQHR = out['popQHR']['QHeff_Perc']
        cQHRvec = pp.calc_QHRpop_x(age_base, popQHR)
        ax.plot(np.log10(agevec), cQHRvec, color='m', label=r'$Q_H$')
        ax.plot(np.log10(agevec), cQHRvec, '.m')

    if plot_axlabels is True:
        ax.set_xlabel(r'$\log \; t_*$', fontsize=10)
        ax.set_ylabel('Cumulative Fraction [%]', fontsize=10)


def plot_fit_complete(out, title=None, figsize=(7.75, 6.5)
                      , out_fig_fname=None, out_dpi=None
                      , legend_loc='best'):
    if type(out) is str:
        try:
            out = stout.read_output_file(out)
        except (ValueError, IndexError, Exception):
            print("Check if the output file is ok.")

    fig = plt.figure(figsize=figsize)

    # Plot spectrum:

    gs1 = gridspec.GridSpec(3, 3)
    gs1.update(bottom=0.47, top=0.96, hspace=0.05, right=0.96)

    gs2 = gridspec.GridSpec(1, 3)
    gs2.update(top=0.38, bottom=0.08, hspace=0.05, right=0.96, wspace=0.02)

    p1 = plt.subplot(gs1[0:2, :])
    p2 = plt.subplot(gs1[2, :], sharex=p1)
    p3 = plt.subplot(gs2[0, 0])
    p4 = plt.subplot(gs2[0, 1:3])

    plot_spec(out, ax=p1, plot_labels=False)
    p1.set_ylabel(r'$F_\lambda/F_{\lambda%i}$' % out['keywords']['l_norm']
                  , fontsize=11)
    plt.tick_params(axis='x', bottom='off'
                    , labelbottom='off')

    p1.set_ylim(0, 2.25)

    if title is None:
        p1.set_title(out['keywords']['arq_synt'])
    else:
        p1.set_title(title)

    # Create legend:
    p1.legend(frameon=False, fontsize=9, ncol=2, loc=legend_loc)

    # Plot residual spectrum:
    plot_residual_spec(out, ax=p2)

    p2.set_ylim(-0.15, 0.15)

    # Plotting SFH:

    plot_sfh(out, p3)

    p3.set_xticks([6, 7, 8, 9, 10])

    p3.legend(frameon=False, fontsize=9)

    p3.set_ylim(0, 110)
    p3.set_xlim(5.8, 11)

    # Annotations:

    # Removing ticks:
    p4.set_xticks([])
    p4.set_yticks([])

    # Calculating mass and light weighted ages:

    if 'popage_base_upp' in out['population'].keys():

        age_base = out['population']['popage_base']
        age_base_upp = out['population']['popage_base_upp']
        Z_base = out['population']['popZ_base']
        popx = out['population']['popx']
        popmu_cor = out['population']['popmu_cor']

        atflux = pp.calc_atflux(age_base, popx, age_base_upp)
        atmass = pp.calc_atmass(age_base, popmu_cor, age_base_upp)
        Zflux = pp.calc_aZflux(Z_base, popx, 0.02)
        Zmass = pp.calc_aZmass(Z_base, popmu_cor, 0.02)

    else:

        age_base = out['population']['popage_base']
        Z_base = out['population']['popZ_base']
        popx = out['population']['popx']
        popmu_cor = out['population']['popmu_cor']

        atflux = pp.calc_atflux(age_base=age_base, popx=popx)
        atmass = pp.calc_atmass(age_base=age_base, popmu=popmu_cor)
        Zflux = pp.calc_aZflux(Z_base, popx, 0.02)
        Zmass = pp.calc_aZmass(Z_base, popmu_cor, 0.02)

    # Annotating stuff:

    annotation_size = 9

    p4.grid(False)

    p4.annotate(r'$\log M_*=$ %0.2f' % np.log10(out['keywords']['Mcor_tot']),
                (0.02, 0.025), size=annotation_size)
    p4.annotate(r'$\sigma_*$=%0.2f' % out['keywords']['vd'],
                (0.02, 0.15), size=annotation_size)
    p4.annotate(r'$N_{\mathrm{Clipped}}$=%i' % out['keywords']['Ntot_clipped'],
                (0.02, 0.3), size=annotation_size)
    p4.annotate(r'$\chi^2/N_{eff}$=%0.2f' % out['keywords']['chi2/N_eff'],
                (0.02, 0.45), size=annotation_size)
    try:
        p4.annotate(r'$\chi^2_{\mathrm{OPT}}$=%0.2f' % out['keywords']['chi2_OPT'],
                    (0.02, 0.6), size=annotation_size)
    except:
        pass
    p4.annotate(r'$\chi^2_{\mathrm{TOT}}$=%0.2f' % out['keywords']['chi2'],
                (0.02, 0.75), size=annotation_size)
    p4.annotate(r'$\overline{\Delta}$=%0.2f%%' % out['keywords']['adev'],
                (0.02, 0.9), size=annotation_size)

    p4.annotate(r'$A_V^Y$ = %0.2f' % (out['keywords']['AV'] + out['keywords']['exAV']),
                (0.27, 0.025), size=annotation_size)
    p4.annotate(r'$A_V$ = %0.2f' % out['keywords']['AV'],
                (0.27, 0.15), size=annotation_size)
    p4.annotate(r'$S/N=$ %0.2f' % out['keywords']['SN_snwin'],
                (0.27, 0.3), size=annotation_size)
    p4.annotate(r'$\langle \log Z_{*} \rangle_M=$%0.2f' % Zmass,
                (0.27, 0.45), size=annotation_size)
    p4.annotate(r'$\langle \log Z_{*} \rangle_L=$%0.2f' % Zflux,
                (0.27, 0.6), size=annotation_size)
    p4.annotate(r'$\langle \log t_{*} \rangle_M=$%0.2f' % atmass,
                (0.27, 0.75), size=annotation_size)
    p4.annotate(r'$\langle \log t_{*} \rangle_L=$%0.2f' % atflux,
                (0.27, 0.9), size=annotation_size)

    p4.annotate('ESM=%s' % out['keywords']['ETC_ESM'],
                (0.51, 0.9), size=annotation_size)

    if out['keywords']['IsPHOcOn'] == 1:
        p4.annotate(r'$\chi^2_{PHO}=%0.2f, k_{PHO}=%0.2f$'\
                    % (out['keywords']['chi2_PHO'], out['keywords']['k_PHO']),
                    (0.51, 0.45), size=annotation_size)
    if out['keywords']['IsPHOcOn'] == -1:
        p4.annotate('Predicting PHO', (0.51, 0.45), size=annotation_size)
    if out['keywords']['IsPHOcOn'] == 0:
        p4.annotate('PHO off', (0.51, 0.45), size=annotation_size)

    if out['keywords']['IsELROn'] == 0:
        p4.annotate(r'$x(A_V^Y) = {} \%$'.format(out['keywords']['x(exAV>0)']),
                    (0.51, 0.3), size=annotation_size)

    if out['keywords']['IsQHRcOn'] == 1:
        p4.annotate(r'$\chi^2_{QHR}=%0.2f, k_{QHR}=%0.2f$'\
                    % (out['keywords']['chi2_QHR'], out['keywords']['k_QHR']),
                    (0.51, 0.6), size=annotation_size)
        if out['keywords']['IsELROn'] == 1:
            exHalpha = np.sum(out['popQHR']['Y_Perc_Line0'][out['population']['popexAV_flag'] == 1])
            epsilonQ = out['ELR']['Err_logR'] / np.sqrt(out['keywords']['k_QHR'])

            p4.annotate(r'$x(A_V^Y) = {} \%, x(A_V^Y)H_\alpha$ = {} %'.format(out['keywords']['x(exAV>0)'], exHalpha),
                        (0.51, 0.15), size=annotation_size)

            p4.annotate(r'$AV_{neb}=%0.2f$' \
                        % (out['keywords']['AV_neb']), (0.51, 0.025), size=annotation_size)

            p4.annotate(r'$\Delta\logR=%0.3f, \epsilon^{eff}_R=%1.2e$' \
                        % (out['ELR']['logR_obs'] - out['ELR']['logR_mod'], epsilonQ),
                        (0.51, 0.3), size=annotation_size)

        if out['keywords']['IsELROn'] == -1:
            p4.annotate('Predicting ELR', (0.51, 0.3), size=annotation_size)
        if out['keywords']['IsELROn'] == 0:
            p4.annotate('ELR off', (0.51, 0.3), size=annotation_size)
    if out['keywords']['IsQHRcOn'] == -1:
        p4.annotate('Predicting QHR', (0.51, 0.6), size=annotation_size)
    if out['keywords']['IsQHRcOn'] == 0:
        p4.annotate('QHRc off', (0.51, 0.6), size=annotation_size)

    if out['keywords']['IsFIRcOn'] == 1:
        p4.annotate(r'$\chi^2_{FIR}=%0.2f, k^2_{FIR}=%0.2f$' \
                    % (out['keywords']['chi2_FIR'], out['keywords']['k_FIR']), \
                    (0.51, 0.75), size=annotation_size)
    if out['keywords']['IsFIRcOn'] == -1:
        p4.annotate(r'$\logL_{FIR}^{mod}=%0.2f, FIR/BOL=%0.2f$' \
                    % (out['keywords']['FIR_logLFIR_mod'], out['keywords']['FIR_BOL_Ratio'])
                    , (0.51, 0.75), size=annotation_size)
    if out['keywords']['IsFIRcOn'] == 0:
        p4.annotate('FIRc off', (0.51, 0.75), size=annotation_size)

    if out_fig_fname is not None:
        plt.savefig(out_fig_fname, dpi=out_dpi)
        plt.close()

    return p1, p2, p3, p4
