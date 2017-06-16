'''
Created on 22/06/2015

@author: andre
'''

from . import flags
import numpy as np


__all__ = ['resample_spectra', 'reshape_cube', 'find_nearest_index',
           'interp1d_spectra', 'gaussian1d_spectra', 'gen_rebin', 'bin_edges', 'hist_resample',
           'age_smoothing_kernel', 'light2mass_ini', 'interp_age', 'vac2air', 'get_subset_slices']


def ReSamplingMatrixNonUniform(lorig, lresam, extrap=False):
    '''
    Compute resampling matrix R_o2r, useful to convert a spectrum sampled at 
    wavelengths lorig to a new grid lresamp. Here, there is no necessity to have constant gris as on :py:func:`ReSamplingMatrix`.
    Input arrays lorig and lresamp are the bin centres of the original and final lambda-grids.
    ResampMat is a Nlresamp x Nlorig matrix, which applied to a vector F_o (with Nlorig entries) returns
    a Nlresamp elements long vector F_r (the resampled spectrum):

        [[ResampMat]] [F_o] = [F_r]

    Warning! lorig and lresam MUST be on ascending order!


    Parameters
    ----------
    lorig : array_like
            Original spectrum lambda array.

    lresam : array_like
             Spectrum lambda array in which the spectrum should be sampled.        

    extrap : boolean, optional
           Extrapolate values, i.e., values for lresam < lorig[0]  are set to match lorig[0] and
                                     values for lresam > lorig[-1] are set to match lorig[-1].


    Returns
    -------
    ResampMat : array_like
                Resample matrix. 

    Examples
    --------
    >>> lorig = np.linspace(3400, 8900, 9000) * 1.001
    >>> lresam = np.linspace(3400, 8900, 5000)
    >>> forig = np.random.normal(size=len(lorig))**2
    >>> matrix = slut.ReSamplingMatrixNonUniform(lorig, lresam)
    >>> fresam = np.dot(matrix, forig)
    >>> print np.trapz(forig, lorig), np.trapz(fresam, lresam)
    '''

    # Init ResampMatrix
    matrix = np.zeros((len(lresam), len(lorig)))

    # Define lambda ranges (low, upp) for original and resampled.
    lo_low = np.zeros(len(lorig))
    lo_low[1:] = (lorig[1:] + lorig[:-1]) / 2
    lo_low[0] = lorig[0] - (lorig[1] - lorig[0]) / 2

    lo_upp = np.zeros(len(lorig))
    lo_upp[:-1] = lo_low[1:]
    lo_upp[-1] = lorig[-1] + (lorig[-1] - lorig[-2]) / 2

    lr_low = np.zeros(len(lresam))
    lr_low[1:] = (lresam[1:] + lresam[:-1]) / 2
    lr_low[0] = lresam[0] - (lresam[1] - lresam[0]) / 2

    lr_upp = np.zeros(len(lresam))
    lr_upp[:-1] = lr_low[1:]
    lr_upp[-1] = lresam[-1] + (lresam[-1] - lresam[-2]) / 2

    # Find out if the original array is simply a subset of the resampled.
    # If it is, this routine already returns the correct resampling matrix
    subset = lresam_subset_lorig(lresam, lorig, matrix)

    if (not subset):

        # Iterate over resampled lresam vector
        for i_r in range(len(lresam)):

            # Find in which bins lresam bin within lorig bin
            bins_resam = np.where(
                (lr_low[i_r] < lo_upp) & (lr_upp[i_r] > lo_low))[0]

            # On these bins, eval fraction of resamled bin is within original
            # bin.
            for i_o in bins_resam:

                aux = 0

                d_lr = lr_upp[i_r] - lr_low[i_r]
                d_lo = lo_upp[i_o] - lo_low[i_o]
                d_ir = lo_upp[i_o] - lr_low[i_r]  # common section on the right
                d_il = lr_upp[i_r] - lo_low[i_o]  # common section on the left

                # Case 1: resampling window is smaller than or equal to the original window.
                # This is where the bug was: if an original bin is all inside the resampled bin, then
                # all flux should go into it, not then d_lr/d_lo fraction.
                # --Natalia@IoA - 21/12/2012
                if (lr_low[i_r] >= lo_low[i_o]) & (lr_upp[i_r] <= lo_upp[i_o]):
                    aux += 1.

                # Case 2: resampling window is larger than the original window.
                if (lr_low[i_r] < lo_low[i_o]) & (lr_upp[i_r] > lo_upp[i_o]):
                    aux += d_lo / d_lr

                # Case 3: resampling window is on the right of the original
                # window.
                if (lr_low[i_r] >= lo_low[i_o]) & (lr_upp[i_r] > lo_upp[i_o]):
                    aux += d_ir / d_lr

                # Case 4: resampling window is on the left of the original
                # window.
                if (lr_low[i_r] < lo_low[i_o]) & (lr_upp[i_r] <= lo_upp[i_o]):
                    aux += d_il / d_lr

                matrix[i_r, i_o] += aux

    # Fix matrix to be exactly = 1 ==> TO THINK
    # print np.sum(matrix), np.sum(lo_upp - lo_low), (lr_upp - lr_low).shape

    # Fix extremes: extrapolate if needed
    if (extrap):
        extrap_ResampMat(matrix, lo_low, lo_upp, lr_low, lr_upp)

    return matrix


def lresam_subset_lorig(lresam, lorig, ResampMat__ro):
    '''
    Finds out if lorig is a subset of lresam. If so,
    create a `diagonal' resampling matrix (in place).
    '''

    subset = False
    tol = 1.e-5

    ff = (lresam >= lorig[0]) & (lresam <= lorig[-1])
    lresam_sub = lresam[ff]

    N_lo = len(lorig)

    if (len(lresam_sub) == N_lo):
        Nok = (abs(lorig - lresam_sub) < tol).sum()

        if (Nok == N_lo):
            subset = True
            ilow = np.where(ff)[0][0]
            iupp = np.where(ff)[0][-1]

            io_diag = np.arange(len(lorig))
            ir_diag = np.arange(ilow, iupp + 1)
            ResampMat__ro[ir_diag, io_diag] = 1.

    return subset


def extrap_ResampMat(ResampMat__ro, lo_low, lo_upp, lr_low, lr_upp):
    '''
    Extrapolate resampling matrix on the borders, i.e.,
    by making it copy the first and the last bins.

    This modifies ResampMat__ro in place.
    '''

    bins_extrapl = np.where((lr_low < lo_low[0]))[0]
    bins_extrapr = np.where((lr_upp > lo_upp[-1]))[0]

    if (len(bins_extrapl) > 0) & (len(bins_extrapr) > 0):
        io_extrapl = np.where((lo_low >= lr_low[bins_extrapl[0]]))[0][0]
        io_extrapr = np.where((lo_upp <= lr_upp[bins_extrapr[0]]))[0][-1]

        ResampMat__ro[bins_extrapl, io_extrapl] = 1.
        ResampMat__ro[bins_extrapr, io_extrapr] = 1.


def resample_spectra(l_orig, l_resam, f_obs, f_err, badpix):
    '''
    Resample IFS wavelength-wise.

    Parameters
    ----------
    l_orig : array
        Original wavelength base of ``spectra``.

    l_resam : array
        Destination wavelength base.

    f_obs : array
        Observed spectra to be resampled.

    f_err : array
        Error spectra to be resampled.

    badpix : array(dtype=bool)
        bad pixel spectra to be resampled.

    Returns
    -------
    spectra_resam : array
        Spectra resampled to ``l_resam``.
    '''
    R = ReSamplingMatrixNonUniform(l_orig, l_resam)
    f_obs = np.tensordot(R, f_obs, (1, 0))
    f_err = np.sqrt(np.tensordot(R, f_err**2, (1, 0)))
    badpix = np.tensordot(R, badpix.astype('float64'), (1, 0)) > 0.0
    f_flag = np.zeros_like(f_obs, dtype='int32')
    out_of_range = (l_resam < l_orig[0]) | (l_resam > l_orig[-1])
    f_flag[out_of_range] |= flags.no_data
    f_flag[badpix] |= flags.bad_pix
    return f_obs, f_err, f_flag


def reshape_cube(f_obs, f_err, badpix, center, new_shape):
    '''
    Reshape IFS into a new spatial shape, putting the given
    photometric center at the center of the new IFS.
    Flag the newly added pixels as ``no_data``.

    Parameters
    ----------
    f_obs : array
        Flux, will remain unchanged.

    f_err : array
        Flux, will remain unchanged.

    badpix : array
        bad pixel spectra to be resampled.

    center : tuple (l, y, x)
        IFS center, or reference pixel.

    new_shape : tuple (Nl, Ny, Nx)
        New shape.

    Returns
    -------
    res_f_obs : array
        Reshaped ``f_obs``.

    res_f_flag : array
        Reshaped ``f_flag``

    res_center : tuple (l, y, x)
        New center or reference pixel of image.
    '''
    shape = f_obs.shape
    y_axis = 1
    x_axis = 2
    Nx = shape[x_axis]
    Ny = shape[y_axis]
    Nx_r = new_shape[x_axis]
    Ny_r = new_shape[y_axis]
    xc = center[x_axis]
    yc = center[y_axis]
    xc_r = int(Nx_r / 2)
    yc_r = int(Ny_r / 2)

    xi_r = xc_r - xc
    if xi_r < 0:
        xi = -xi_r
        xi_r = 0
    else:
        xi = 0

    yi_r = yc_r - yc
    if yi_r < 0:
        yi = -yi_r
        yi_r = 0
    else:
        yi = 0

    xf_r = xi_r + Nx
    if xf_r > Nx_r:
        xf = xi + Nx_r
        xf_r = Nx_r
    else:
        xf = Nx

    yf_r = yi_r + Ny
    if yf_r > Ny_r:
        yf = yi + Ny_r
        yf_r = Ny_r
    else:
        yf = Ny

    res_f_obs = np.zeros((new_shape))
    res_f_err = np.zeros((new_shape))
    res_badpix = np.ones_like(res_f_obs, dtype='bool')
    res_f_obs[:, yi_r:yf_r, xi_r:xf_r] = f_obs[:, yi:yf, xi:xf]
    res_f_err[:, yi_r:yf_r, xi_r:xf_r] = f_err[:, yi:yf, xi:xf]
    res_badpix[:, yi_r:yf_r, xi_r:xf_r] = badpix[:, yi:yf, xi:xf]

    res_center = (center[0], yc_r, xc_r)

    return res_f_obs, res_f_err, res_badpix, res_center


def fwhm2sigma(fwhm):
    return fwhm / (2 * np.sqrt(2 * np.log(2)))


def interp1d_spectra(l, flux, flags=None):
    '''
    Interpolate linearly 1-d spectra to fill all the gaps and
    extend the limits (copies the boundaries).

    Parameters
    ----------
    l : array
        Wavelength array.

    flux : array
        [Masked] array with gaps.

    flags : array, optional
        Array with flags as integers. Bad values
        are greater than zero. Must be set if
        ``flux`` is not a masked array. Default: ``None``.

    Returns
    -------
    flux_interp : array
        The same as ``flux``, wit gaps replaced by linear interpolation.
    '''
    if not isinstance(flux, np.ma.MaskedArray):
        if flags is None:
            raise Exception('flux must be a masked array if flags is not set.')
        flux = np.ma.array(flux, mask=flags > 0)
    lc = l[~flux.mask]
    fc = flux.compressed()
    return np.interp(l, lc, fc)


def gaussian1d_spectra(fwhm, l, flux, flags=None):
    '''
    Filter 1-d spectra using a Gaussian kernel. Interpolate linearly
    the flagged wavelengths before applying the filter.

    Parameters
    ----------
    fwhm : float
        Full width at half maximum of the gaussian.

    l : array
        Wavelength array.

    flux : array
        [Masked] array with gaps.

    flags : array, optional
        Array with flags as integers. Bad values
        are greater than zero. Must be set if
        ``flux`` is not a masked array. Default: ``None``.

    Returns
    -------
    flux_gauss : array
        The same as ``flux``, interpolated by a gaussian kernel.
    '''
    flux_i = interp1d_spectra(l, flux, flags)
    dl = (l[1] - l[0])
    sig = fwhm2sigma(fwhm)
    l_gauss = np.arange(
        min(-2 * dl, -5 * sig), max(2 * dl + 1, 5 * sig + 1), dl)
    amp = 1.0 / (sig * np.sqrt(2.0 * np.pi))
    gauss = dl * amp * np.exp(-0.5 * (l_gauss / sig)**2)
    return np.convolve(flux_i,  gauss, 'same')


def gen_rebin(a, e, bin_e, mean=True):
    '''
    Rebinning function. Given the value array `a`, with generic
    positions `e` such that `e.shape == a.shape`, return the sum
    or the mean values of `a` inside the bins defined by `bin_e`.

    Parameters
    ----------
    a : array like
        The array values to be rebinned.

    e : array like
        Generic positions of the values in `a`.

    bin_e : array like
        Bins of `e` to be used in the rebinning.

    mean : boolean
        Divide the sum by the number of points inside the bin.
        Defaults to `True`.

    Returns
    -------
    a_e : array
        An array of length `len(bin_e)-1`, containing the sum or
        mean values of `a` inside each bin.

    Examples
    --------

    TODO: Add examples for gen_rebin.
    '''
    a_e = np.histogram(e.ravel(), bin_e, weights=a.ravel())[0]
    if mean:
        N = np.histogram(e.ravel(), bin_e)[0]
        mask = N > 0
        a_e[mask] /= N[mask]
    return a_e


def bin_edges(bin_center):
    '''
    Computes the bin edges as the bissection
    of the bins, expanding the borders accordingly.

    Parameters
    ----------
    bin_center : array
        Bin centers.

    Returns
    -------
    bin_edges : array
        Bin edges, with length ``len(bin_center) + 1``.
    '''
    bin_edges = np.empty(len(bin_center) + 1)
    bin_edges[1:-1] = (bin_center[:-1] + bin_center[1:]) / 2.0
    bin_edges[0] = bin_center[0] - (bin_center[1] - bin_center[0]) / 2.0
    bin_edges[-1] = bin_center[-1] + (bin_center[-1] - bin_center[-2]) / 2.0
    return bin_edges


def hist_resample(bins_o, bins_r, v, density=False):
    '''
    Resample an histogram into another set of bins.

    Parameters
    ----------
    bins_o : array
        The original bin edges.

    bins_r : array
        Bin edges to resample.

    v : array
        The heights of the histogram.
        Must be of length ``len(bins_o) - 1``.
    density : boolean, optional
        ``v`` is a density, not a histogram.
        Default: ``False``.

    Returns
    -------
    v_r : array
        The heights of the resampled histogram
        (or resampled densities if ``density=True``),
        with length ``len(bins_r) - 1``.
        If ``bins_r`` completely overlaps ``bins_o``,
        then the sums ``v.sum()`` and ``v_r.sum()``
        will be equal (within numeric error).
    '''
    n_v_r = len(bins_r) - 1
    v_r = np.zeros(n_v_r, dtype=v.dtype)
    i_0 = 0
    j = 0
    lbo = bins_o[i_0]
    rbo = bins_o[i_0 + 1]
    lbr = bins_r[j]
    rbr = bins_r[j + 1]
    dbr = rbr - lbr

    if lbo < lbr:
        # Skip leading original bins that do not overlap the resampled bins.
        while lbr > rbo:
            i_0 += 1
            lbo = bins_o[i_0]
            rbo = bins_o[i_0 + 1]
            dbo = rbo - lbo
            # print 'dbo', dbo
        frac = (rbo - lbr) / dbo
        if density:
            v_r[j] += frac * v[i_0] * dbo
        else:
            v_r[j] += frac * v[i_0]
        i_0 += 1
        last_edge = rbo
        v_r[j] /= dbr
    else:
        # Skip leading resampled bins that do not overlap the original bins.
        while rbr < lbo:
            j += 1
            lbr = bins_r[j]
            rbr = bins_r[j + 1]
        last_edge = rbr

    for i in xrange(i_0, len(v)):
        rbo = bins_o[i + 1]
        lbo = bins_o[i]
        dbo = rbo - lbo
        if rbo < rbr:
            # This original bin is entirely contained in the resampled bin.
            if density:
                v_r[j] += v[i] * dbo
            else:
                v_r[j] += v[i]
            last_edge = rbo
            continue
        # print 'dbo', dbo
        while rbr < rbo:
            # Compute the slices of the original bin going to the resampled
            # bins.
            frac = (rbr - last_edge) / dbo
            if density:
                v_r[j] += frac * v[i] * dbo
                v_r[j] /= dbr
            else:
                v_r[j] += frac * v[i]
            last_edge = rbr
            j += 1
            if j >= n_v_r:
                break
            lbr = bins_r[j]
            rbr = bins_r[j + 1]
            dbr = rbr - lbr
        if j >= n_v_r:
            break
        # Compute the last slice of the original bin.
        frac = (rbo - last_edge) / dbo
        if density:
            # print 'dbr', dbr
            v_r[j] += frac * v[i] * dbo
        else:
            v_r[j] += frac * v[i]
        last_edge = rbo
    return v_r


def age_smoothing_kernel(log_age_base, log_tc, logtc_FWHM=0.1):
    '''
    Given an array of logAgeBase = log10(base ages), and an array logTc of "continuous" log ages 
    (presumably uniformly spaced), computes the smoothing kernel s_bc, which resamples and smooths
    any function X(logAgeBase) to Y(logTc). 
    If X = X_b, where b is the logAgeBase index, and Y = Y_c, with c as the index in logTc, then 

    Y_c = sum-over-all-b's of X_b s_bc

    gives the smoothed & resampled version of X. 
    The smoothing is performed in log-time, with gaussian function of FWHM = logtc_FWHM. 
    Conservation of X is ensured (ie, X_b = sum-over-all-c's of Y_c s_bc).

    Notice that logTc and logtc_FWHM are given default values, in case of lazy function callers...
    However, I do NOT know how to call the function using the default logTc but setting a FWHM different from default!!!
    ???Andre????

    Input:  logAgeBase = array of log base ages [in log yr]
            logTc = array of log "continous" ages [in log yr]
            logtc_FWHM = width of smoothing filter [in dex]
    Output: s__bc = (len(logAgeBase) , len(logTc)) matrix [adimensional]

    ElCid@Sanchica - 18/Mar/2012
    '''
    s__bc = np.zeros((len(log_age_base), len(log_tc)))
    logtc_sig = logtc_FWHM / (np.sqrt(8 * np.log(2)))
    for i_b, a_b in enumerate(log_age_base):
        aux1 = np.exp(-0.5 * ((log_tc - a_b) / logtc_sig)**2)
        s__bc[i_b, :] = aux1 / aux1.sum()
    return s__bc


def interp_age(prop, log_age_base, log_age_interp):
    '''
    Interpolate linearly Mstars or fbase_norm in log time.

    Parameters
    ----------
    prop : array
        Array containing ``Mstars`` or ``fbase_norm``.

    log_age_base : array
        The age base, the same legth as the first
        dimension of ``prop``.

    log_age_interp : array
        The age to which interpolate ``prop``. The
        returned value will have the same length in
        the first dimension.

    Returns
    -------
    propI_interp : array
        The same ``prop``, interpolated to ``logAgeInterp``.
    '''
    n_met = prop.shape[1]
    n_age_interp = len(log_age_interp)
    prop_interp = np.empty((n_age_interp, n_met), dtype='>f8')
    for z in xrange(n_met):
        prop_interp[:, z] = np.interp(log_age_interp, log_age_base, prop[:, z])
    return prop_interp


def light2mass_ini(popx, fbase_norm, Lobs_norm, q_norm, A_V):
    '''
    Compute the initial mass from popx (and other parameters).
    The most important thing to remember is that popx (actually luminosity)
    must be "dereddened" before converting it to mass using fbase_norm
    (popmu_ini).

    Based on the subroutine ConvertLight2Mass (starlight source code).

    Parameters
    ----------
    popx: array
        The light fractions (in percents).

    fbase_norm: array
        Light to mass ratio for each population.

    Lobs_norm : array
        Luminosity norm of ``popx``.

    q_norm : float 
        Ratio between the extinction in l_norm (where Lobs_norm
        is calculated) and ``A_V``.

    A_V : array 
        Extinction in V band.

    Returns
    -------
    Mini : array
        The initial mass in each population.

    '''
    Lobs = popx / 100.0
    Lobs *= Lobs_norm
    Lobs *= 10.0**(0.4 * q_norm * A_V)

    Mini = Lobs / fbase_norm[..., np.newaxis, np.newaxis]
    return Mini


def vac2air(wave):
    '''
    As given in Morton (1991, ApJS, 77, 119).
    '''
    wave_air = wave / \
        (1.0 + 2.735182e-4 + 131.4182 / wave**2 + 2.76249e8 / wave**4)
    return wave_air


def find_nearest_index(array, value):
    '''
    Return the array index that is closest to the valued provided. Note that
    this is intended for use with coordinates array.
    '''
    if np.isscalar(value):
        value = np.array([value])
    else:
        value = np.array(value)
    idx = (np.abs(array - value[:, np.newaxis])).argmin(axis=1)
    if len(idx) == 1:
        idx = np.asscalar(idx)
    return idx


def get_subset_slices(l_obs_d, l_obs_o):
    if l_obs_d[0] > l_obs_o[0]:
        l1_d = 0
        l1_o = find_nearest_index(l_obs_o, l_obs_d[0])
    else:
        l1_d = find_nearest_index(l_obs_d, l_obs_o[0])
        l1_o = 0

    if l_obs_d[-1] < l_obs_o[-1]:
        l2_d = len(l_obs_d)
        l2_o = find_nearest_index(l_obs_o, l_obs_d[-1]) + 1
    else:
        l2_d = find_nearest_index(l_obs_d, l_obs_o[-1]) + 1
        l2_o = len(l_obs_o)
    return slice(l1_d, l2_d), slice(l1_o, l2_o)
