import os
import gzip
import bz2
import numpy as np
from astropy.table import Table



def read_output_table(filename):
    '''
    Read STARLIGHT output file - version PANcMExStarlightChains_v03.for ...

    ElCid@Sanchica - 13/Feb/2012

    --------------------------------------------------------------------
             Some notes on the structure of STARLIGHT output files
    --------------------------------------------------------------------

     Considering the 1st line to be number 1 (ATT: **subtract** 1 in python!!): 
         * The pop-vector block starts at line 64 and has N_base entries from n1 to n2:
            n1 = 64
            n2 = n1 + N_base - 1

        * The Average & Chains light fractions block starts at line 64 + N_base + 5 and has N_par entries
            n1 = 64 + N_base + 5 
            n2 = n1 + N_par - 1, where N_par = N_base + 2 + N_exAV
        * then comes the chain-LAx-pop-vectors block, whose N_base elements go from
            n1 = 64 + N_base + 5 + N_par + 2
            n2 = n1 + N_base - 1
        * and the chain-mu_cor-pop-vectors block, whose N_base elements go from
             n1 = 64 + N_base + 5 + N_par + 2 + N_base + 2
            n2 = n1 + N_base - 1
        * The chain chi2's and masses are in lines    
             64 + N_base + 5 + N_par + 2 + N_base + 2 + N_base + 2, and
             64 + N_base + 5 + N_par + 2 + N_base + 2 + N_base + 3, respectively, 
             followed by 2 lines with v_0_before_EX0s and v_d_before_EX0s.

        * The specral block starts with new line containing Nl_obs, index_Best_SSP, 
         and i_SaveBestSingleCompFit at 
              64 + N_base + 5 + N_par + 2 + N_base + 2 + N_base + 10
        * The l_obs , f_obs , f_syn , f_wei , Best_f_SSP info has Nl_obs entries, running from 
             n1 = 64 + N_base + 5 + N_par + 2 + N_base + 2 + N_base + 11
            n2 = n1 + Nl_obs -1 

        * The FIRc/QHRc/PHOc-related ouput is given at the end of the file, and 
        I still have to implement their reading here!
    --------------------------------------------------------------------
'''
    if not os.path.exists(filename):
        raise Exception('File not found: %s' % filename)

    if filename.endswith('.gz'):
        f = gzip.GzipFile(filename)
    elif filename.endswith('.bz2'):
        f = bz2.BZ2File(filename)
    else:
        f = open(filename)
    data = f.read().splitlines()
    f.close()

    fileVersion = data[1].split()[3]
    # FIXME: How to handle starlight file versions?
#    if fileVersion != _getStarlightVersion():
#        raise Exception('Incorrect starlight version. Expected: %s, found: %' (_getStarlightVersion(), fileVersion))
    keywords = {}
    tables = {}
    tables['keywords'] = keywords
    keywords['file_version'] = fileVersion
    keywords['arq_synt'] = os.path.basename(filename)

    #--------------------------------------------------------------------
    # Read "header": general info plus pop-vectors
    #--------------------------------------------------------------------

    # Some input info
    keywords['arq_spec'] = data[5].split()[0]
    keywords['obs_dir'] = data[5].split()[2]
    keywords['arq_base'] = data[6].split()[0]
    keywords['arq_masks'] = data[7].split()[0]
    keywords['arq_config'] = data[8].split()[0]
    keywords['N_base'] = int(data[9].split()[0])
    keywords['N_exAV_components'] = int(data[10].split()[0])
    keywords['N_exAV'] = int(data[10].split()[1])
    keywords['IsFIRcOn'] = int(data[10].split()[2])
    keywords['IsQHRcOn'] = int(data[10].split()[3])
    keywords['IsPHOcOn'] = int(data[10].split()[4])
    keywords['iFitPowerLaw'] = int(data[11].split()[0])
    keywords['alpha_PowerLaw'] = float(data[12].split()[0])
    keywords['red_law_option'] = data[13].split()[0]
    keywords['q_norm'] = float(data[14].split()[0])
    keywords['flux_unit'] = float(data[14].split()[1])

    # (Re)Sampling Parameters
    keywords['l_ini'] = float(data[17].split()[0])
    keywords['l_fin'] = float(data[18].split()[0])
    keywords['dl'] = float(data[19].split()[0])
    keywords['dl_cushion'] = float(data[19].split()[1])

    # Normalization info
    keywords['l_norm'] = float(data[22].split()[0])
    keywords['llow_norm'] = float(data[23].split()[0])
    keywords['lupp_norm'] = float(data[24].split()[0])
    keywords['fobs_norm'] = float(data[25].split()[0])
    keywords['Lobs_norm'] = float(data[25].split()[1])
    keywords['LumDistInMpc'] = float(data[25].split()[2])

    # S/N
    keywords['llow_SN'] = float(data[28].split()[0])
    keywords['lupp_SN'] = float(data[29].split()[0])
    keywords['SN_snwin'] = float(data[30].split()[0])
    keywords['SN_normwin'] = float(data[31].split()[0])
    keywords['SNerr_snwin'] = float(data[32].split()[0])
    keywords['SNerr_normwin'] = float(data[33].split()[0])
    keywords['fscale_chi2'] = float(data[34].split()[0])
    keywords['IsOptimize_fn_OPT'] = float(data[34].split()[1])

    # etc...
    keywords['idum_orig'] = int(data[37].split()[0])
    keywords['NOl_eff'] = int(data[38].split()[0])
    keywords['Nl_eff'] = int(data[39].split()[0])
    keywords['Ntot_clipped'] = int(data[40].split()[0])
    keywords['clip_method'] = data[40].split()[1]
    keywords['Nglobal_steps'] = int(data[41].split()[0])
    keywords['N_chains'] = int(data[42].split()[0])
    keywords['NEX0s_base'] = int(data[43].split()[0])
    keywords['iCLIPBUG_flag'] = int(data[44].split()[0])
    keywords['i_RC_CRASH_FLAG'] = int(data[44].split()[1])
    keywords['IsBurInOver_BurnIn'] = int(data[44].split()[2])
    keywords['n_censored_weights'] = int(data[44].split()[3])
    keywords['wei_nsig_threshold'] = float(data[44].split()[4])
    keywords['wei_limit'] = float(data[44].split()[5])
    keywords['idt_all'] = int(data[45].split()[0])
    keywords['wdt_TotTime'] = float(data[45].split()[1])
    keywords['wdt_UsrTime'] = float(data[45].split()[2])
    keywords['wdt_SysTime'] = float(data[45].split()[3])

    ## Synthesis Results - Best model ##
    keywords['chi2'] = float(data[49].split()[0])
    keywords['adev'] = float(data[50].split()[0])
    keywords['chi2_TOT'] = float(data[51].split()[0])
    keywords['chi2_Opt'] = float(data[51].split()[1])
    keywords['chi2_FIR'] = float(data[51].split()[2])
    keywords['chi2_QHR'] = float(data[51].split()[3])
    keywords['chi2_PHO'] = float(data[51].split()[4])

    keywords['sum_x'] = float(data[52].split()[0])
    keywords['Flux_tot'] = float(data[53].split()[0])
    keywords['Mini_tot'] = float(data[54].split()[0])
    keywords['Mcor_tot'] = float(data[55].split()[0])

    keywords['v_0'] = float(data[57].split()[0])
    keywords['v_d'] = float(data[58].split()[0])
    keywords['A_V'] = float(data[59].split()[0])

    # Reset populations lists
    popx = []    # column 2
    popmu_ini = []    # column 3
    popmu_cor = []    # column 4
    popage_base = []    # column 5
    popZ_base = []    # column 6
    popfbase_norm = []    # column 7
    popexAV_flag = []    # column 8
    popMstars = []    # column 9
    aFe = []    # column 11
    SSP_chi2r = []    # column 12
    SSP_adev = []    # column 13
    SSP_AV = []    # column 14
    SSP_x = []    # column 15
    popAV_tot = []    # column 16
    popLAx = []    # column 17

    # j     x_j(%)      Mini_j(%)     Mcor_j(%)     age_j(yr)     Z_j      (L/M)_j   exAV?  Mstars   component_j        new/Fe...    |  SSP_chi2r SSP_adev(%)   SSP_AV   SSP_x(%)    |  AV_tot   <LAx>_j(%)
    # Reads all these things (as lists) from lines _n1 to _n2
    _n1 = 63
    _n2 = _n1 + keywords['N_base']
    for i in range(_n1, _n2):
        popx.append(float(data[i].split()[1]))
        popmu_ini.append(float(data[i].split()[2]))
        popmu_cor.append(float(data[i].split()[3]))
        popage_base.append(float(data[i].split()[4]))
        popZ_base.append(float(data[i].split()[5]))
        popfbase_norm.append(float(data[i].split()[6]))
        popexAV_flag.append(float(data[i].split()[7]))
        popMstars.append(float(data[i].split()[8]))
        aFe.append(float(data[i].split()[10]))
        SSP_chi2r.append(float(data[i].split()[11]))
        SSP_adev.append(float(data[i].split()[12]))
        SSP_AV.append(float(data[i].split()[13]))
        SSP_x.append(float(data[i].split()[14]))
        popAV_tot.append(float(data[i].split()[15]))
        popLAx.append(float(data[i].split()[16]))

    # ATT: In my SM macros I fix PL-age & Z, and renormalize light fractions to 100%.
    #      This could be done like this:
    #          aux = np.array(popx)   ; popx   = np.array(100. * aux / aux.sum())
    #          aux = np.array(popLAx) ; popLAx = np.array(100. * aux / aux.sum())
    #      But I'm **NOT** doing this here since all this function does is to replicate the
    #       STARLIGHT output to new python dictionary! Furthermore, renomrmalization may
    #        be inconvenient to create luminosity x age 3D-CALIFA-images...

    # ATT2: Ignoring Power-law fixes.
    # OBS: PL have age = 0 in the Starlight output file:(
    #      Here I change it so that age_PL = 5e5 yr... & Z_PL = solar
    #      This is all obsolete anyway. The built-in PL is not used anymore.
    # if (int(StarlightOut['iFitPowerLaw']) > 0):
    #    print '@@> [Warning!] ...Fixing PL age & Z ...Why are u using PLs??'
    #    popage_base][StarlightOut['N_base'] - 1] = 5e5
    #    popZ_base[StarlightOut['N_base'] - 1]    = 0.02
    cols = [popx,
            popmu_ini,
            popmu_cor,
            popage_base,
            popZ_base,
            popfbase_norm,
            popexAV_flag,
            popMstars,
            aFe,
            SSP_chi2r,
            SSP_adev,
            SSP_AV,
            SSP_x,
            popAV_tot]

    names = ['popx',
             'popmu_ini',
             'popmu_cor',
             'popage_base',
             'popZ_base',
             'popfbase_norm',
             'popexAV_flag',
             'popMstars',
             'aFe',
             'p_chi2r',
             'p_adev',
             'p_AV',
             'p_x',
             'popAV_tot']
    tables['population'] = Table(cols, names=names)

    #--------------------------------------------------------------------
    # Read chain-related info (in arrays!)
    # ATT: The average solution is counted as an extra chain entry - the 1st one (Chain_*[0])!
    #--------------------------------------------------------------------
    ## Synthesis Results - Average & Chains ##

    # j    par_j: min, <> & last-chain-values for 0 ... N_chains chains (ave is 0'th chain!)
    # Reading chain par-vector (light fractions at l_norm, + extinctions and etc)
    # Notice that Chain_Par contains AV (maybe more than 1!) and fn, as well
    # as x!
    N_par = keywords['N_base'] + 1 + keywords['N_exAV'] + 1
    _n1 = 63 + keywords['N_base'] + 6 - 1
    _n2 = _n1 + N_par - 1 + 1
    Best_Par = []
    Ave_Par = []
    Chain_Par = []
    for unused, i in enumerate(range(_n1, _n2)):
        Best_Par.append(np.float(data[i].split()[1]))
        Ave_Par.append(np.float(data[i].split()[2]))
        x_ = [np.float(x) for x in data[i].split()[3:3 + keywords['N_chains']]]
        Chain_Par.append(x_)

    # j Lambda-Averaged pop-vectors <LAx_*>_j: min, <> & last-chain-values for 0 ... N_chains chains (ave is 0'th chain!)
    # Reading chain LAx pop-vectors
    _n1 = 63 + keywords['N_base'] + 6 - 1 + N_par - 1 + 1 + 2
    _n2 = _n1 + keywords['N_base']
    Best_LAx = []
    Ave_LAx = []
    Chain_LAx = []
    for unused, i in enumerate(range(_n1, _n2)):
        Best_LAx.append(np.float(data[i].split()[1]))
        Ave_LAx.append(np.float(data[i].split()[2]))
        x_ = [np.float(x) for x in data[i].split()[3:3 + keywords['N_chains']]]
        Chain_LAx.append(x_)

    # j   Mcor_j: min, <> & last-chain-values for 0 ... N_chains chains (ave is 0'th chain!)
    # Reading chain mu_cor pop-vectors
    _n1 = 63 + keywords['N_base'] + 6 - 1 + \
        N_par - 1 + 1 + 2 + keywords['N_base'] + 2
    _n2 = _n1 + keywords['N_base']
    Best_mu_cor = []
    Ave_mu_cor = []
    Chain_mu_cor = []
    for unused, i in enumerate(range(_n1, _n2)):
        Best_mu_cor.append(np.float(data[i].split()[1]))
        Ave_mu_cor.append(np.float(data[i].split()[2]))
        x_ = [np.float(x) for x in data[i].split()[3:3 + keywords['N_chains']]]
        Chain_mu_cor.append(x_)

    # chi2/Nl_eff & Mass for min, <> & i_chain = 0 ...  N_chains chains (ave is 0'th chain!)
    # Read Chain chi2/Nl_eff's , as well as kinematics before_EX0s
    i = 63 + keywords['N_base'] + 6 - 1 + N_par - 1 + 1 + \
        2 + keywords['N_base'] + 2 + keywords['N_base'] + 2
    keywords['best_chi2'] = np.float(data[i].split()[1])
    keywords['ave_chi2'] = np.float(data[i].split()[2])
    keywords['cha_chi2'] = [
        np.float(x) for x in data[i].split()[3:3 + keywords['N_chains']]]

    keywords['best_Mcor'] = np.float(data[i + 1].split()[1])
    keywords['ave_Mcor'] = np.float(data[i + 1].split()[2])
    keywords['cha_Mcor'] = [
        np.float(x) for x in data[i + 1].split()[3:3 + keywords['N_chains']]]

    keywords['v_0_before_EX0s'] = float(data[i + 2].split()[0])
    keywords['v_d_before_EX0s'] = float(data[i + 3].split()[0])

    # Store chains in tables.
    cols = [Best_Par, Ave_Par, Chain_Par]
    chains_names = ['best', 'average', 'chains']
    tables['chains_par'] = Table(cols, names=chains_names)

    cols = [Best_LAx, Ave_LAx, Chain_LAx]
    tables['chains_LAx'] = Table(cols, names=chains_names)

    cols = [Best_mu_cor, Ave_mu_cor, Chain_mu_cor]
    tables['chains_mu_cor'] = Table(cols, names=chains_names)
    #--------------------------------------------------------------------

    N_par = keywords['N_base'] + 1 + keywords['N_exAV'] + 1

    # Synthetic spectrum (Best Model) ##l_obs f_obs f_syn wei Best_f_SSP

    i = 63 + keywords['N_base'] + 6 - 1 + N_par - 1 + 1 + \
        2 + keywords['N_base'] + 2 + keywords['N_base'] + 2 + 8
    keywords['Nl_obs'] = int(data[i].split()[0])
    keywords['index_Best_SSP'] = int(data[i].split()[1])
    keywords['i_SaveBestSingleCompFit'] = int(data[i].split()[2])

    # Reset & read spectral arrays (later turned into numpy.arrays)
    l_obs = []
    f_obs = []
    f_syn = []
    f_wei = []
    Best_f_SSP = []

    # Read spectra. Notice that new 5th column (with Best_f_SSP) is only read
    # & returned if it actually exists!
    _n1 = i + 1
    _n2 = _n1 + keywords['Nl_obs']
    for i in range(_n1, _n2):
        l_obs.append(float(data[i].split()[0]))
        f_obs.append(float(data[i].split()[1]))
        f_syn.append(float(data[i].split()[2]))
        f_wei.append(float(data[i].split()[3]))
        if (keywords['i_SaveBestSingleCompFit'] == 1):
            Best_f_SSP.append(float(data[i].split()[4]))

    cols = [l_obs, f_obs, f_syn, f_wei]
    names = ['l_obs', 'f_obs', 'f_syn', 'f_wei']
    if (keywords['i_SaveBestSingleCompFit'] == 1):
        cols.append(Best_f_SSP)
        names.append('Best_f_SSP')

    tables['spectra'] = Table(cols, names=names)

    #--------------------------------------------------------------------
    # Reading FIR-related output
    #--------------------------------------------------------------------
    if (keywords['IsFIRcOn'] != 0):
        N_par = keywords['N_base'] + 1 + keywords['N_exAV'] + 1

        # Skip spectra
        i = 63 + keywords['N_base'] + 6 - 1 + N_par - 1 + 1 + 2 + \
            keywords['N_base'] + 2 + keywords['N_base'] + 2 + 8
        keywords['Nl_obs'] = int(data[i].split()[0])
        _n1 = i + 1
        _n2 = _n1 + keywords['Nl_obs']
        _n3 = _n2 + 8

        keywords['FIR_arq_ETCinfo'] = (data[_n3].split()[0])
        _n3 += 1
        keywords['FIR_LumDistInMpc'] = float(data[_n3].split()[0])
        _n3 += 1
        keywords['FIR_logLFIR_TOTInLsun'] = float(data[_n3].split()[0])
        _n3 += 1
        keywords['FIR_LFIRFrac2Model'] = float(data[_n3].split()[0])
        _n3 += 1
        keywords['FIR_logLFIR_obsInLsun'] = float(data[_n3].split()[0])
        _n3 += 1
        keywords['FIR_ErrlogLFIRInDex'] = float(data[_n3].split()[0])
        _n3 += 1
        keywords['FIR_RangelogLFIRInDex'] = float(data[_n3].split()[0])
        _n3 += 1
        keywords['FIRChi2ScaleFactor'] = float(data[_n3].split()[0])
        _n3 += 1
        keywords['FIR_logLFIR_lowInLsun'] = float(data[_n3].split()[0])
        keywords['FIR_logLFIR_uppInLsun'] = float(data[_n3].split()[1])
        _n3 += 1
        keywords['FIRbeta_D'] = float(data[_n3].split()[0])
        keywords['FIRbeta_I'] = float(data[_n3].split()[1])
        _n3 += 1
        keywords['log_LFIR/LOpt_rough'] = float(data[_n3].split()[0])
        _n3 += 1

        _n3 += 2

        keywords['FIRModlogLFIRInLsun'] = float(data[_n3].split()[0])
        keywords['FIRModObsRatio'] = float(data[_n3].split()[1])
        _n3 += 1
        keywords['FIRModlogLBOLInLsun'] = float(data[_n3].split()[0])
        keywords['FIR_BOL_Ratio'] = float(data[_n3].split()[1])
        _n3 += 1
        keywords['chi2_FIR'] = float(data[_n3].split()[0])
        keywords['chi2_Opt'] = float(data[_n3].split()[1])
        _n3 += 1
        keywords['chi2_FIR/TOT_Perc'] = float(data[_n3].split()[0])
        keywords['chi2_Opt/TOT_Perc'] = float(data[_n3].split()[1])
        _n3 += 1

        _n3 += 2

        # Reset & read FIR-related SSP arrays
        x_FIR = []
        x_BOL = []
        BolCor = []
        FracLion = []
        Lbol_M = []
        Rmat = []
        R_opt = []
        R_Lya = []
        R_LCE = []

        # Read FIR-related SSP arrays
        _n1 = _n3
        _n2 = _n1 + keywords['N_base']
        for i in range(_n1, _n2):
            x_FIR.append(float(data[i].split()[6]))
            x_BOL.append(float(data[i].split()[7]))
            BolCor.append(float(data[i].split()[8]))
            FracLion.append(float(data[i].split()[9]))
            Lbol_M.append(float(data[i].split()[10]))
            Rmat.append(float(data[i].split()[11]))
            R_opt.append(float(data[i].split()[12]))
            R_Lya.append(float(data[i].split()[13]))
            R_LCE.append(float(data[i].split()[14]))

        cols = [x_FIR,
                x_BOL,
                BolCor,
                FracLion,
                Lbol_M,
                Rmat,
                R_opt,
                R_Lya,
                R_LCE]

        names = ['x_FIR',
                 'x_BOL',
                 'BolCor',
                 'FracLion',
                 'Lbol_M',
                 'Rmat',
                 'R_opt',
                 'R_Lya',
                 'R_LCE']

        tables['FIR'] = Table(cols, names=names)

    #--------------------------------------------------------------------
    # Reading QHR-related output
    #--------------------------------------------------------------------
    if (keywords['IsQHRcOn'] != 0):
        N_par = keywords['N_base'] + 1 + keywords['N_exAV'] + 1

        # Skip spectra
        i = 63 + keywords['N_base'] + 6 - 1 + N_par - 1 + 1 + 2 + \
            keywords['N_base'] + 2 + keywords['N_base'] + 2 + 8
        keywords['Nl_obs'] = int(data[i].split()[0])
        _n1 = i + 1
        _n2 = _n1 + keywords['Nl_obs']
        _n3 = _n2 + 8

        # Skip FIR
        if (keywords['IsFIRcOn'] != 0):
            _n3 = _n3 + 27 + keywords['N_base']

        keywords['QHRbeta_I'] = float(data[_n3].split()[0])
        _n3 += 1
        keywords['QHR_arq_ETCinfo'] = (data[_n3].split()[0])
        _n3 += 1
        keywords['QHR_LumDistInMpc'] = float(data[_n3].split()[0])
        _n3 += 1
        keywords['QHR_GlobalChi2ScaleFactor'] = float(
            data[_n3].split()[0])
        _n3 += 1
        keywords['NQHR_Ys'] = int(data[_n3].split()[0])
        _n3 += 1

        _n3 += 2

        # Reset & read QHR observed
        QHR_lambda = []
        QHR_frecomb = []
        QHR_logY_TOT = []
        QHR_YFrac2Model = []
        QHR_ErrlogY = []
        QHR_RangelogY = []
        QHR_Chi2ScaleFactor = []
        QHR_logY_obs = []
        QHR_logY_low = []
        QHR_logY_upp = []

        # Read QHR observed
        _n1 = _n3
        _n2 = _n1 + keywords['NQHR_Ys']
        for i in range(_n1, _n2):
            QHR_lambda.append(float(data[i].split()[1]))
            QHR_frecomb.append(float(data[i].split()[2]))
            QHR_logY_TOT.append(float(data[i].split()[3]))
            QHR_YFrac2Model.append(float(data[i].split()[4]))
            QHR_ErrlogY.append(float(data[i].split()[5]))
            QHR_RangelogY.append(float(data[i].split()[6]))
            QHR_Chi2ScaleFactor.append(float(data[i].split()[7]))
            QHR_logY_obs.append(float(data[i].split()[8]))
            QHR_logY_low.append(float(data[i].split()[9]))
            QHR_logY_upp.append(float(data[i].split()[10]))

        cols = [QHR_lambda,
                QHR_frecomb,
                QHR_logY_TOT,
                QHR_YFrac2Model,
                QHR_ErrlogY,
                QHR_RangelogY,
                QHR_Chi2ScaleFactor,
                QHR_logY_obs,
                QHR_logY_low,
                QHR_logY_upp]

        names = ['lambda',
                 'frecomb',
                 'logY_TOT',
                 'YFrac2Model',
                 'ErrlogY',
                 'RangelogY',
                 'Chi2ScaleFactor',
                 'logY_obs',
                 'logY_low',
                 'logY_upp']

        tables['QHR_Obs'] = Table(cols, names=names)

        # Read Emission Line Ratio-related things
        _n3 = _n2 + 2

        keywords['IsELROn'] = int(data[_n3].split()[0])
        keywords['ELR_lambda_A'] = float(data[_n3].split()[1])
        keywords['ELR_lambda_B'] = float(data[_n3].split()[2])
        keywords['ELR_ind_A'] = int(data[_n3].split()[3])
        keywords['ELR_ind_B'] = int(data[_n3].split()[4])
        keywords['ELR_logRint'] = float(data[_n3].split()[5])
        keywords['ELR_AV_neb'] = float(data[_n3].split()[6])
        keywords['ELR_errAV_neb'] = float(data[_n3].split()[7])
        _n3 += 1

        if (keywords['IsELROn'] == 1):

            keywords['ELR_Err_logR'] = float(data[_n3].split()[0])
            keywords['ELR_RangelogR'] = float(data[_n3].split()[1])
            keywords['ELR_logR_low'] = float(data[_n3].split()[2])
            keywords['ELR_logR_upp'] = float(data[_n3].split()[3])
            keywords['ELR_Chi2ScaleFactor'] = float(data[_n3].split()[4])
            _n3 += 1

            # WARNING!!! Starlight output labels logRobs & ModlogR in the wrong
            # order!
            keywords['ELR_ModlogR'] = float(data[_n3].split()[0])
            keywords['ELR_logRobs'] = float(data[_n3].split()[1])
            keywords['ELR_chi2_ELR'] = float(data[_n3].split()[2])
            keywords[
                'ELR_chi2_ELR/chi2_QHR'] = float(data[_n3].split()[3])
            keywords[
                'ELR_chi2_ELR/chi2_TOT'] = float(data[_n3].split()[4])
            _n3 += 1
        else:
            _n3 += 2

        _n3 += 2

        keywords['log_QH0_PhotPerSec'] = float(data[_n3].split()[0])
        keywords['log_QHeff_PhotPerSec'] = float(data[_n3].split()[1])
        _n3 += 1
        keywords['chi2_QHR/TOT_Perc'] = float(data[_n3].split()[0])
        keywords['chi2_Opt/TOT_Perc'] = float(data[_n3].split()[1])
        _n3 += 1
        keywords['chi2_QHR'] = float(data[_n3].split()[0])
        keywords['chi2_Opt'] = float(data[_n3].split()[1])
        keywords['chi2_TOT'] = float(data[_n3].split()[2])
        _n3 += 1

        _n3 += 2

        # Reset & read QHR model
        QHR_lambda = []
        QHR_q_lambda = []
        QHR_logY_obs = []
        QHR_ModlogY = []
        QHR_chi2_Y = []
        QHR_chi2Y_chi2Opt = []
        QHR_chi2Y_chi2TOT = []

        # Read QHR model
        _n1 = _n3
        _n2 = _n1 + keywords['NQHR_Ys']
        for i in range(_n1, _n2):
            QHR_lambda.append(float(data[i].split()[1]))
            QHR_q_lambda.append(float(data[i].split()[2]))
            QHR_logY_obs.append(float(data[i].split()[3]))
            QHR_ModlogY.append(float(data[i].split()[4]))
            QHR_chi2_Y.append(float(data[i].split()[5]))
            QHR_chi2Y_chi2Opt.append(float(data[i].split()[6]))
            QHR_chi2Y_chi2TOT.append(float(data[i].split()[7]))

        cols = [QHR_lambda,
                QHR_q_lambda,
                QHR_logY_obs,
                QHR_ModlogY,
                QHR_chi2_Y,
                QHR_chi2Y_chi2Opt,
                QHR_chi2Y_chi2TOT]

        names = ['lambda',
                 'q_lambda',
                 'logY_obs',
                 'ModlogY',
                 'chi2_Y',
                 'chi2Y_chi2Opt',
                 'chi2Y_chi2TOT']

        tables['QHR_Mod'] = Table(cols, names=names)

        _n3 = _n2 + 2

        # Reset & read FIR-related SSP arrays
        qH__40 = []
        QH2Lnorm__40 = []
        QH0_Perc = []
        QHeff_Perc = []
        Y_Perc = dict()
        for il in range(0, keywords['NQHR_Ys']):
            Y_Perc[il] = []

        # Read FIR-related SSP arrays
        _n1 = _n3
        _n2 = _n1 + keywords['N_base']
        for i in range(_n1, _n2):
            qH__40.append(float(data[i].split()[6]))
            QH2Lnorm__40.append(float(data[i].split()[7]))
            QH0_Perc.append(float(data[i].split()[8]))
            QHeff_Perc.append(float(data[i].split()[9]))
            for il in range(0, keywords['NQHR_Ys']):
                Y_Perc[il].append(float(data[i].split()[10 + il]))

        cols = [qH__40, QH2Lnorm__40, QH0_Perc, QHeff_Perc]
        names = ['qH__40', 'QH2Lnorm__40', 'QH0_Perc', 'QHeff_Perc']
        for il in range(0, keywords['NQHR_Ys']):
            cols.append(Y_Perc[il])
            names.append('Y_Perc_Line' + str(il))
        tables['QHR'] = Table(cols, names=names)

    #--------------------------------------------------------------------
    # Reading PHO-related output
    #--------------------------------------------------------------------
    if (keywords['IsPHOcOn'] != 0):
        N_par = keywords['N_base'] + 1 + keywords['N_exAV'] + 1

        # Skip spectra
        i = 63 + keywords['N_base'] + 6 - 1 + N_par - 1 + 1 + 2 + \
            keywords['N_base'] + 2 + keywords['N_base'] + 2 + 8
        keywords['Nl_obs'] = int(data[i].split()[0])
        _n1 = i + 1
        _n2 = _n1 + keywords['Nl_obs']
        _n3 = _n2 + 8

        # Skip FIR
        if (keywords['IsFIRcOn'] != 0):
            _n3 = _n3 + 27 + keywords['N_base']

        # Skip QHR
        if (keywords['IsQHRcOn'] != 0):
            _n3 = _n3 + 29 + 2 * \
                keywords['NQHR_Ys'] + keywords['N_base']

        keywords['PHO_arq_ETCinfo'] = (data[_n3].split()[0])
        _n3 += 1
        keywords['PHO_LumDistInMpc'] = float(data[_n3].split()[0])
        _n3 += 1
        keywords['PHO_Redshift'] = float(data[_n3].split()[0])
        _n3 += 1
        keywords['PHO_GlobalChi2ScaleFactor'] = float(
            data[_n3].split()[0])
        _n3 += 1
        keywords['NPHO_Ys'] = int(data[_n3].split()[0])
        _n3 += 1

        _n3 += 2

        # Reset & read PHO observed
        PHO_name = []
        PHO_logY_TOT = []
        PHO_YFrac2Model = []
        PHO_ErrlogY = []
        PHO_RangelogY = []
        PHO_Chi2ScaleFactor = []
        PHO_logY_obs = []
        PHO_logY_low = []
        PHO_logY_upp = []

        # Read PHO observed
        _n1 = _n3
        _n2 = _n1 + keywords['NPHO_Ys']
        for i in range(_n1, _n2):
            PHO_name.append(data[i].split()[0])
            PHO_logY_TOT.append(float(data[i].split()[2]))
            PHO_YFrac2Model.append(float(data[i].split()[3]))
            PHO_ErrlogY.append(float(data[i].split()[4]))
            PHO_RangelogY.append(float(data[i].split()[5]))
            PHO_Chi2ScaleFactor.append(float(data[i].split()[6]))
            PHO_logY_obs.append(float(data[i].split()[7]))
            PHO_logY_low.append(float(data[i].split()[8]))
            PHO_logY_upp.append(float(data[i].split()[9]))

        cols = [PHO_name,
                PHO_logY_TOT,
                PHO_YFrac2Model,
                PHO_ErrlogY,
                PHO_RangelogY,
                PHO_Chi2ScaleFactor,
                PHO_logY_obs,
                PHO_logY_low,
                PHO_logY_upp]

        names = ['name',
                 'logY_TOT',
                 'YFrac2Model',
                 'ErrlogY',
                 'RangelogY',
                 'Chi2ScaleFactor',
                 'logY_obs',
                 'logY_low',
                 'logY_upp']

        tables['PHO_Obs'] = Table(cols, names=names)

        _n3 = _n2 + 2

        keywords['chi2_PHO/TOT_Perc'] = float(data[_n3].split()[0])
        keywords['chi2_Opt/TOT_Perc'] = float(data[_n3].split()[1])
        _n3 += 1
        keywords['chi2_PHO'] = float(data[_n3].split()[0])
        keywords['chi2_Opt'] = float(data[_n3].split()[1])
        keywords['chi2_TOT'] = float(data[_n3].split()[2])
        _n3 += 1

        _n3 += 2

#  name/code  MeanLamb StdDevLamb  q_MeanLamb   logY_obs      ModlogY       chi2_Y   chi2_Y/chi2_Opt  chi2_Y/chi2_TOT
        # Reset & read PHO model
        PHO_name = []
        PHO_MeanLamb = []
        PHO_StdDevLamb = []
        PHO_q_MeanLamb = []
        PHO_logY_obs = []
        PHO_ModlogY = []
        PHO_chi2_Y = []
        PHO_chi2Y_chi2Opt = []
        PHO_chi2Y_chi2TOT = []

        # Read PHO model
        _n1 = _n3
        _n2 = _n1 + keywords['NPHO_Ys']
        for i in range(_n1, _n2):
            PHO_name.append(data[i].split()[1])
            PHO_MeanLamb.append(float(data[i].split()[2]))
            PHO_StdDevLamb.append(float(data[i].split()[3]))
            PHO_q_MeanLamb.append(float(data[i].split()[4]))
            PHO_logY_obs.append(float(data[i].split()[5]))
            PHO_ModlogY.append(float(data[i].split()[6]))
            PHO_chi2_Y.append(float(data[i].split()[7]))
            PHO_chi2Y_chi2Opt.append(float(data[i].split()[8]))
            PHO_chi2Y_chi2TOT.append(float(data[i].split()[9]))

        cols = [PHO_name,
                PHO_MeanLamb,
                PHO_StdDevLamb,
                PHO_q_MeanLamb,
                PHO_logY_obs,
                PHO_ModlogY,
                PHO_chi2_Y,
                PHO_chi2Y_chi2Opt,
                PHO_chi2Y_chi2TOT]
        names = ['name',
                 'MeanLamb',
                 'StdDevLamb',
                 'q_MeanLamb',
                 'logY_obs',
                 'ModlogY',
                 'chi2_Y',
                 'chi2Y_chi2Opt',
                 'chi2Y_chi2TOT']

        tables['PHO_Mod'] = Table(cols, names=names)

        _n3 = _n2 + 2

        # Reset & read PHO-related SSP arrays
        Y_Perc = []

        # Read PHO-related SSP arrays
        _n1 = _n3
        _n2 = _n1 + keywords['N_base']
        for i in range(_n1, _n2):
            Y_Perc.append([float(x) for x in data[i].split()[6:]])

        cols = [Y_Perc]
        names = [' Y_Perc']
        tables['PHO'] = Table(cols, names=names)

    return tables
