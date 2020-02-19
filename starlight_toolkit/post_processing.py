import numpy as np
from starlight_toolkit.synphot import resampler


def convert_x_lambda(popx, wl, wl_0, base_wl, base_f):
#FIXME: lacking implementation for exAV fits.
    
    #Resample and normalize:
    base_wl_res = np.arange(np.round(base_wl[0]), np.round(base_wl[-1]), 1)
    base_f_res  = np.array([resampler(base_wl, base_f[i], base_wl_res) for i in range(len(base_f))])
    
    base_f_norm = [base_f_res[i][base_wl_res==wl]/base_f_res[i][base_wl_res==wl_0]  for i in range(len(base_f))]
    base_f_norm = np.array([base_f_norm[i][0]  for i in range(len(base_f))])


    #Convert x to lambda:
    popx_wl =  np.array([popx[i]*base_f_norm/np.sum(popx[i]*base_f_norm) for i in range(len(popx))])
    
    return popx_wl

def calc_sfh(age_base, popmu, base_type='SSP'):
    if base_type=='SSP':
        #A vector with unique ages:
        agevec = np.unique(age_base)
        #The SFH:
        sfh  = [np.sum(popmu[age_base==agevec[i]]) for i in range(len(agevec))]
        sfh  /= popmu.sum() 
        sfh *= 100
        #The cumulative SFH:
        csfh = np.cumsum(sfh[::-1])
    return agevec, sfh, csfh[::-1]


def calc_sfh_x(age_base, popx):
    #A vector with unique ages:
    agevec = np.unique(age_base)
    #The SFH:
    sfh  = [np.sum(popx[age_base==agevec[i]]) for i in range(len(agevec))]
    sfh  /= popx.sum() 
    sfh *= 100
    #The cumulative SFH:
    csfh = np.cumsum(sfh[::-1])
    return agevec, sfh, csfh[::-1]


def calc_Zfh(age_base, Z_base, popmu, Z_sun):
    #A vector with unique ages:
    agevec = np.unique(age_base)
    #The SFH:
    Zfh  = [(popmu[age_base==agevec[i]] * np.log10(Z_base[age_base==agevec[i]]/Z_sun)).sum()/popmu[age_base==agevec[i]].sum() for i in range(len(agevec))]
    for i in range(len(Zfh)):
        if np.isnan(Zfh[i]):
            Zfh[i] = 0
    #The cumulative ZFH:
    cZfh = np.cumsum(Zfh[::-1])
    return agevec, Zfh, cZfh[::-1]



def calc_QHRpop_x(age_base, popQHR):
    #A vector with unique ages:
    agevec = np.unique(age_base)
    #The SFH:
    QHRvec  = [np.sum(popQHR[age_base==agevec[i]]) for i in range(len(agevec))]
    QHRvec  /= popQHR.sum() 
    QHRvec *= 100
    #The cumulative SFH:
    cQHRvec = np.cumsum(QHRvec[::-1])
    return cQHRvec[::-1]


def calc_atflux(age_base, age_base_upp, popx):
    if age_base_upp is not None:
        log_t1 = np.log10(age_base)
        log_t2 = np.log10(age_base_upp)
        log_t  = (log_t1 + log_t2) / 2.0
    else:
        log_t = np.log10(age_base)
    return np.sum(log_t * popx) / popx.sum()


def calc_atmass(age_base, age_base_upp, popmu):
    if age_base_upp is not None:
        log_t1 = np.log10(age_base)
        log_t2 = np.log10(age_base_upp)
        log_t  = (log_t1 + log_t2) / 2.0
    else:
        log_t  = np.log10(age_base)        
    return np.sum(log_t * popmu) / popmu.sum()

   
def calc_aZflux(Z_base, popx, Z_sun): 
    return (popx * np.log10(Z_base/Z_sun)).sum()/ popx.sum()


def calc_aZmass(Z_base, popmu, Z_sun): 
    return (popmu * np.log10(Z_base/Z_sun)).sum()/popmu.sum()




