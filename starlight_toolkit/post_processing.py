import numpy as np

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


def calc_atflux(age_base, popx):
    return np.sum(np.log10(age_base) * popx)/popx.sum()


def calc_atmass(age_base, popmu):
    return np.sum(np.log10(age_base) * popmu)/popmu.sum()

   
def calc_aZflux(Z_base, popx, Z_sun): 
    return (popx * np.log10(Z_base/Z_sun)).sum()/ popx.sum()


def calc_aZmass(Z_base, popmu, Z_sun): 
    return (popmu * np.log10(Z_base/Z_sun)).sum()/popmu.sum()


#def calc_aZflux(Z_base, popx, Z_sun): 
    #return np.log10((popx * Z_base/Z_sun).sum()/popx.sum())


#def calc_aZmass(Z_base, popmu, Z_sun): 
    #return np.log10((popmu * Z_base/Z_sun).sum()/popmu.sum())



