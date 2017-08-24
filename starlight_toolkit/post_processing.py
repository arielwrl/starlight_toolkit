import numpy as np

def calc_sfh(age_base, popmu):
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
   
def calc_meanZflux(Z_base, popx, Z_sun): 
    return (popx * np.log10(Z_base/0.02)).sum()/ popx.sum()

def calc_meanZmass(Z_base, popmu, Z_sun): 
    return (popmu * np.log10(Z_base/0.02)).sum()/popmu.sum()
