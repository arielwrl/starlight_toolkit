import numpy as np

#A function to calculate star-formation histories:
def calc_sfh(age_base, popmu):
    #A vector with unique ages:
    agevec = np.unique(age_base)
    #The SFH:
    sfh  = [np.sum(popmu[age_base==agevec[i]]) for i in range(len(agevec))]
    #The cumulative SFH, which is what we will plot:
    csfh = np.cumsum(sfh[::-1])
    return agevec, sfh, csfh[::-1]

#A function to calculate mean stellar ages:
def calc_atflux(age_base, popx):
    return np.sum(np.log10(age_base) * popx)/popx.sum()

def calc_atmass(age_base, popx):
    return np.sum(np.log10(age_base) * popmu)/popmu.sum()
 
