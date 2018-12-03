import numpy as np
from astropy.table import Table

def read_base_ascii(base_file, base_dir='./', out_fname=None, save_fits=False):

    # Third column may or may not be the upper age limit for CSP bases
    try:
        dt = np.dtype([('sspfile', '|S60'), ('age_base', '>f8'), ('Z_base', '>f8')
                     , ('age_base_upp', '>f8'), ('Mstars', '>f8'), ('YA_V', '>i4')
                     , ('aFe', '>f8') ])
        bdata = np.genfromtxt(base_file, dtype=dt, skip_header=1, usecols=(0,1,2,3,4,5,6))
    except Exception:
        dt = np.dtype([('sspfile', '|S60'), ('age_base', '>f8'), ('Z_base', '>f8'),
                       ('Mstars', '>f8'), ('YA_V', '>i4'), ('aFe', '>f8') ])
        bdata = np.genfromtxt(base_file, dtype=dt, skip_header=1, usecols=(0,1,2,4,5,6))

    base = Table(bdata)

    wl_ssp = np.genfromtxt(base['sspfile'][0]).transpose()[0]

    base['base_wl']    = np.full((len(base),len(wl_ssp)), wl_ssp)
    base['base_spec']  = np.empty((len(base),len(wl_ssp)))

    for i in range(len(base)):
        base['base_spec'][i] = np.genfromtxt(base['sspfile'][i]).transpose()[1]

    if save_fits==True:
        base.write(out_fname, overwrite=True)

    return base
