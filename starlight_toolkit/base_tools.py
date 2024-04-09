from os import path

import numpy as np
from astropy.table import Table

def read_base_ascii(base_file, base_dir='./', out_fname=None, save_fits=False):
    """
    Reads data from a Starlight base file into an astropy table, with the option of saving
    the table to a fits file.

    Parameters:
        base_file (str): The name of the base ASCII file.
        base_dir (str, optional): The directory path where the base file is located. Default is './'.
        out_fname (str, optional): The filename for saving the output FITS file. Default is None.
        save_fits (bool, optional): Whether to save the output as a FITS file. Default is False.

    Returns:
        astropy.table.Table: A table containing the data read from the base ASCII file.

    Raises:
        FileNotFoundError: If the specified base file does not exist.
    """

    bfile = path.join(base_dir, base_file)
    
    try:
        dt = np.dtype([('sspfile', '|S60'), ('age_base', '>f8'), ('Z_base', '>f8')
                     , ('age_base_upp', '>f8'), ('Mstars', '>f8'), ('YA_V', '>i4')
                     , ('aFe', '>f8') ])
        bdata = np.genfromtxt(bfile, dtype=dt, skip_header=1, usecols=(0, 1, 2, 3, 4, 5, 6))
    except:
        dt = np.dtype([('sspfile', '|S60'), ('age_base', '>f8'), ('Z_base', '>f8'),
                       ('Mstars', '>f8'), ('YA_V', '>i4'), ('aFe', '>f8') ])
        bdata = np.genfromtxt(bfile, dtype=dt, skip_header=1, usecols=(0, 1, 2, 4, 5, 6))

    base = Table(bdata)

    wl_ssp = np.genfromtxt(path.join(base_dir, base['sspfile'][0])).transpose()[0]

    base['base_wl'] = np.full((len(base),len(wl_ssp)), wl_ssp)
    base['base_spec'] = np.empty((len(base),len(wl_ssp)))

    for i in range(len(base)):
        base['base_spec'][i] = np.genfromtxt(path.join(base_dir, base['sspfile'][i])).transpose()[1]

    if save_fits:
        base.write(out_fname, overwrite=True)

    return base


