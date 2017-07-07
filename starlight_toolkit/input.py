import numpy as np

#TODO: Function to write grid files.


def write_pho_input(file_names, pho_dir, redshifts, globalchi2scalefactor
                    filter_names, filter_files, mags, errors
                    , ranges=np.zeros_like(filter_names)
                    , chi2scalefactors=np.full(1./len(filter_names))):

    N_PHO = len(filter_names)

    for i in range(len(file_names)):
        etc = open(pho_dir + file_names[i], 'w')
        etc.write('PHO\n' + str(redshifts[i]) + '\n' + str(N_PHO) + '\n' + \
        str(globalchi2scalefactor) + '\n')

        for filter in range(N_PHO):
            etc.write(filter_names[filter] + ' ' + filter_files[filter] + ' ' +\
            str(magnitudes[i,N_PHO]) + ' 1. ' + ' ' + \
            str(magnitudes[i,filter]) + ' ' + str(ranges[filter]) + ' ' \
            str(ranges[filter] + ' ' + str(chi2scalefactors[filter]) + '\n')

        etc.write('#PHO_name(iY) PHO_Filter_file(iY) PHO_logY_TOT(iY) PHO_YFrac2Model(iY) ')
        etc.write('PHO_ErrlogY(iY) PHO_RangelogY(iY) PHO_Chi2ScaleFactor(iY)')
        etc.close
