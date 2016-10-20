import numpy as np
import sys, os
from info import MAXTIER 
from settings import basedir

sample = 'W'
looptypes = ['int', 'vac']
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *
elif sample == 'W':
    from W_parameters import *
R = 20.
datadir = basedir+"%s_R%d/"%(sample, R)


_list = [(looptype, tier, i_ori, ori) for looptype in looptypes for tier in range(1, MAXTIER+1) \
            for i_ori, ori in enumerate(orientations)]
for looptype, tier, i_ori, ori in _list:
    amplitudes = np.load(datadir+"%s_atoms_amplitude_poly_%s_T%d_R%d_ori%d_0000.npy"%(\
                                sample, looptype, tier, R, i_ori))
    #amplitudes1 = np.load(datadir+"%s_atoms_amplitude1_poly_%s_T%d_R%d_ori%d_0000.npy"%(\
    #                            sample, looptype, tier, R, i_ori))
    #amplitudes2 = np.load(datadir+"%s_atoms_amplitude2_poly_%s_T%d_R%d_ori%d_0000.npy"%(\
    #                            sample, looptype, tier, R, i_ori))
    i_file = 1 
    while True:
        try:
            amplitudes += np.load(datadir+"%s_atoms_amplitude_poly_%s_T%d_R%d_ori%d_%04d.npy"%(\
                                sample, looptype, tier, R, i_ori, i_file))
            #amplitudes1 += np.load(datadir+"%s_atoms_amplitude1_poly_%s_T%d_R%d_ori%d_%04d.npy"%(\
            #                    sample, looptype, tier, R, i_ori, i_file))
            #amplitudes2 += np.load(datadir+"%s_atoms_amplitude2_poly_%s_T%d_R%d_ori%d_%04d.npy"%(\
            #                    sample, looptype, tier, R, i_ori, i_file))
            i_file += 1
        except IOError:
            break

    np.save(datadir+"%s_atoms_amplitude_poly_%s_T%d_R%d_ori%d_combined.npy"%(\
                        sample, looptype, tier, R, i_ori), amplitudes)
    #np.save(datadir+"%s_atoms_amplitude1_poly_%s_T%d_R%d_ori%d_combined.npy"%(\
    #                    sample, looptype, tier, R, i_ori), amplitudes1)
    #np.save(datadir+"%s_atoms_amplitude2_poly_%s_T%d_R%d_ori%d_combined.npy"%(\
    #                    sample, looptype, tier, R, i_ori), amplitudes2)
    print "done %s, tier %d, orientation %d, total %d files" % (\
                looptype, tier, i_ori, i_file)

val = raw_input("Delete uncombined files? ")
if val in ['y', 'Y', 'yes']:
    for looptype, tier, i_ori, ori in _list:
        os.system("rm "+datadir+"%s_atoms_amplitude_poly_%s_T%d_R%d_ori%d_????.npy"%(\
                                sample, looptype, tier, R, i_ori))
        print "Deleted %s, tier %d, orientation %d, total %d files" % (\
                    looptype, tier, i_ori, i_file)
