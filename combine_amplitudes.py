import numpy as np
import sys, os
from info import *

sample = 'W'
looptypes = ['int', 'vac']
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *
elif sample == 'W':
    from W_parameters import *
datadir = "%s_R%d/" % (sample, R)


_list = [(looptype, tier, i_ori, ori) for looptype in looptypes for tier in range(1, MAXTIER+1) \
            for i_ori, ori in enumerate(orientations)]
for looptype, tier, i_ori, ori in _list:
    amplitudes = np.load(datadir+"%s_atoms_amplitude_%s_T%d_R%d_ori%d_0000.npy"%(\
                                sample, looptype, tier, R, i_ori))
    #amplitudes1 = np.load(datadir+"%s_atoms_amplitude1_%s_R%d_ori%d_0000.npy"%(\
    #                            sample, looptype, R, i_ori))
    #amplitudes2 = np.load(datadir+"%s_atoms_amplitude2_%s_R%d_ori%d_0000.npy"%(\
    #                            sample, looptype, R, i_ori))
    i_file = 1 
    while True:
        try:
            amplitudes[:,1] += np.load(datadir+"%s_atoms_amplitude_%s_T%d_R%d_ori%d_%04d.npy"%(\
                                sample, looptype, tier, R, i_ori, i_file))[:,1]
            #amplitudes1[:,1] += np.load(data_dir+"%s_atoms_amplitude1_%s_R%d_ori%d_%04d.npy"%(\
            #                    sample, looptype, R, i_ori, i_file))[:,1]
            #amplitudes2[:,1] += np.load(data_dir+"%s_atoms_amplitude2_%s_R%d_ori%d_%04d.npy"%(\
            #                    sample, looptype, R, i_ori, i_file))[:,1]
            i_file += 1
        except IOError:
            break

    np.save(datadir+"%s_atoms_amplitude_%s_T%d_R%d_ori%d_combined.npy"%(\
                        sample, looptype, tier, R, i_ori), amplitudes)
    #np.save(datadir+"%s_atoms_amplitude1_%s_R%d_ori%d_combined.npy"%(\
    #                    sample, looptype, R, i_ori), amplitudes1)
    #np.save(datadir+"%s_atoms_amplitude2_%s_R%d_ori%d_combined.npy"%(\
    #                    sample, looptype, R, i_ori), amplitudes2)
    print "done %s, tier %d, orientation %d" % (looptype, tier, i_ori)
