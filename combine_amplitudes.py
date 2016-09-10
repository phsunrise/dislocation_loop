import numpy as np
import sys, os

data_dir = "W_R10_D_3.0R/"

sample = 'W'
looptypes = ['int', 'vac']
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *
elif sample == 'W':
    from W_parameters import *

_list = [(looptype, i_ori, ori) for looptype in looptypes for i_ori, ori in enumerate(orientations)]
for looptype, i_ori, ori in _list:
    amplitudes = np.load(data_dir+"%s_atoms_amplitude_%s_R%d_ori%d_0000.npy"%(\
                                sample, looptype, R, i_ori))
    amplitudes1 = np.load(data_dir+"%s_atoms_amplitude1_%s_R%d_ori%d_0000.npy"%(\
                                sample, looptype, R, i_ori))
    amplitudes2 = np.load(data_dir+"%s_atoms_amplitude2_%s_R%d_ori%d_0000.npy"%(\
                                sample, looptype, R, i_ori))
    i_file = 1 
    while True:
        try:
            amplitudes[:,1] += np.load(data_dir+"%s_atoms_amplitude_%s_R%d_ori%d_%04d.npy"%(\
                                sample, looptype, R, i_ori, i_file))[:,1]
            amplitudes1[:,1] += np.load(data_dir+"%s_atoms_amplitude1_%s_R%d_ori%d_%04d.npy"%(\
                                sample, looptype, R, i_ori, i_file))[:,1]
            amplitudes2[:,1] += np.load(data_dir+"%s_atoms_amplitude2_%s_R%d_ori%d_%04d.npy"%(\
                                sample, looptype, R, i_ori, i_file))[:,1]
            i_file += 1
        except IOError:
            break

    np.save(data_dir+"%s_atoms_amplitude_%s_R%d_ori%d_combined.npy"%(\
                        sample, looptype, R, i_ori), amplitudes)
    np.save(data_dir+"%s_atoms_amplitude1_%s_R%d_ori%d_combined.npy"%(\
                        sample, looptype, R, i_ori), amplitudes1)
    np.save(data_dir+"%s_atoms_amplitude2_%s_R%d_ori%d_combined.npy"%(\
                        sample, looptype, R, i_ori), amplitudes2)
