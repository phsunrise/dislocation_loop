import numpy as np
import sys, os
from settings import basedir, MAXTIER
import glob

sample = 'W'
looptypes = ['int', 'vac']
from W_parameters import R, orientations 
i_dir = 0
while os.path.isdir(basedir + "%s_R%d_amp2d_%d/"%(sample, R, i_dir)):
    print "folder %s exists, containing %d files" % (\
            basedir+"%s_R%d_amp2d_%d/"%(sample, R, i_dir), \
            len(glob.glob(basedir+"%s_R%d_amp2d_%d/%s_atoms_amplitude_2d_*_T*_R*_ori*_????.npy"%(\
                sample, R, i_dir, sample))))
    i_dir += 1
val = raw_input("enter folder number: ")
i_dir = int(val)
datadir = basedir+"%s_R%d_amp2d_%d/"%(sample, R, i_dir)

_list = [(looptype, tier, i_ori, ori) for looptype in looptypes for tier in range(1, MAXTIER+1) \
            for i_ori, ori in enumerate(orientations)]
for looptype, tier, i_ori, ori in _list:
    amplitudes = np.load(datadir+"%s_atoms_amplitude_2d_%s_T%d_R%d_ori%d_0000.npy"%(\
                                sample, looptype, tier, R, i_ori))
    i_file = 1 
    while True:
        try:
            amplitudes += np.load(datadir+"%s_atoms_amplitude_2d_%s_T%d_R%d_ori%d_%04d.npy"%(\
                                sample, looptype, tier, R, i_ori, i_file))
            i_file += 1
        except IOError:
            break

    np.save(datadir+"%s_atoms_amplitude_2d_%s_T%d_R%d_ori%d_combined.npy"%(\
                        sample, looptype, tier, R, i_ori), amplitudes)
    print "done %s, tier %d, orientation %d, total %d files" % (\
                looptype, tier, i_ori, i_file)

#val = raw_input("Delete uncombined files? ")
#if val in ['y', 'Y', 'yes']:
#    for looptype, tier, i_ori, ori in _list:
#        os.system("rm "+datadir+"%s_atoms_amplitude_%s_T%d_R%d_ori%d_????.npy"%(\
#                                sample, looptype, tier, R, i_ori))
#        print "Deleted %s, tier %d, orientation %d, total %d files" % (\
#                    looptype, tier, i_ori, i_file)
