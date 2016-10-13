import numpy as np
from info import *

R = 40

n = 32 # number of subfiles to combine

for looptype in ['vac', 'int']:
    for i_file in xrange(300):
        try:
            sdata = np.load("W_R%d/W_atoms_s_%s_T1_R%d_%04d.npy"%(\
                    R, looptype, R, i_file*n))
        except IOError:
            print "done %s!" % looptype
            break
        for j in range(1,n):
            try:
                sdata = np.vstack((sdata, np.load(\
                    "W_R%d/W_atoms_s_%s_T1_R%d_%04d.npy"%(\
                    R, looptype, R, i_file*n+j))))
            except IOError:
                break

        np.save("W_R%d/W_atoms_s_%s_T1_R%d_%04d_combined.npy"%(\
            R, looptype, R, i_file), sdata)
        print "done file %d" % i_file
