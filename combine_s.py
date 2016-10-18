import numpy as np
from info import MAXTIER
from settings import basedir


R = 20.
datadir = basedir + "W_R%d/"%(R)

n = 64 # number of subfiles to combine

_list = [(looptype, tier) for looptype in ['vac', 'int'] \
            for tier in range(1, MAXTIER+1)]
for looptype, tier in _list:
    for i_file in xrange(1000):
        try:
            sdata = np.load(datadir+"W_atoms_s_%s_T%d_R%d_%04d.npy"%(\
                    looptype, tier, R, i_file*n))
        except IOError:
            print "done %s!" % looptype
            break
        for j in range(1,n):
            try:
                sdata = np.vstack((sdata, np.load(\
                    datadir+"W_atoms_s_%s_T%d_R%d_%04d.npy"%(\
                    looptype, tier, R, i_file*n+j))))
            except IOError:
                break

        np.save(datadir+"W_atoms_s_%s_T%d_R%d_%04d_combined.npy"%(\
            looptype, tier, R, i_file), sdata)
        print "done file %d" % i_file
