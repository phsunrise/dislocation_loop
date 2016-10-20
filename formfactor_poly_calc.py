import numpy as np
import sys, os
from settings import basedir, ff_dir
from formfactor import formfactor
import time

sample = 'W'
from W_parameters import *
looptypes = ['vac', 'int']
R = 20.
datadir = basedir + "W_R%d/"%R

# get uniform distributed points on sphere
sph_array = np.load("uniformsphere.npy")

print "sample:", sample
print "R=%.1f"%(R)
q_array = np.load(datadir+"q_array.npy")

for i_ori, ori in enumerate(orientations):
    print "starting orientation %d..." % i_ori

    # get the form factor
    ff_tier2 = []
    ff_tier3 = []
    for i_qval, qval in enumerate(q_array):
        time_start = time.time()
        q = qval * np.einsum('ij,j', ori, eq)
        K = np.einsum('ij,j', ori, h) + q
        Knorm = np.linalg.norm(K) 
        hloop = np.einsum('ij,jk,k', rot, ori, h)

        ff_tier2.append([])
        ff_tier3.append([])
        for th, phi in sph_array:
            Kloop = Knorm*np.array([np.sin(th)*np.cos(phi), np.sin(th)*np.sin(phi), np.cos(th)]) 
            qloop = Kloop-hloop
            ff_tier2[-1].append(formfactor(qloop, 2))
            ff_tier3[-1].append(formfactor(qloop, 3))
        time_end = time.time()
        print "done q=%f, used time %f s" % (qval, time_end-time_start)

    np.save(ff_dir+"ff_tier2_ori%d.npy"%(i_ori), ff_tier2)
    np.save(ff_dir+"ff_tier3_ori%d.npy"%(i_ori), ff_tier3)
