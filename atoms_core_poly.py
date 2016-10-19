import numpy as np
import sys, os
from settings import basedir
from getopt import getopt
import time

do_debug = False
do_mpi = True # when -n and -p commands are supplied, will do "fake mpi"

opts, args = getopt(sys.argv[1:], "dn:p:")
for o, a in opts:
    if o == '-d':
        do_debug = True
    elif o == '-n':
        nprocs = int(a)
        do_mpi = False
    elif o == '-p':
        rank = int(a)

if do_mpi:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nprocs = comm.Get_size()

sample = 'W'
from W_parameters import *
looptypes = ['vac', 'int']

print "sample:", sample
print "R=%.1f, D=%.1f*R"%(R, D/R)
datadir = basedir+"%s_R%d/"%(sample, R)
q_array = np.load(datadir+"q_array.npy")
# get uniform distributed points on sphere
sph_array = np.load("uniformsphere.npy")

_list = [(looptype, i_ori, ori) for looptype in looptypes \
                        for i_ori, ori in enumerate(orientations)]
for _i_list, (looptype, i_ori, ori) in enumerate(_list):
    if _i_list % nprocs != rank:
        continue
    print "starting looptype %s, orientation %d..." % (looptype, i_ori)
    amplitudes = [] 
    for i_qval, qval in enumerate(q_array):
        time_start = time.time()
        q = qval * np.einsum('ij,j', ori, eq)
        K = np.einsum('ij,j', ori, h) + q
        Knorm = np.linalg.norm(K) 
        
        amplitudes.append([])
        for th, phi in sph_array:
            amplitude = 0.
            Kloop = Knorm*np.array([np.sin(th)*np.cos(phi), \
                            np.sin(th)*np.sin(phi), np.cos(th)]) 
            '''
            calculate the laue term
            The dislocation loop is inside the region r <= R
            '''
            for x_p in np.arange(-np.ceil(2.*R/a1), np.ceil(2.*R/a1)):
                for y_p in np.arange(-np.ceil(2.*R/a1), np.ceil(2.*R/a1)):
                    r = x_p*ex_p+y_p*ey_p
                    if np.linalg.norm(r) <= R:
                        if looptype == 'int':
                            amplitude += np.exp(\
                                -0.5*r.dot(r)/D**2)*np.cos(Kloop.dot(r))
                        elif looptype == 'vac':
                            amplitude -= np.exp(\
                                -0.5*r.dot(r)/D**2)*np.cos(Kloop.dot(r))
            amplitudes[-1].append(amplitude)
        time_end = time.time()
        print "done q=%f, used time %f s" % (qval, time_end-time_start)
        if do_debug and i_qval == 5:
            sys.exit(0)

    np.save(datadir+"%s_atoms_core_poly_%s_R%d_ori%d.npy"%(\
                sample, looptype, R, i_ori), amplitudes)
