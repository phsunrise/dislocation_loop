import numpy as np
import sys, os
from info import MAXTIER, NFILES 
from settings import basedir
import time

sample = 'W'
looptypes = ['int', 'vac']
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *
elif sample == 'W':
    from W_parameters import *

datadir = basedir + "%s_R%d/"%(sample, R)
q_array = np.linspace(-1.0, 1.0, 201)

D = 4.*R ## parameter in gaussian used to smooth out the boundary

if __name__ == '__main__':
    from getopt import getopt
    do_debug = False

    opts, args = getopt(sys.argv[1:], "d")
    for o, a in opts:
        if o == '-d':
            do_debug = True

    _list = [(looptype, tier, i_ori, ori) for looptype in looptypes for tier in range(1, MAXTIER+1) 
                            for i_ori, ori in enumerate(orientations)]
    # first get list of unprocessed data
    filelist = []
    for looptype, tier, i_ori, ori in _list:
        i_file = 0
        while os.path.isfile(datadir+"%s_atoms_s_%s_T%d_R%d_%04d_combined.npy"%(\
                  sample, looptype, tier, R, i_file)):
            if not os.path.isfile(\
                  datadir+"%s_atoms_amplitude_%s_T%d_R%d_ori%d_%04d.npy"%(\
                    sample, looptype, tier, R, i_ori, i_file)):
                filelist.append([looptype, tier, i_ori, ori, i_file])
            i_file += 1
    if do_debug:
        print filelist
        print "total %d files" % len(filelist)
        sys.exit(0)

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nprocs = comm.Get_size()

    for i_i_file, [looptype, tier, i_ori, ori, i_file] in enumerate(filelist):
        ## here is the parallelism criterion
        if i_i_file % nprocs != rank:
            continue

        ## read s data
        sdata = np.load(datadir+"%s_atoms_s_%s_T%d_R%d_%04d_combined.npy"%(sample, looptype, tier, R, i_file))

        amplitudes = [] 
        amplitudes1 = []  # amplitude inside rho<50, |z|<100 for comparison
        amplitudes2 = []  # amplitude inside rho<100, |z|<200 for comparison
        #time0 = time.time()
        for i_qval, qval in enumerate(q_array):
            q = qval * np.einsum('ij,j', ori, eq)
            qloop = np.einsum('ij,j', rot, q)
            K = np.einsum('ij,j', ori, h) + q
            Kloop = np.einsum('ij,j', rot, K)

            ''' 
            no laue term or form factor calculation in this script; 
            these are taken care of in the intensity script
            '''
            r = sdata[:, 0:3]
            s = sdata[:, 3:6]
            r1 = r+s
            qr = r.dot(qloop)
            Ks = s.dot(Kloop)
            _temp = (np.cos(qr)*(np.cos(Ks)-1.)-\
                np.sin(qr)*np.sin(Ks))*np.exp(-0.5*np.sum(r1**2, axis=1)/D**2)
            if tier == 1:
                _temp = _temp*sdata[:, 6]

            amplitude = np.sum(_temp)
            amplitude1 = np.sum(_temp[np.logical_and(
                            (r[:,0]**2+r[:,1]**2)**0.5<=50., abs(r[:,2])<=100.)])
            amplitude2 = np.sum(_temp[np.logical_and(
                            (r[:,0]**2+r[:,1]**2)**0.5<=100., abs(r[:,2])<=200.)])

            amplitudes.append([qval, amplitude])
            amplitudes1.append([qval, amplitude1])
            amplitudes2.append([qval, amplitude2])

            #print "done q=%.2f, used time %f" % (qval, time.time()-time0)
            #time0=time.time()

        # save data
        np.save(datadir+"%s_atoms_amplitude_%s_T%d_R%d_ori%d_%04d.npy"%(\
                            sample, looptype, tier, R, i_ori, i_file), amplitudes)
        np.save(datadir+"%s_atoms_amplitude1_%s_T%d_R%d_ori%d_%04d.npy"%(\
                            sample, looptype, tier, R, i_ori, i_file), amplitudes1)
        np.save(datadir+"%s_atoms_amplitude2_%s_T%d_R%d_ori%d_%04d.npy"%(\
                            sample, looptype, tier, R, i_ori, i_file), amplitudes2)
