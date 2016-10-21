import numpy as np
import sys, os
from info import MAXTIER 
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

# see if q_array exists
try:
    q_array = np.load(datadir+"q_array.npy")
except IOError:
    q_array = np.linspace(-1.0, 1.0, 201)
    np.save(datadir+"q_array.npy", q_array)

D = 4.0*R ## parameter in gaussian used to smooth out the boundary

if __name__ == '__main__':
    from getopt import getopt
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

    _list = [(looptype, tier, i_ori, ori) for looptype in looptypes for tier in range(1, MAXTIER+1) 
                            for i_ori, ori in enumerate(orientations)]
    # first get list of unprocessed data
    filelist = []
    for looptype, tier, i_ori, ori in _list:
        i_file = 0
        while os.path.isfile(datadir+"%s_atoms_s_%s_T%d_R%d_%04d_combined.npy"%(\
                  sample, looptype, tier, R, i_file)):
            if not os.path.isfile(\
                  datadir+"%s_atoms_amplitude_poly_%s_T%d_R%d_ori%d_%04d.npy"%(\
                    sample, looptype, tier, R, i_ori, i_file)):
                filelist.append([looptype, tier, i_ori, ori, i_file])
            i_file += 1
    if do_debug:
        print filelist
        print "%d files" % len(filelist)

    # get uniform distributed points on sphere, save it in data folder
    sph_array = np.load("uniformsphere.npy")
    os.system("cp uniformsphere.npy %s"%(datadir))

    for i_i_file, [looptype, tier, i_ori, ori, i_file] in enumerate(filelist):
        ## here is the parallelism criterion
        if i_i_file % nprocs != rank:
            continue

        ## read s data
        sdata = np.load(datadir+"%s_atoms_s_%s_T%d_R%d_%04d_combined.npy"%(sample, looptype, tier, R, i_file))
        r = sdata[:, 0:3]
        s = sdata[:, 3:6]
        r1 = r+s

        amplitudes = [] 
        #amplitudes1 = []  # amplitude inside rho<50, |z|<100 for comparison
        #mask1 = np.logical_and((r[:,0]**2+r[:,1]**2)**0.5<=50., abs(r[:,2])<=100.)
        #amplitudes2 = []  # amplitude inside rho<100, |z|<200 for comparison
        #mask2 = np.logical_and((r[:,0]**2+r[:,1]**2)**0.5<=100., abs(r[:,2])<=200.)
        for i_qval, qval in enumerate(q_array):
            time_start = time.time()
            q = qval * np.einsum('ij,j', ori, eq)
            K = np.einsum('ij,j', ori, h) + q
            Knorm = np.linalg.norm(K) 
            hloop = np.einsum('ij,jk,k', rot, ori, h)

            amplitudes.append([])
            #amplitudes1.append([])
            #amplitudes2.append([])
            for th, phi in sph_array:
                Kloop = Knorm*np.array([np.sin(th)*np.cos(phi), np.sin(th)*np.sin(phi), np.cos(th)]) 
                Ks = s.dot(Kloop) 
                qr = r.dot(Kloop-hloop)
                ''' 
                no laue term or form factor calculation in this script; 
                these are taken care of in the intensity script
                '''
                _temp = (np.cos(qr)*(np.cos(Ks)-1.)-\
                    np.sin(qr)*np.sin(Ks))*np.exp(-0.5*np.sum(r1**2, axis=1)/D**2)
                if tier == 1:
                    _temp = _temp*sdata[:, 6]

                amplitudes[-1].append(np.sum(_temp))
                #amplitudes1[-1].append(np.sum(_temp[mask1]))
                #amplitudes2[-1].append(np.sum(_temp[mask2]))

            time_end = time.time()
            print "used time %f for q=%f" % (time_end-time_start, qval)
            if do_debug and i_qval == 5:
                sys.exit(0)

        # save data
        np.save(datadir+"%s_atoms_amplitude_poly_%s_T%d_R%d_ori%d_%04d.npy"%(\
                            sample, looptype, tier, R, i_ori, i_file), amplitudes)
        #np.save(datadir+"%s_atoms_amplitude1_poly_%s_T%d_R%d_ori%d_%04d.npy"%(\
        #                    sample, looptype, tier, R, i_ori, i_file), amplitudes1)
        #np.save(datadir+"%s_atoms_amplitude2_poly_%s_T%d_R%d_ori%d_%04d.npy"%(\
        #                    sample, looptype, tier, R, i_ori, i_file), amplitudes2)
