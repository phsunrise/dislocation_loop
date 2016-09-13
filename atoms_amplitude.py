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
q_array = np.linspace(-0.5, 0.5, 101)

D = 4.0*R ## parameter in gaussian used to smooth out the boundary

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
        while os.path.isfile(datadir+"%s_atoms_s_%s_T%d_R%d_%04d.npy"%(\
                  sample, looptype, tier, R, i_file)):
            if not os.path.isfile(\
                  "data/%s_atoms_amplitude_%s_T%d_R%d_ori%d_%04d.npy"%(\
                    sample, looptype, tier, R, i_ori, i_file)):
                filelist.append([looptype, tier, i_ori, ori, i_file])
            i_file += 1
    if do_debug:
        print filelist
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
        sdata = np.load("data/%s_atoms_s_%s_T%d_R%d_%04d.npy"%(sample, looptype, tier, R, i_file))
        if looptype == 'vac':
            pass
        elif looptype == 'int':
            sdata[:, 3:6] = -sdata[:, 3:6]

        amplitudes = [] 
        #amplitudes1 = []  # amplitude inside rho<5R, |z|<10R for comparison
        #amplitudes2 = []  # amplitude inside rho<10R, |z|<20R for comparison
        for i_qval, qval in enumerate(q_array):
            q = qabs * np.einsum('ij,j', ori, eq)
            qloop = np.einsum('ij,j', rot, q)
            K = np.einsum('ij,j', ori, h) + q
            Kloop = np.einsum('ij,j', rot, K)

            amplitude = 0.
            #amplitude1 = 0.
            #amplitude2 = 0.
            ''' 
            no laue term or form factor calculation in this script; these are taken care of 
            in the intensity script
            '''
            for line in sdata:
                r = line[0:3]
                s = line[3:6]
                r1 = r+s
                qr = qloop.dot(r)
                Ks = Kloop.dot(s)

                _temp = (np.cos(qr)*(np.cos(Ks)-1.)-\
                         np.sin(qr)*np.sin(Ks))*np.exp(-0.5*r1.dot(r1)/D**2)
                amplitude += _temp
                #if np.linalg.norm(r[0:2])<=50. and abs(r[2])<=100.:
                #    amplitude1 += _temp
                #if np.linalg.norm(r[0:2])<=100. and abs(r[2])<=200.:
                #    amplitude2 += _temp

            amplitudes.append([qval, amplitude])
            #amplitudes1.append([qval, amplitude1])
            #amplitudes2.append([qval, amplitude2])

        # save data
        np.save("data/%s_atoms_amplitude_%s_T%d_R%d_ori%d_%04d.npy"%(\
                            sample, looptype, tier, R, i_ori, i_file), amplitudes)
        #np.save("data/%s_atoms_amplitude1_%s_T%d_R%d_ori%d_%04d.npy"%(\
        #                    sample, looptype, tier, R, i_ori, i_file), amplitudes1)
        #np.save("data/%s_atoms_amplitude2_%s_T%d_R%d_ori%d_%04d.npy"%(\
        #                    sample, looptype, tier, R, i_ori, i_file), amplitudes2)
