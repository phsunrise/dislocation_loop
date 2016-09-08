'''
This code applies only to fcc crystals with {111} dislocation loops
The script is mostly written in loop coordinates
'''
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from info import *

sample = 'Cu'
looptypes = ['int', 'vac']
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *
elif sample == 'W':
    from W_parameters import *

qR_array = R*np.linspace(-0.5, 0.5, 101)

D = 4.0*R ## parameter in gaussian used to smooth out the boundary

if __name__ == '__main__':
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nprocs = comm.Get_size()

    _list = [(looptype, i_ori, ori) for looptype in looptypes for i_ori, ori in enumerate(orientations)]
    filelists = []
    for looptype, i_ori, ori in _list:
        # first get list of unprocessed data
        filelists.append([])
        for i_file in xrange(NFILES_max):
            if os.path.isfile("data/%s_atoms_s_%s_R%d_%04d.npy"%(\
                  sample, looptype, R, i_file)) and not os.path.isfile(\
                  "data/%s_atoms_amplitude_%s_R%d_ori%d_%04d.npy"%(\
                    sample, looptype, R, i_ori, i_file)):
                filelists[-1].append(i_file)
    #print filelists
    #sys.exit(0)

    for (looptype, i_ori, ori), filelist in zip(_list, filelists):
        print "Processing looptype %s, orientation %d" % (looptype, i_ori)
        for i_i_file, i_file in enumerate(filelist):
            ## here is the parallelism criterion
            if i_i_file % nprocs != rank:
                continue

            ## read s data
            sdata = np.load("data/%s_atoms_s_%s_R%d_%04d.npy"%(sample, looptype, R, i_file))

            amplitudes = [] 
            amplitudes1 = []  # amplitude inside rho<5R, |z|<10R for comparison
            amplitudes2 = []  # amplitude inside rho<10R, |z|<20R for comparison
            for i_qR, qR in enumerate(qR_array):
                q = qR/R * np.einsum('ij,j', ori, eq)
                qloop = np.einsum('ij,j', rot, q)
                K = np.einsum('ij,j', ori, h) + q
                Kloop = np.einsum('ij,j', rot, K)

                amplitude = 0.
                amplitude1 = 0.
                amplitude2 = 0.
                ''' 
                no laue term calculation in this script; this is taken care of 
                in the intensity script
                '''
                for x, y, z, sx, sy, sz in sdata:
                    r = np.array([x,y,z])
                    qr = qloop.dot(r)
                    if looptype == 'vac':
                        Ks = Kloop.dot([sx, sy, sz])
                    elif looptype == 'int':
                        Ks = Kloop.dot([-sx, -sy, -sz])

                    _temp = (np.cos(qr)*(np.cos(Ks)-1.)-\
                             np.sin(qr)*np.sin(Ks))*np.exp(-0.5*r.dot(r)/D**2)
                    amplitude += _temp
                    if np.sqrt(x**2+y**2)<=50. and abs(z)<=100.:
                        amplitude1 += _temp
                    if np.sqrt(x**2+y**2)<=100. and abs(z)<=200.:
                        amplitude2 += _temp

                amplitudes.append([qR/R, amplitude])
                amplitudes1.append([qR/R, amplitude1])
                amplitudes2.append([qR/R, amplitude2])

            # save data
            np.save("data/%s_atoms_amplitude_%s_R%d_ori%d_%04d.npy"%(\
                                sample, looptype, R, i_ori, i_file), amplitudes)
            np.save("data/%s_atoms_amplitude1_%s_R%d_ori%d_%04d.npy"%(\
                                sample, looptype, R, i_ori, i_file), amplitudes1)
            np.save("data/%s_atoms_amplitude2_%s_R%d_ori%d_%04d.npy"%(\
                                sample, looptype, R, i_ori, i_file), amplitudes2)
