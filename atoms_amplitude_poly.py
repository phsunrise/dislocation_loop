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
                  datadir+"%s_atoms_amplitude_poly_%s_T%d_R%d_ori%d_%04d.npy"%(\
                    sample, looptype, tier, R, i_ori, i_file)):
                filelist.append([looptype, tier, i_ori, ori, i_file])
            i_file += 1
    if do_debug:
        print filelist
        print "%d files" % len(filelist)
        sys.exit(0)

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nprocs = comm.Get_size()

    # get uniform distributed points on sphere
    sph_array = np.load("uniformsphere.npy")

    for i_i_file, [looptype, tier, i_ori, ori, i_file] in enumerate(filelist):
        ## here is the parallelism criterion
        if i_i_file % nprocs != rank:
            continue

        ## read s data
        sdata = np.load(datadir+"%s_atoms_s_%s_T%d_R%d_%04d_combined.npy"%(sample, looptype, tier, R, i_file))

        amplitudes = [] 
        amplitudes1 = []  # amplitude inside rho<50, |z|<100 for comparison
        amplitudes2 = []  # amplitude inside rho<100, |z|<200 for comparison
        for i_qval, qval in enumerate(q_array):
            q = qval * np.einsum('ij,j', ori, eq)
            #qloop = np.einsum('ij,j', rot, q)
            K = np.einsum('ij,j', ori, h) + q
            Kloop = np.einsum('ij,j', rot, K) # this K is held fixed for different crystallines;
                                              # h will be rotated and q will change accordingly

            amplitudes.append([])
            amplitudes1.append([])
            amplitudes2.append([])
            for th, phi in sph_array:
                amplitude = 0.
                amplitude1 = 0.
                amplitude2 = 0.
                ''' 
                no laue term or form factor calculation in this script; 
                these are taken care of in the intensity script
                '''
                for line in sdata:
                    r = line[0:3]
                    s = line[3:6]
                    r1 = r+s
                    Ks = Kloop.dot(s)
                    hloop = hnorm*np.array([np.sin(th)*np.cos(phi), np.sin(th)*np.sin(phi), np.cos(th)]) 
                    qr = (Kloop-hloop).dot(r)

                    _temp = (np.cos(qr)*(np.cos(Ks)-1.)-\
                        np.sin(qr)*np.sin(Ks))*np.exp(-0.5*r1.dot(r1)/D**2)
                    if tier == 1:
                        _temp *= line[6]
                    amplitude += _temp
                    if np.linalg.norm(r[0:2])<=50. and abs(r[2])<=100.:
                        amplitude1 += _temp
                    if np.linalg.norm(r[0:2])<=100. and abs(r[2])<=200.:
                        amplitude2 += _temp

                amplitudes[-1].append(amplitude)
                amplitudes1[-1].append(amplitude1)
                amplitudes2[-1].append(amplitude2)

        # save data
        np.save(datadir+"%s_atoms_amplitude_poly_%s_T%d_R%d_ori%d_%04d.npy"%(\
                            sample, looptype, tier, R, i_ori, i_file), amplitudes)
        np.save(datadir+"%s_atoms_amplitude1_poly_%s_T%d_R%d_ori%d_%04d.npy"%(\
                            sample, looptype, tier, R, i_ori, i_file), amplitudes1)
        np.save(datadir+"%s_atoms_amplitude2_poly_%s_T%d_R%d_ori%d_%04d.npy"%(\
                            sample, looptype, tier, R, i_ori, i_file), amplitudes2)
