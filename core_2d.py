import numpy as np
import sys, os
from settings import basedir
from getopt import getopt
import glob
import pickle
from dampingfunc import dampingfunc
from ampint import ampint

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
from W_parameters import R

print "sample:", sample

if do_debug:
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

    del sys.modules["W_parameters"]
    sys.path.insert(0, datadir)
    from W_parameters import orientations 
    looptypes = ['vac', 'int']
    _list = [(looptype, i_ori, ori) for looptype in looptypes \
                            for i_ori, ori in enumerate(orientations)]
    _list.append(datadir)
    with open("core_2d_filelist.pickle", 'wb') as fout:
        pickle.dump(_list, fout)
    sys.exit(0)

else:
    _list = pickle.load(open("core_2d_filelist.pickle", 'rb'))
    datadir = _list[-1]
    _list = _list[:-1]

    del sys.modules["W_parameters"]
    sys.path.insert(0, datadir)
    from W_parameters import R, funcform, funcparams, rot, a0, a1, ex_p, ey_p, C12, C44, Bloop 


_file = np.load(datadir+"qarray_2d.npz")
q_array = _file['q_array'] 
K_array = _file['K_array']
h = _file['h']
print "h = 2pi/a0 *", h*a0/2./np.pi 

for _i_list, (looptype, i_ori, ori) in enumerate(_list):
    if _i_list % nprocs != rank:
        continue
    print "starting looptype %s, orientation %d..." % (looptype, i_ori)
    #if not (i_ori==2 and looptype=='int'):
    #    continue
    print "starting orientation %d..." % i_ori
    amplitudes = np.zeros((q_array.shape[0], q_array.shape[1])).astype(complex)
    amplitude_ints = np.zeros((q_array.shape[0], q_array.shape[1])).astype(complex)
    counter = 0
    for index, _ in np.ndenumerate(amplitudes):
        q = q_array[index]
        qloop = np.einsum('ij,jk,k', ori, rot, q)
        K = K_array[index] 
        Kloop = np.einsum('ij,jk,k', ori, rot, K)
        
        amplitude = 0.
        '''
        calculate the laue term
        The dislocation loop is inside the region r <= R
        '''
        for x_p in np.arange(-np.ceil(2.*R/a1), np.ceil(2.*R/a1)):
            for y_p in np.arange(-np.ceil(2.*R/a1), np.ceil(2.*R/a1)):
                r = x_p*ex_p+y_p*ey_p
                if np.linalg.norm(r) <= R:
                    if looptype == 'int':
                        amplitude += np.exp(1.j*Kloop.dot(r))\
                                    *dampingfunc(r, R, funcform, funcparams)
                    elif looptype == 'vac':
                        amplitude -= np.exp(1.j*Kloop.dot(r))\
                                    *dampingfunc(r, R, funcform, funcparams)
        amplitudes[index] = amplitude
        ## add the integral to infinity
        rcutoff = 400.
        if looptype == 'int':
            amplitude_ints[index] = -ampint(qvec=qloop, Kvec=Kloop, \
                    bvec=Bloop, R=R, r0=rcutoff, C12=C12, C44=C44, a0=a0)
        else:
            amplitude_ints[index] = ampint(qvec=qloop, Kvec=Kloop, \
                    bvec=Bloop, R=R, r0=rcutoff, C12=C12, C44=C44, a0=a0)

        counter += 1
        if counter % 200 == 0:
            print "done %.2f%%\r" % (counter*100./np.prod(amplitudes.shape)),
            sys.stdout.flush()

    print "done 100.00%"
    #print amplitudes
    np.save(datadir+"%s_atoms_core_2d_%s_R%d_ori%d.npy"%(\
                sample, looptype, R, i_ori), amplitudes)
    np.save(datadir+"%s_atoms_ampint_2d_%s_R%d_ori%d.npy"%(\
                sample, looptype, R, i_ori), amplitude_ints)
