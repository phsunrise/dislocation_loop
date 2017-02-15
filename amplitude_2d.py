import numpy as np
import sys, os
from settings import basedir, MAXTIER, NFILES 
import time
from getopt import getopt
import glob
import pickle
from dampingfunc import dampingfunc

sample = 'W'
if sample == 'W':
    from W_parameters import R, rot, funcform, funcparams, orientations

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

if do_debug:
    ## choose data folder
    i_dir = 0
    while os.path.isdir(basedir + "%s_R%d_amp2d_%d/"%(sample, R, i_dir)):
        print "folder %s exists, containing %d files" % (\
                basedir+"%s_R%d_amp2d_%d/"%(sample, R, i_dir), \
                len(glob.glob(basedir+"%s_R%d_amp2d_%d/%s_atoms_amplitude_2d_*_T*_R*_ori*_????.npy"%(\
                    sample, R, i_dir, sample))))
        i_dir += 1
    if i_dir > 0:
        val = raw_input("enter folder number: ")
        i_dir = int(val)

    datadir = basedir + "%s_R%d_amp2d_%d/"%(sample, R, i_dir)
    if not os.path.isdir(datadir): ## if folder doesn't exist
        ## create folder
        os.system("mkdir %s" % datadir)
        print "created folder %s" % datadir
        ## save parameters
        os.system("cp %s_parameters.py %s" % (sample, datadir))
        print "copied file to %s%s_parameters.py" % (datadir, sample)
        os.system("cp qarray_2d.npz %s" % (datadir))
        print "copied file to %sqarray_2d.npz" % (datadir)
    else:
        import filecmp 
        if not filecmp.cmp("%s_parameters.py"%sample, datadir+"%s_parameters.py"%sample):
            print "%s_parameters.py file is different! Aborting..." % sample 
            sys.exit(1)
        if not filecmp.cmp("qarray_2d.npz", datadir+"qarray_2d.npz"):
            print "qarray_2d.npz file is different! Aborting..."
            sys.exit(1)

    looptypes = ['int', 'vac']
    _list = [(looptype, tier, i_ori, ori) for looptype in looptypes for tier in range(1, MAXTIER+1) 
                            for i_ori, ori in enumerate(orientations)]

    # get list of unprocessed data
    filelist = []
    s_datadir = basedir + "%s_R%d/"%(sample, R)
    for looptype, tier, i_ori, ori in _list:
        i_file = 0
        while os.path.isfile(s_datadir+"%s_atoms_s_%s_T%d_R%d_%04d_combined.npy"%(\
                  sample, looptype, tier, R, i_file)):
            if not os.path.isfile(\
                  datadir+"%s_atoms_amplitude_2d_%s_T%d_R%d_ori%d_%04d.npy"%(\
                    sample, looptype, tier, R, i_ori, i_file)):
                filelist.append([looptype, tier, i_ori, ori, i_file])
            i_file += 1
    #print filelist
    print "total %d files" % len(filelist)
    filelist.append(s_datadir)
    filelist.append(datadir)
    with open("amplitude_2d_filelist.pickle", 'wb') as fout:
        pickle.dump(filelist, fout)
    sys.exit(0)
else:
    filelist = pickle.load(open("amplitude_2d_filelist.pickle", 'rb'))
    s_datadir = filelist[-2] 
    datadir = filelist[-1] 
    filelist = filelist[:-2]

## load q array
_file = np.load("qarray_2d.npz")
q_array = _file['q_array']
K_array = _file['K_array']
h = _file['h']
for i_i_file, [looptype, tier, i_ori, ori, i_file] in enumerate(filelist):
    ## here is the parallelism criterion
    if i_i_file % nprocs != rank:
        continue

    ## read s data
    sdata = np.load(s_datadir+"%s_atoms_s_%s_T%d_R%d_%04d_combined.npy"%(sample, looptype, tier, R, i_file))

    amplitudes = np.zeros((q_array.shape[0], q_array.shape[1])).astype(complex) 
    #time0 = time.time()
    for index, _ in np.ndenumerate(amplitudes):
        #q = q_array[index]
        #qloop = np.einsum('ij,jk,k', ori, rot, q)
        K = K_array[index] 
        Kloop = np.einsum('ij,jk,k', ori, rot, K)

        ''' 
        no laue term or form factor calculation in this script
        '''
        r = sdata[:, 0:3]
        s = sdata[:, 3:6]
        r1 = r+s
        _temp = (np.exp(1.j*r1.dot(Kloop))-np.exp(1.j*r.dot(Kloop)))\
                    *dampingfunc(r, R, funcform, funcparams)
        #qr = r.dot(qloop)
        #Ks = s.dot(Kloop)
        ##hs = s.dot(hloop)
        #_temp = (np.cos(qr)*(np.cos(Ks)-1.)-\
        #    np.sin(qr)*np.sin(Ks))*dampingfunc(r1, R, funcform, funcparams)

        amplitudes[index] = np.sum(_temp)

        #print "done q=[%.2f, %.2f, %.2f], used time %f" % (
        #            q[0], q[1], q[2], time.time()-time0)
        #time0=time.time()

    # save data
    np.save(datadir+"%s_atoms_amplitude_2d_%s_T%d_R%d_ori%d_%04d.npy"%(\
                        sample, looptype, tier, R, i_ori, i_file), amplitudes)
