import numpy as np
from scipy.integrate import dblquad 
import sys, os
from datetime import datetime
from getopt import getopt
import pickle
from settings import MAXTIER, NFILES, basedir, preproc_dir 
from displacement_iso import disp 

do_debug = False
do_mpi = True

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
looptypes = ['vac', 'int']
from W_parameters import R, a0, C12, C44, Bloop 

datadir = basedir+"%s_R%d/" % (sample, R)
if not os.path.isdir(datadir):
    print "creating directory %s ..."%datadir
    os.system("mkdir %s"%datadir)

# get uncalculated files
if do_debug:
    filelist = []
    for looptype in looptypes:
        for tier in range(1, MAXTIER+1):
            for i_file in xrange(NFILES):
                if os.path.isfile(preproc_dir+"%s_atoms_s_%s_pre_T%d_%04d.npy"%(\
                   sample, looptype, tier, i_file)) and not os.path.isfile(datadir+"%s_atoms_s_%s_T%d_R%d_%04d.npy"%(sample, looptype, tier, R, i_file)):
                    filelist.append([looptype, tier, i_file])
    print filelist
    print "%d files" % len(filelist)
    with open("s_iso_filelist.pickle", 'wb') as fout:
        pickle.dump(filelist, fout)
    sys.exit(0)
else:
    filelist = pickle.load(open("s_iso_filelist.pickle", 'rb'))

for i_i_file, [looptype, tier, i_file] in enumerate(filelist):
    if i_i_file % nprocs != rank:
        continue

    xyz_list = np.load(preproc_dir+"%s_atoms_s_%s_pre_T%d_%04d.npy"%(\
                            sample, looptype, tier, i_file))
    data = []
    for i_xyz, xyz in enumerate(xyz_list):
        x, y, z = xyz[0:3]/R
        if looptype == 'vac' and abs(z) <= 1.e-6: 
            if abs(x**2+y**2) <= 1.: # neglect atoms in the core
                continue
            else:
                s = 0.5*(disp(Bloop, x, y, 1.e-6, 1., C12, C44)\
                        +disp(Bloop, x, y, -1.e-6, 1., C12, C44))
        else:
            s = disp(Bloop, x, y, z, 1., C12, C44)
        data.append([xyz[0], xyz[1], xyz[2], s[0], s[1], s[2]])
        
        if i_xyz % 10 == 0:
            print "done %s, tier %d, file %04d, entry %04d, at %s" % (\
                        looptype, tier, i_file, i_xyz, datetime.now())

    data = np.array(data)
    if looptype == 'int':
        data[:,3:6] = -data[:,3:6]
    np.save(datadir+"%s_atoms_s_%s_T%d_R%d_%04d.npy"%(sample, looptype, tier, R, i_file), data)
    print "done %s, tier %d, file %04d, at %s" % (\
                looptype, tier, i_file, datetime.now())
