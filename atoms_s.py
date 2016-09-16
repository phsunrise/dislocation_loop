'''
Adaptive quadrature method
'''
import numpy as np
from scipy.integrate import dblquad 
import sys, os
from datetime import datetime
from getopt import getopt
from info import *

do_debug = False

opts, args = getopt(sys.argv[1:], "d")
for o, a in opts:
    if o == '-d':
        do_debug = True

sample = 'W'
looptypes = ['vac', 'int']
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *
elif sample == 'W':
    from W_parameters import *

datadir = "%s_R%d/" % (sample, R)

# get uncalculated files
filelist = []
for looptype in looptypes:
    for tier in range(1, MAXTIER+1):
        for i_file in xrange(NFILES):
            if os.path.isfile("preproc/%s_atoms_s_%s_pre_T%d_%04d.npy"%(\
               sample, looptype, tier, i_file)) and not os.path.isfile(datadir+"%s_atoms_s_%s_T%d_R%d_%04d.npy"%(sample, looptype, tier, R, i_file)):
                filelist.append([looptype, tier, i_file])
if do_debug:
    print filelist
    print "%d files" % len(filelist)
    sys.exit(0)

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
for i_i_file, [looptype, tier, i_file] in enumerate(filelist):
    if i_i_file % nprocs != rank:
        continue

    xyz_list = np.load("preproc/%s_atoms_s_%s_pre_T%d_%04d.npy"%(\
                            sample, looptype, tier, i_file))
    data = []
    for i_xyz, xyz in enumerate(xyz_list):
        x, y, z = xyz[0:3]/R

        def func(th, th1, ind):
            a = -(x*np.cos(th1)+y*np.sin(th1)+1.)/z
            b = -(x*np.cos(th1)+y*np.sin(th1)-1.)/z
            t = (b+a)/2.+(b-a)*np.cos(th)/2.
            e = (np.cos(th1)*ex+np.sin(th1)*ey+t*ez)/(1.+t**2)**0.5
            eloop = np.array([np.cos(th1),np.sin(th1),t])/(1.+t**2)**0.5
            g = np.zeros((3,3))
            sum1 = sum(e**2/(C44+d*e**2))
            for i in xrange(3):
                for j in xrange(3):
                    g[i,j] = ((1 if i==j else 0)/(C44+d*e[i]**2) - \
                        (e[i]*e[j]/(C44+d*e[i]**2)/(C44+d*e[j]**2)*\
                         (C44+C12)/(1.+(C44+C12)*sum1)
                        )
                      )
            gloop =  np.einsum('ij,kl,jl', rot, rot, g)
            return (eloop.dot(np.array([x,y,z]))*\
                np.einsum('ij,i,jk', Ploop, eloop, gloop))[ind]

        integrand0, err = dblquad(func, 0., 2.*np.pi, \
                lambda x: -np.pi, lambda x: np.pi, args=(0,))
        integrand1, err = dblquad(func, 0., 2.*np.pi, \
                lambda x: -np.pi, lambda x: np.pi, args=(1,))
        integrand2, err = dblquad(func, 0., 2.*np.pi, \
                lambda x: -np.pi, lambda x: np.pi, args=(2,))

        s = -1./(8.*np.pi**3*R**2*abs(z))*np.array(\
                [integrand0, integrand1, integrand2])

        if tier == 1:
            data.append([xyz[0], xyz[1], xyz[2], s[0], s[1], s[2], xyz[3]])
        else:
            data.append([xyz[0], xyz[1], xyz[2], s[0], s[1], s[2]])
        
        if i_file % 10 == 0:
            print "done %s, tier %d, file %04d, entry %04d, at %s" % (\
                        looptype, tier, i_file, i_xyz, datetime.now())

    data = np.array(data)
    if looptype == 'int':
        data[:,3:6] = -data[:,3:6]
    np.save(datadir+"%s_atoms_s_%s_T%d_R%d_%04d.npy"%(sample, looptype, tier, R, i_file), data)
