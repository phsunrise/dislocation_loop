import numpy as np
import sys, os
from getopt import getopt

do_debug = False
MAXTIER = 3

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

# get number of files in "preproc"
NFILES = 0 
while os.path.isfile("preproc/%s_atoms_s_vac_pre_T1_%04d.npy"%(\
                        sample, NFILES)):
    NFILES += 1

# get uncalculated files
filelist = []
for looptype in looptypes:
    for tier in range(1, MAXTIER+1):
        for i_file in xrange(NFILES):
            if not os.path.isfile("data/%s_atoms_s_%s_T%d_R%d_%04d.npy"%(\
                  sample, looptype, tier, R, i_file)):
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
    for xyz in xyz_list:
        x, y, z = xyz[0:3]/R
        n = (30 if x**2+y**2+z**2>4. else 60)
        n1 = n
        integrand = 0.
        for jSum1, th1 in enumerate(np.linspace(0., 2.*np.pi, n1, \
                                endpoint=False)):
            coeff = (2 if jSum1%2==0 else 4)*2*np.pi/n1/3
                    # coefficient for th1 integration
            a = -(x*np.cos(th1)+y*np.sin(th1)+1.)/z
            b = -(x*np.cos(th1)+y*np.sin(th1)-1.)/z

            for jSum in np.arange(1, n+1):
                t = (b+a)/2.+(b-a)*np.cos((2.*jSum-1)*
                                        np.pi/2./n)/2.
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
                integrand += (coeff*eloop.dot(np.array([x,y,z]))*\
                    np.einsum('ij,i,jk', Ploop, eloop, gloop))
        s = -1./(4.*np.pi**2*R**2*n*abs(z))*integrand
        if tier == 1:
            data.append([xyz[0], xyz[1], xyz[2], s[0], s[1], s[2], f])
        else:
            data.append([xyz[0], xyz[1], xyz[2], s[0], s[1], s[2]])

    np.save("data/%s_atoms_s_%s_T%d_R%d_%04d.npy"%(sample, looptype, tier, R, i_file), data)
