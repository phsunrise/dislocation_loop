import numpy as np
import sys, os
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
            if not os.path.isfile(datadir+"%s_atoms_s_%s_T%d_R%d_%04d.npy"%(\
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

        s_prev = []
        n = 1 
        n1 = 1 
        th1_array = [2.*np.pi]
        th_array = [np.pi]
        integrand = 0.
        flag = 0 # 0 to double n1 first, 1 to double n first
        justchangedflag = False
        while True:
            for th1 in th1_array: 
                a = -(x*np.cos(th1)+y*np.sin(th1)+1.)/z
                b = -(x*np.cos(th1)+y*np.sin(th1)-1.)/z
                for th in th_array: 
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
                    integrand += (eloop.dot(np.array([x,y,z]))*\
                        np.einsum('ij,i,jk', Ploop, eloop, gloop))
            s = -1./(4.*np.pi**2*R**2*n*abs(z))*integrand*2.*np.pi/n1

            if s_prev != [] and max(abs((s-s_prev)/s)) < 1.e-10:
                if justchangedflag:
                    break
                s_prev = s
                flag = (flag+1) % 2
                justchangedflag = True
                if flag == 0:
                    n1 *= 2 
                    th1_array = np.arange(2*np.pi/n1, 2.*np.pi*(1.+1./n1), 4.*np.pi/n1)
                    th_array = np.arange(np.pi/n, np.pi*(1.+1./n), np.pi/n)
                elif flag == 1:
                    n *= 2
                    th1_array = np.arange(2*np.pi/n1, 2.*np.pi*(1.+1./n1), 2.*np.pi/n1)
                    th_array = np.arange(np.pi/n, np.pi*(1.+1./n), 2.*np.pi/n)
                continue

            justchangedflag = False
            s_prev = s
            if flag == 0:
                n1 *= 2 
                th1_array = np.arange(2*np.pi/n1, 2.*np.pi*(1.+1./n1), 4.*np.pi/n1)
            elif flag == 1:
                n *= 2
                th_array = np.arange(np.pi/n, np.pi*(1.+1./n), 2.*np.pi/n)

        if tier == 1:
            data.append([xyz[0], xyz[1], xyz[2], s[0], s[1], s[2], xyz[3]])
        else:
            data.append([xyz[0], xyz[1], xyz[2], s[0], s[1], s[2]])

    data = np.array(data)
    if looptype = 'int':
        data[:,3:6] = -data[:,3:6]
    np.save(datadir+"%s_atoms_s_%s_T%d_R%d_%04d.npy"%(sample, looptype, tier, R, i_file), data)
