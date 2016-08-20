'''
Script to calculate s(r) vector, parallelized
To run:
    bsub -W %TIME -o %OUT -e %ERR -q %QUEUE mpiexec -n %NPROC \
        python s_parallel.py 
The output is several files, each containing a n*6 table; 
each line is [r/R, th, z/R, sx, sy, sz]
'''
import os, sys
import numpy as np

## if 2 more arguments are provided, use "fake parallelism"
if len(sys.argv) == 3:
    nprocs = int(sys.argv[1])
    rank = int(sys.argv[2])
else:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nprocs = comm.Get_size()

sample = 'Al'
do_save = True # check the range before changing this to True!!

if sample == 'Cu':
    from Cu_parameters import *
elif sample == 'Al':
    from Al_parameters import *

P = np.zeros((3, 3))
for i in xrange(3):
    for j in xrange(3):
        P[i,j] = ((C12*F.dot(B) + d*F[i]*B[i]) * \
                    (1 if i==j else 0) + \
                    C44*(F[i]*B[j]+F[j]*B[i])
            )
Ploop =  np.einsum('ij,kl,jl', rot, rot, P)

for i_r, r in enumerate(np.linspace(0., 5., 51)):
    if i_r % nprocs != rank:
        continue
    if r == 0.:
        th_list = np.array([0.])
    else:
        th_list = np.linspace(0., 2*np.pi, 20, endpoint=False)

    for th in th_list: 
        x = r*np.cos(th)
        y = r*np.sin(th)
        if do_save:
            try:
                data = np.load("data/%s_s_%d.npy"%(sample, rank))
            except IOError:
                data = []

        for z in np.linspace(0.1, 10., 100):
            n = (30 if abs(r)+abs(z)>2. else 60)
                # used in Chebyshev-Gauss numerical quadrature
            n1 = n # must be even, used in th1 integration
            integrand = 0
            for jSum1, th1 in enumerate(np.linspace(0., 2.*np.pi, n1, endpoint=False)):
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
            s = -1./(4.*np.pi**2*n*abs(z))*integrand/np.linalg.norm(B)
            if do_save:
                try:
                    data = np.vstack((data, [r,th,z,s[0],s[1],s[2]]))
                except ValueError:
                    data = [r,th,z,s[0],s[1],s[2]] 
            print r,th,z,s

        if do_save:
            np.save("data/%s_s_%d.npy"%(sample, rank), data)

#if do_save:
#    with open("%s_s.csv"%sample, 'wb') as f:
#        wr = csv.writer(f, delimiter=' ', quoting=csv.QUOTE_NONNUMERIC)
#        data = np.load("%s_s.npy"%sample)
#        for i_row in xrange(len(data)):
#            wr.writerow(list(data[i_row]))
