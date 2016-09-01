'''
This code applies only to fcc crystals with {111} dislocation loops
The script is mostly written in loop coordinates
'''
import numpy as np
import sys, os
import matplotlib.pyplot as plt

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

sample = 'Cu'
looptypes = ['int', 'vac']
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

D = 2.5*R ## parameter in gaussian used to smooth out the boundary
_list = [(looptype, i_ori, ori) for looptype in looptypes for i_ori, ori in enumerate(orientations)]
for looptype, i_ori, ori in _list:
    ## read s data
    sdata = np.load("data/%s_atoms_s_%s_R%d_0000.npy"%(sample, looptype, R))
    i_file = 1
    while os.path.isfile("data/%s_atoms_s_%s_R%d_%04d.npy"%(sample,\
                                looptype, R, i_file)):
        sdata = np.vstack((sdata, \
                    np.load("data/%s_atoms_s_%s_R%d_%04d.npy"%(sample, \
                                looptype, R, i_file))))
        i_file += 1
    print "read %d files" % i_file

    amplitudes = [] 
    amplitudes1 = []  # amplitude inside rho<5R, |z|<10R for comparison
    for i_qR, qR in enumerate(np.linspace(-5., 5., 101)):
        ## here is the parallelism criterion
        if i_qR % nprocs != rank:
            continue
        q = qR/R * np.einsum('ij,j', ori, eq)
        qloop = np.einsum('ij,j', rot, q)
        K = np.einsum('ij,j', ori, h) + q
        Kloop = np.einsum('ij,j', rot, K)

        amplitude = 0.
        amplitude1 = 0.
        ## first calculate the laue term
        '''
        The dislocation loop is inside the region r <= R, which
            in p coordinates is x^2+y^2+xy <= (R/a1)^2. The constraint on
            x is 3/4*x^2 <= (R/a1)^2
        '''
        _xlim = np.floor(R / a1 * 2./np.sqrt(3))
        for x in np.arange(-_xlim, _xlim+1):
            _ylims = np.roots([1., x, x**2-(R/a1)**2])
            for y in np.arange(np.ceil(np.min(_ylims)), np.floor(np.max(_ylims))+1):
                if looptype == 'int':
                    r = x*ex_p+y*ey_p
                    amplitude += np.exp(-0.5*r.dot(r)/D**2)*np.cos(qloop.dot(r))
                elif looptype == 'vac':
                    ''' According to the configuration defined in the 
                            parameter file, the z=0 plane should be type C
                    '''
                    r = (x+orig_C[0])*ex_p+(y+orig_C[1])*ey_p
                    amplitude -= np.exp(-0.5*r.dot(r)/D**2)*np.cos(qloop.dot(r))
        amplitude1 = amplitude

        ## then calculate the other atoms 
        for x, y, z, sx, sy, sz in sdata:
            r = np.array([x,y,z])
            qr = qloop.dot(r)
            if looptype == 'vac':
                Ks = Kloop.dot([sx, sy, sz])
            elif looptype == 'int':
                Ks = Kloop.dot([-sx, -sy, -sz])

            amplitude += (np.cos(qr)*(np.cos(Ks)-1.)-\
                    np.sin(qr)*np.sin(Ks))*np.exp(-0.5*r.dot(r)/D**2)
            if np.sqrt(x**2+y**2)<=5.*R and abs(z)<=10.*R:
                amplitude1 += (np.cos(qr)*(np.cos(Ks)-1.)-\
                    np.sin(qr)*np.sin(Ks))*np.exp(-0.5*r.dot(r)/D**2)

        print "q = %f, amplitude = %f" % (qR/R, amplitude)
        amplitudes.append([qR/R, amplitude])
        amplitudes1.append([qR/R, amplitude1])

    # plot and save data
    amplitudes = np.array(amplitudes)
    amplitudes1 = np.array(amplitudes1)
    #plt.plot(amplitudes[:,0], amplitudes[:,0]**4/R**2*amplitudes[:,1]**2,\
    #                label='rho<10R, |z|<20R')
    #plt.plot(amplitudes1[:,0], amplitudes1[:,0]**4/R**2*amplitudes1[:,1]**2,\
    #                label='rho<5R, |z|<10R')
    #plt.legend()
    #plt.show()

    np.save("data/%s_atoms_%s_R%d_amplitudes_ori%d_%04d.npy"%(\
                        sample, looptype, R, i_ori, rank), amplitudes)
    np.save("data/%s_atoms_%s_R%d_amplitudes1_ori%d_%04d.npy"%(\
                        sample, looptype, R, i_ori, rank), amplitudes1)
