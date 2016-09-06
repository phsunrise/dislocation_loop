'''
This code applies only to fcc crystals with {111} dislocation loops
The script is mostly written in loop coordinates
'''
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from info import *
from atoms_amplitude import *

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(1, 1, 1)

for looptype in looptypes:
    print "starting looptype %s..." % looptype
    intensities = [] 
    intensities1 = []
    intensities2 = []
    for i_ori, ori in enumerate(orientations):
        print "starting orientation %d..." % i_ori
        amplitudes = [] 
        ## first calculate the laue term
        for i_qR, qR in enumerate(qR_array):
            q = qR/R * np.einsum('ij,j', ori, eq)
            qloop = np.einsum('ij,j', rot, q)
            K = np.einsum('ij,j', ori, h) + q
            Kloop = np.einsum('ij,j', rot, K)

            amplitude = 0.
            amplitude1 = 0.
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
            amplitudes.append([qR/R, amplitude])

        amplitudes = np.array(amplitudes)
        amplitudes1 = np.copy(amplitudes)  # amplitude inside rho<5R, |z|<10R for comparison
        amplitudes2 = np.copy(amplitudes)  # amplitude inside rho<10R, |z|<20R for comparison

        i_file = 0 
        while os.path.isfile("data/%s_atoms_s_%s_R%d_%04d.npy"%(sample,\
                                    looptype, R, i_file)):
            try:
                amplitudes[:,1] += np.load("data/%s_atoms_amplitude_%s_R%d_ori%d_%04d.npy"%(\
                                    sample, looptype, R, i_ori, i_file))[:,1]
                amplitudes1[:,1] += np.load("data/%s_atoms_amplitude1_%s_R%d_ori%d_%04d.npy"%(\
                                    sample, looptype, R, i_ori, i_file))[:,1]
                amplitudes2[:,1] += np.load("data/%s_atoms_amplitude2_%s_R%d_ori%d_%04d.npy"%(\
                                    sample, looptype, R, i_ori, i_file))[:,1]
                i_file += 1
            except IOError:
                print "file data/%s_atoms_amplitude_%s_R%d_ori%d_%04d.npy does not exist!"%(\
                                    sample, looptype, R, i_ori, i_file)
                sys.exit(0)
        intensities.append(amplitudes[:,1]**2)
        intensities1.append(amplitudes1[:,1]**2)
        intensities2.append(amplitudes2[:,1]**2)
        print "read %d files, done orientation %d" % (i_file, i_ori)

    # average over all orientations
    intensities = np.mean(np.array(intensities), axis=0)
    intensities1 = np.mean(np.array(intensities1), axis=0)
    intensities2 = np.mean(np.array(intensities2), axis=0)

    # plot q^4/R^2*I 
    if looptype == 'vac':
        linestyle = '--'
    elif looptype == 'int':
        linestyle = '-'
    ax.plot(qR_array, intensities*(qR_array/R)**4/R**2, color='r', ls=linestyle)
    ax.plot(qR_array, intensities1*(qR_array/R)**4/R**2, color='b', ls=linestyle)
    ax.plot(qR_array, intensities2*(qR_array/R)**4/R**2, color='g', ls=linestyle)

plt.show()
