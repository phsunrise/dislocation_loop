import numpy as np
import sys, os
import matplotlib.pyplot as plt
from info import *
from atoms_amplitude import *
from formfactor import formfactor

R = 10.
D = 4.*R
print "sample:", sample
print "R=%.1f, D=%.1f*R"%(R, D/R)
datadir = '%s_R%d/'%(sample, R)

fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(1, 1, 1)

for looptype in looptypes:
    print "starting looptype %s..." % looptype
    intensities = [] 
    #intensities1 = []
    #intensities2 = []
    for i_ori, ori in enumerate(orientations):
        print "starting orientation %d..." % i_ori
        amplitudes = [] 
        ff_tier2 = []
        ff_tier3 = []
        for i_qval, qval in enumerate(q_array):
            q = qval * np.einsum('ij,j', ori, eq)
            qloop = np.einsum('ij,j', rot, q)
            K = np.einsum('ij,j', ori, h) + q
            Kloop = np.einsum('ij,j', rot, K)
            
            ## get the form factor
            ff_tier2.append(formfactor(qloop, 2))
            ff_tier3.append(formfactor(qloop, 3))

            amplitude = 0.
            '''
            first calculate the laue term
            The dislocation loop is inside the region r <= R
            '''
            for x_p in np.arange(-np.ceil(2.*R/a1), np.ceil(2.*R/a1)):
                for y_p in np.arange(-np.ceil(2.*R/a1), np.ceil(2.*R/a1)):
                    r = x_p*ex_p+y_p*ey_p
                    if np.linalg.norm(r) <= R:
                        if looptype == 'int':
                            amplitude += np.exp(\
                                -0.5*r.dot(r)/D**2)*np.cos(qloop.dot(r))
                        elif looptype == 'vac':
                            amplitude -= np.exp(\
                                -0.5*r.dot(r)/D**2)*np.cos(qloop.dot(r))
            amplitudes.append([qval, amplitude])

        amplitudes = np.array(amplitudes)
        #amplitudes1 = np.copy(amplitudes)  # amplitude inside rho<5R, |z|<10R for comparison
        #amplitudes2 = np.copy(amplitudes)  # amplitude inside rho<10R, |z|<20R for comparison
        ff_tier2 = np.array(ff_tier2)
        ff_tier3 = np.array(ff_tier3)

        amplitudes[:,1] += np.load(datadir+"%s_atoms_amplitude_%s_T1_R%d_ori%d_combined.npy"%(\
                            sample, looptype, R, i_ori))[:,1]*2.
        amplitudes[:,1] += np.load(datadir+"%s_atoms_amplitude_%s_T2_R%d_ori%d_combined.npy"%(\
                            sample, looptype, R, i_ori))[:,1]*ff_tier2*2.
        amplitudes[:,1] += np.load(datadir+"%s_atoms_amplitude_%s_T3_R%d_ori%d_combined.npy"%(\
                            sample, looptype, R, i_ori))[:,1]*ff_tier3*2.

        intensities.append(amplitudes[:,1]**2)
        #intensities1.append(amplitudes1[:,1]**2)
        #intensities2.append(amplitudes2[:,1]**2)

    # average over all orientations
    intensities = np.mean(np.array(intensities), axis=0)
    #intensities1 = np.mean(np.array(intensities1), axis=0)
    #intensities2 = np.mean(np.array(intensities2), axis=0)

    # plot q^4/R^2*I 
    if looptype == 'vac':
        linestyle = '--'
    elif looptype == 'int':
        linestyle = '-'
    ax.plot(q_array, intensities*(q_array)**4/R**2, color='r', ls=linestyle)
    #ax.plot(q_array, intensities1*(q_array)**4/R**2, color='b', ls=linestyle)
    #ax.plot(q_array, intensities2*(q_array)**4/R**2, color='g', ls=linestyle)

ax.set_xlim(-0.3, 0.3)
ax.set_ylim(0., 3.)
plt.show()
