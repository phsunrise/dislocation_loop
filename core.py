import numpy as np
import sys, os
from settings import basedir

sample = 'W'
from W_parameters import *
looptypes = ['vac', 'int']

print "sample:", sample
print "R=%.1f, D=%.1f*R"%(R, D/R)
datadir = basedir+"%s_R%d/"%(sample, R)
q_array = np.load(datadir+"q_array.npy")

for looptype in looptypes:
    print "starting looptype %s..." % looptype
    for i_ori, ori in enumerate(orientations):
        print "starting orientation %d..." % i_ori
        amplitudes = [] 
        for i_qval, qval in enumerate(q_array):
            q = qval * np.einsum('ij,j', ori, eq)
            qloop = np.einsum('ij,j', rot, q)
            K = np.einsum('ij,j', ori, h) + q
            Kloop = np.einsum('ij,j', rot, K)
            
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
                            amplitude += np.exp(\
                                -0.5*r.dot(r)/D**2)*np.cos(Kloop.dot(r))
                        elif looptype == 'vac':
                            amplitude -= np.exp(\
                                -0.5*r.dot(r)/D**2)*np.cos(Kloop.dot(r))
            amplitudes.append([qval, amplitude])

        amplitudes = np.array(amplitudes)
        np.save(datadir+"%s_atoms_core_%s_R%d_ori%d.npy"%(\
                    sample, looptype, R, i_ori), amplitudes)
