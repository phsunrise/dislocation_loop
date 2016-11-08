import numpy as np
import sys, os
from settings import basedir
import glob

sample = 'W'
from W_parameters import R 
looptypes = ['vac', 'int']

print "sample:", sample
i_dir = 0
while os.path.isdir(basedir + "%s_R%d_amp2d_%d/"%(sample, R, i_dir)):
    print "folder %s exists, containing %d files" % (\
            basedir+"%s_R%d_amp2d_%d/"%(sample, R, i_dir), \
            len(glob.glob(basedir+"%s_R%d_amp2d_%d/%s_atoms_amplitude_2d_*_T*_R*_ori*_????.npy"%(\
                sample, R, i_dir, sample))))
    i_dir += 1
val = raw_input("enter folder number: ")
i_dir = int(val)
datadir = basedir+"%s_R%d_amp2d_%d/"%(sample, R, i_dir)

sys.path.insert(0, datadir)
from W_parameters import * 
print "R=%.1f, D=%.1f*R"%(R, D/R)

_file = np.load(datadir+"qarray_2d.npz")
q_array = _file['q_array'] 
K_array = _file['K_array']
h = _file['h']
print "h = 2pi/a0 *", h*a0/2./np.pi 

for looptype in looptypes:
    print "starting looptype %s..." % looptype
    for i_ori, ori in enumerate(orientations):
        print "starting orientation %d..." % i_ori
        amplitudes = np.zeros((q_array.shape[0], q_array.shape[1]))
        counter = 0
        for index, _ in np.ndenumerate(amplitudes):
            q = q_array[index]
            qloop = np.einsum('ij,jk,k', ori, rot, q)
            K = K_array[index] 
            Kloop = np.einsum('ij,jk,k', ori, rot, K)
            
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
            amplitudes[index] = amplitude
            counter += 1
            if counter % 200 == 0:
                print "done %.2f%%\r" % (counter*100./np.prod(amplitudes.shape)),
                sys.stdout.flush()

        print "done 100.00%%"
        np.save(datadir+"%s_atoms_core_2d_%s_R%d_ori%d.npy"%(\
                    sample, looptype, R, i_ori), amplitudes)
