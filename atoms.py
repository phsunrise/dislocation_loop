'''
This code applies only to fcc crystals with {111} dislocation loops
The script is mostly written in loop coordinates
'''
import numpy as np
import sys, os
import matplotlib.pyplot as plt

sample = 'Al'
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

## read s data
sdata = np.load("data/%s_atoms_s_0.npy"%sample)
i_file = 1
while os.path.isfile("data/%s_atoms_s_%d.npy"%(sample, i_file)):
    sdata = np.vstack((sdata, \
                np.load("data/%s_atoms_s_%d.npy"%(sample, i_file))))
    i_file += 1
print "read %d files" % i_file

'''
In loop coordinates, there should be 3 different planes,
    each of which is a triangular lattice. The ex_p and ey_p vectors
    indicate the orientation of the two sides of the triangle
The lattice orientation is the same for all 3 planes; they only
    differ by a constant translation
'''
## side length of triangle
a1 = a0*np.sqrt(2)/2.
## distance between planes
a2 = a0*np.sqrt(3)/3.
ex_p_th = np.pi/6.
ey_p_th = np.pi/2.
ex_p = a1*(ex*np.cos(ex_p_th)+ey*np.sin(ex_p_th))
ey_p = a1*(ex*np.cos(ey_p_th)+ey*np.sin(ey_p_th))
ez_p = a2*ez
## From now on, the coordinates will be in units of ex_p, ey_p, ez_p
'''
Here we assume that the dislocation loop has the configuration of plane A
In this case, starting from the origin, the planes will be (A), C, A, B, C, A, ... 
    with z coordinates (0), 0.5, 1.5, 2.5, 3.5, 4.5, ... in units of a2
The following is the coordinates for one atom on each plane 
'''
orig_A = np.array([0., 0., 1.5])
orig_B = np.array([1./3, 1./3, 2.5])
orig_C = np.array([-1./3, 2./3, 0.5])

## the summation goes to rho = 5R in xy plane, and z from -10R to 10R
amplitudes = [] 
for qR in np.linspace(-5., 5., 51):
    if qR == 0:
        continue
    q = qR/R * eq
    qloop = np.einsum('ij,j', rot, q)
    K = h+q
    Kloop = np.einsum('ij,j', rot, K)

    amplitude = 0.
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
            amplitude += np.cos(qloop.dot(x*ex_p+y*ey_p))

    for x, y, z, sx, sy, sz in sdata:
        qr = qloop.dot([x, y, z])
        Ks = Kloop.dot([sx, sy, sz])
        amplitude += np.cos(qr)*(np.cos(Ks)-1.)-np.sin(qr)*np.sin(Ks)

    print "qR = %f, amplitude = %f" % (qR, amplitude)
    amplitudes.append([qR, amplitude])

# plot and save data
amplitudes = np.array(amplitudes)
plt.plot(amplitudes[:,0], amplitudes[:,1]**2)
plt.show()

np.save("data/%s_amplitudes.npy"%sample, amplitudes)
