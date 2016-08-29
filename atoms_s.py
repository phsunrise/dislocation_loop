'''
This code applies only to fcc crystals with {111} dislocation loops
The script is mostly written in loop coordinates
'''
import numpy as np
import sys, os

#sample = sys.argv[1]
sample = 'Al'
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

P = np.zeros((3, 3))
for i in xrange(3):
    for j in xrange(3):
        P[i,j] = ((C12*F.dot(B) + d*F[i]*B[i]) * \
                    (1 if i==j else 0) + \
                    C44*(F[i]*B[j]+F[j]*B[i])
            )
Ploop =  np.einsum('ij,kl,jl', rot, rot, P)


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

rhomax = 5.*R
zmax = 10.*R
## the summation goes to rho = 5R in xy plane, and z from -10R to 10R
##     (both roughly, since the origins will give some offset)
for zind in np.arange(np.ceil(-zmax/(3.*a2)), np.floor(zmax/(3.*a2))):
    _xlim = np.floor(rhomax / a1 * 2./np.sqrt(3))
    for xind in np.arange(-_xlim, _xlim+1):
        _ylims = np.roots([1., x, x**2-(rhomax/a1)**2])
        for yind in np.arange(np.ceil(np.min(_ylims)), np.floor(np.max(_ylims))+1):
