'''
This script generates the positions of atoms that we wish to calculate
It is a list of (x/R,y/R,z/R) values
'''
import numpy as np
import os, sys

sample = sys.argv[1] 
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

## first check if the preproc folder exists
if not os.path.isdir("preproc/"):
    print "Creating folder preproc/ ..."
    os.system("mkdir preproc")

# remove existing files
os.system("rm preproc/%s_atoms_s_pre_????.npy"%sample)

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
ex_p = a1*np.array([np.cos(ex_p_th), np.sin(ex_p_th), 0.])
ey_p = a1*np.array([np.cos(ey_p_th), np.sin(ey_p_th), 0.])
ez_p = a2*np.array([0., 0., 1.])
## From now on, the coordinates will be in units of ex_p, ey_p, ez_p
'''
Here we assume that the dislocation loop has the configuration of plane A
In this case, starting from the origin, the planes will be (A), C, A, B, C, A, ... 
    with z coordinates (0), 0.5, 1.5, 2.5, 3.5, 4.5, ... in units of a2
The following is the coordinates for one atom on each plane 
'''
orig_A = np.array([0., 0., 1.5])
orig_B = np.array([-1./3, -1./3, 2.5])
orig_C = np.array([1./3, 1./3, 0.5])

## the summation goes to rho = 5R in xy plane, and z from -10R to 10R
##     (both roughly, since the origins will give some offset)
rhomax = 5.*R
zmax = 10.*R

xyz_list = []
for z_p in 0.5+np.arange(np.ceil(-zmax/a2), np.floor(zmax/a2)+1):
    if z_p % 3 == 0.5: 
        x0, y0 = orig_C[0:2]
    elif z_p % 3 == 1.5:
        x0, y0 = orig_A[0:2]
    elif z_p % 3 == 2.5:
        x0, y0 = orig_B[0:2]

    _xlim = np.floor(rhomax / a1 * 2./np.sqrt(3))
    for x_p in np.arange(-_xlim, _xlim+1):
        _ylims = np.roots([1., x_p, x_p**2-(rhomax/a1)**2])
        for y_p in np.arange(np.ceil(np.min(_ylims)), \
                                 np.floor(np.max(_ylims))+1):
            xyz_list.append(((x_p+x0)*ex_p + (y_p+y0)*ey_p + z_p*ez_p)/R)

nproc = raw_input("Total length %d, number of processors? "%len(xyz_list))
nproc = int(nproc)

sublistlen = int(np.ceil(len(xyz_list)*1./nproc))
for i in xrange(nproc):
    np.save("preproc/%s_atoms_s_pre_%04d.npy"%(sample, i), \
        np.array(xyz_list[i*sublistlen:(i+1)*sublistlen]))
    print "file %d saved"%i
