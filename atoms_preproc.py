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
