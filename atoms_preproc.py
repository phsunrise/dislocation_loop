'''
This script generates the positions of atoms that we wish to calculate
It is a list of (x,y,z) values
'''
import numpy as np
import os, sys

sample = sys.argv[1] 
looptype = 'int'
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *
elif sample == 'W':
    from W_parameters import *

## first check if the preproc folder exists
if not os.path.isdir("preproc/"):
    print "Creating folder preproc/ ..."
    os.system("mkdir preproc")

## get existing files
try:
    xyz_list_orig = np.load("preproc/%s_atoms_s_%s_pre_0000.npy"%(\
                                sample, looptype))
    i_file = 1
    while os.path.isfile("preproc/%s_atoms_s_%s_pre_%04d.npy"%(\
                                sample, looptype, i_file)):
        xyz_list_orig = np.vstack((xyz_list_orig, \
                  np.load("preproc/%s_atoms_s_%s_pre_%04d.npy"%(\
                                sample, looptype, i_file))))
        i_file += 1
    nfiles = i_file
except IOError:
    xyz_list_orig = []
    nfiles = 0

## the summation goes to rhomax in xy plane, and z from -zmax to zmax 
##     (both roughly, since the origins will give some offset)
rhomax = 300.
zmax = 300.

xyz_list = []
for z_p in np.arange(np.ceil(-zmax/a2), np.floor(zmax/a2)+1):
    if looptype == 'vac' and z_p == 0.:
        continue
    if crystaltype == 'FCC':
        if z_p % 3 == 0.: 
            x0, y0 = orig_C[0:2]
        elif z_p % 3 == 1.:
            x0, y0 = orig_A[0:2]
        elif z_p % 3 == 2.:
            x0, y0 = orig_B[0:2]
    elif crystaltype == 'BCC':
        if z_p % 2 == 0.:
            x0, y0 = orig_A[0:2]
        elif z_p % 2 == 1.:
            x0, y0 = orig_B[0:2]

    _xlim = np.floor(rhomax / a1 * 2./np.sqrt(3))
    for x_p in np.arange(-_xlim, _xlim+1):
        _ylims = np.roots([1., x_p, x_p**2-(rhomax/a1)**2])
        for y_p in np.arange(np.ceil(np.min(_ylims)), \
                                 np.floor(np.max(_ylims))+1):
            if looptype == 'vac':
                xyz = (x_p+x0)*ex_p + (y_p+y0)*ey_p + z_p*ez_p
            elif looptype == 'int':
                xyz = (x_p+x0)*ex_p + (y_p+y0)*ey_p + (z_p+0.5)*ez_p

            xyz_list.append(xyz)

#xyz_list_orig = set(tuple(xyz) for xyz in xyz_list_orig)
#xyz_list = set(tuple(xyz) for xyz in xyz_list)
#xyz_list = np.array(list(xyz_list - xyz_list_orig))

nproc = raw_input("Total length %d, number of processors? "%len(xyz_list))
nproc = int(nproc)

sublistlen = int(np.ceil(len(xyz_list)*1./nproc))
for i in xrange(nproc):
    np.save("preproc/%s_atoms_s_%s_pre_%04d.npy"%(\
        sample, looptype, i+nfiles), \
        np.array(xyz_list[i*sublistlen:(i+1)*sublistlen]))
    print "file %d saved"%i
