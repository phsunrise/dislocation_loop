import sys, os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
sdata_inside = [] 
for line in sdata:
    if np.linalg.norm(line[0:3]) < 2.*R and line[2]>=-a0/3.: 
        sdata_inside.append(line)

sdata_inside = np.array(sdata_inside).T
#ax.quiver(sdata_inside[0], sdata_inside[1], sdata_inside[2],\
#          sdata_inside[3], sdata_inside[4], sdata_inside[5])
ax.scatter(sdata_inside[0]-sdata_inside[3], \
           sdata_inside[1]-sdata_inside[4], \
           sdata_inside[2]-sdata_inside[5], \
           c='b')

## atoms inside loop
loopatoms = []
_xlim = np.floor(R / a1 * 2./np.sqrt(3))
for x in np.arange(-_xlim, _xlim+1):
    _ylims = np.roots([1., x, x**2-(R/a1)**2])
    for y in np.arange(np.ceil(np.min(_ylims)), np.floor(np.max(_ylims))+1):
        loopatoms.append(x*ex_p+y*ey_p)
loopatoms = np.array(loopatoms).T
ax.scatter(loopatoms[0], loopatoms[1], 0., c='r')
plt.show()
