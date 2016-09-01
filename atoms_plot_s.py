import sys, os
if sys.platform == 'darwin':
    import matplotlib
    matplotlib.use("Qt4Agg")
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

sample = 'Cu'
looptype = 'vac'

if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

## read s data
sdata = np.load("data/%s_atoms_s_%s_R%d_0000.npy"%(sample, looptype, R))
i_file = 1
while os.path.isfile("data/%s_atoms_s_%s_R%d_%04d.npy"%(sample, looptype, R, i_file)):
    sdata = np.vstack((sdata, \
                np.load("data/%s_atoms_s_%s_R%d_%04d.npy"%(sample, looptype, R, i_file))))
    i_file += 1
print "read %d files" % i_file

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
sdata_inside1 = [] 
sdata_inside2 = [] 
for line in sdata:
    if np.linalg.norm(line[0:3]) < 2.*R and line[2]>=0.: 
        sdata_inside1.append(line)
    elif np.linalg.norm(line[0:3]) < 2.*R and line[2]>=-a0: 
        sdata_inside2.append(line)

sdata_inside1 = np.array(sdata_inside1).T
sdata_inside2 = np.array(sdata_inside2).T
if looptype == 'int':
    sdata_inside1[3:6] = -1.*sdata_inside1[3:6]
    sdata_inside2[3:6] = -1.*sdata_inside2[3:6]
elif looptype == 'vac':
    sdata_inside1[3:6] = 1.*sdata_inside1[3:6]
    sdata_inside2[3:6] = 1.*sdata_inside2[3:6]
#ax.quiver(sdata_inside[0], sdata_inside[1], sdata_inside[2],\
#          sdata_inside[3], sdata_inside[4], sdata_inside[5])
ax.scatter(sdata_inside1[0]+sdata_inside1[3], \
           sdata_inside1[1]+sdata_inside1[4], \
           sdata_inside1[2]+sdata_inside1[5], \
           c='b')
ax.scatter(sdata_inside2[0]+sdata_inside2[3], \
           sdata_inside2[1]+sdata_inside2[4], \
           sdata_inside2[2]+sdata_inside2[5], \
           c='g')

if looptype == 'int':
    ## atoms inside loop
    loopatoms = []
    _xlim = np.floor(R / a1 * 2./np.sqrt(3))
    for x in np.arange(-_xlim, _xlim+1):
        _ylims = np.roots([1., x, x**2-(R/a1)**2])
        for y in np.arange(np.ceil(np.min(_ylims)), np.floor(np.max(_ylims))+1):
            loopatoms.append(x*ex_p+y*ey_p)
    loopatoms = np.array(loopatoms).T
    ax.scatter(loopatoms[0], loopatoms[1], 0., c='r')
elif looptype == 'vac':
    ### atoms outside of loop
    #loopatoms = []
    #_xlim = np.ceil(R / a1 * 2./np.sqrt(3))
    #for x in np.arange(-_xlim, _xlim+1):
    #    _ylims = np.roots([1., x, x**2-(R/a1)**2])
    #    for y in np.arange(np.ceil(np.min(_ylims)), np.floor(np.max(_ylims))+1):
    #        loopatoms.append(x*ex_p+y*ey_p)
    pass

plt.show()
