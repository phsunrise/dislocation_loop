import sys, os
if sys.platform == 'darwin':
    import matplotlib
    matplotlib.use("Qt4Agg")
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

sample = 'W'
from W_parameters import *
looptype = 'vac'

datadir = "%s_R%d/" % (sample, R)

if len(sys.argv) == 3:
    i_filemin, i_filemax = int(sys.argv[1]), int(sys.argv[2])
else:
    i_filemin, i_filemax = 0, np.iinfo(np.int32).max

## range of data to be shown
rhocut = 1.5*R
zcuthigh = 3*a0
zcutlow = -2*a0
xpycuthigh = 0*a0 
xpycutlow = -R

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
sdata_inside1 = [] 
sdata_inside2 = [] 
i_file = i_filemin 
while os.path.isfile(datadir+"%s_atoms_s_%s_T1_R%d_%04d_combined.npy"%(\
        sample, looptype, R, i_file)) and i_file <= i_filemax:
    sdata = np.load(datadir+"%s_atoms_s_%s_T1_R%d_%04d_combined.npy"%(\
            sample, looptype, R, i_file))
    sdata = np.vstack((sdata, -sdata))
    counter = 0
    for line in sdata:
        if not xpycutlow <= line[0]+line[1]+line[3]+line[4] <= xpycuthigh:
            continue
        if np.linalg.norm(line[0:2])<=rhocut and 0.<line[2]<=zcuthigh:
            sdata_inside1.append(line) ## lattice above z=0
            #print np.einsum('ij,j', rot.T, line[0:3]/a0)
            counter += 1
        elif np.linalg.norm(line[0:2])<=rhocut and zcutlow<=line[2]<0.: 
            sdata_inside2.append(line) ## lattice below z=0
            counter += 1
    if counter > 0:
        print "processed file %d, added %d points" % (i_file, counter)
    i_file += 1

sdata_inside1 = np.array(sdata_inside1).T
sdata_inside2 = np.array(sdata_inside2).T
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
    for xpy in np.arange(-np.ceil(5.*R/a1), np.ceil(5.*R/a1)):
        for xmy in np.arange(-np.ceil(5.*R/a1), np.ceil(5.*R/a1)):
            if (xpy+xmy) % 2 == 1:
                x_p = (xpy+xmy)/2.
                y_p = (xpy-xmy)/2.
                r = x_p*ex_p + y_p*ey_p
                if xpycutlow<=r[0]+r[1]<=xpycuthigh and np.linalg.norm(r) <= min(rhocut, R):
                    loopatoms.append(r)
elif looptype == 'vac':
    ## atoms outside of loop
    loopatoms = []
    for xpy in np.arange(-np.ceil(4.*rhocut/a1), np.ceil(4.*rhocut/a1)):
        for xmy in np.arange(-np.ceil(4.*rhocut/a1), np.ceil(4.*rhocut/a1)):
            if (xpy+xmy) % 2 == 0:
                x_p = (xpy+xmy)/2.
                y_p = (xpy-xmy)/2.
                r = x_p*ex_p + y_p*ey_p
                if xpycutlow<=r[0]+r[1]<=xpycuthigh and 0 < np.linalg.norm(r) <= rhocut:
                    loopatoms.append(r)

loopatoms = np.array(loopatoms).T
if len(loopatoms) > 0:
    ax.scatter(loopatoms[0], loopatoms[1], 0., c='r')

## now set axis scale to be the same
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
zmin, zmax = ax.get_zlim()
xmid = (xmin+xmax) / 2.
ymid = (ymin+ymax) / 2.
zmid = (zmin+zmax) / 2.
max_range = max(xmax-xmin, ymax-ymin, zmax-zmin)
ax.set_xlim(xmid-max_range/2., xmid+max_range/2.)
ax.set_ylim(ymid-max_range/2., ymid+max_range/2.)
ax.set_zlim(zmid-max_range/2., zmid+max_range/2.)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.show()
