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
    if np.linalg.norm(line[0:3]) < 3.*R: 
        sdata_inside.append(line)

sdata_inside = np.array(sdata_inside).T
ax.quiver(sdata_inside[0], sdata_inside[1], sdata_inside[2],\
          sdata_inside[3], sdata_inside[4], sdata_inside[5])
plt.show()
