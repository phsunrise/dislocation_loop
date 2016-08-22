import numpy as np
import os, sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import interpolate

sample = 'Al'
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

I_data = np.load("data/%s_intensities.npy"%sample)
I_datax = []
for i in xrange(len(I_data)):
    if abs(I_data[i, 1]) < 1.e-6:
        I_datax.append(I_data[i])
I_datax = np.array(I_datax)

xi = np.append(
        -np.logspace(np.log10(5./R), np.log10(0.1/R), 20),
        np.logspace(np.log10(0.1/R), np.log10(5./R), 20))
zi = np.append(
        -np.logspace(np.log10(5./R), np.log10(0.1/R), 20),
        np.logspace(np.log10(0.1/R), np.log10(5./R), 20))
xi = np.linspace(-5./R, 5./R, 50)
zi = np.linspace(-5./R, 5./R, 50)
Ii = interpolate.griddata((I_datax[:,0], I_datax[:,2]), I_datax[:,3],
                        (xi[None,:], zi[:,None]), method='linear')
Ii = np.log10(Ii)

fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(1,1,1)
CS = ax.contour(xi, zi, Ii, 15, colors='k')
CS = ax.contourf(xi, zi, Ii, 15, cmap=plt.cm.jet)
fig.colorbar(CS, fraction=0.05, pad=0.025)

#ax.scatter(I_datax[:,0], I_datax[:,2], marker='o', c='b', s=5)
ax.set_xlim(-5./R, 5./R)
ax.set_ylim(-5./R, 5./R)
plt.show()
