import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

sample = 'Al'
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

I_data = np.load("data/%s_intensities.npy"%sample)
xx = np.append(
        -np.logspace(np.log10(2.887/R), np.log10(0.1/R), 20),
        np.logspace(np.log10(0.1/R), np.log10(2.887/R), 20))
yy = np.append(
        -np.logspace(np.log10(2.887/R), np.log10(0.1/R), 20),
        np.logspace(np.log10(0.1/R), np.log10(2.887/R), 20))
zz = np.append(
        -np.logspace(np.log10(2.887/R), np.log10(0.1/R), 20),
        np.logspace(np.log10(0.1/R), np.log10(2.887/R), 20))
grid_x, grid_y, grid_z = np.meshgrid(xx,yy,zz, indexing='ij')
grid_I = interpolate.griddata(
                I_data[:,0:3], I_data[:,3],
                (grid_x, grid_y, grid_z), method='linear'
            )
I_func = interpolate.RegularGridInterpolator(
                (xx,yy,zz), grid_I
            )

pts1 = [[x, 0., 0.] for x in xx]
pts2 = [[0., x, 0.] for x in xx]
pts3 = [[0., 0., x] for x in xx]
pts4 = [[0., x, x] for x in xx/np.sqrt(2)]
exponent = 1.
plt.plot(xx, xx**exponent*I_func(pts1), label='100') 
plt.plot(xx, xx**exponent*I_func(pts2), label='010') 
plt.plot(xx, xx**exponent*I_func(pts3), label='001') 
plt.plot(xx, xx**exponent*I_func(pts4), label='011') 

plt.legend(loc='upper left')
plt.savefig("plots/Al_intensities.pdf")
plt.show()
