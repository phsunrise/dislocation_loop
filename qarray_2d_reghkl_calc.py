import numpy as np

sample = 'W'
if sample == 'W':
    from W_parameters import *

E = 10.0 # keV
wavelength = 12.39842 / E # Angstrom
k = 2*np.pi/wavelength

h = 2*np.pi/a0*np.array([2.,2.,0.]) # NOTE: this h follows the notation in Larson's papers; it is NOT the reciprocal lattice h

## we are taking the slice on the (1 -1 0) plane
hk_low, hk_high, N_hk = 1.75, 2.25, 128
l_low, l_high, N_l = -0.0125, 0.0125, 5
hk_array = (np.linspace(hk_low, hk_high, N_hk, endpoint=False) \
                + (hk_high-hk_low)/N_hk/2.)
l_array = (np.linspace(l_low, l_high, N_l, endpoint=False) \
                + (l_high-l_low)/N_l/2.)

q_array = np.zeros((len(hk_array), len(l_array), 3))
K_array = np.zeros((len(hk_array), len(l_array), 3))
q_array_norm = np.zeros((len(hk_array), len(l_array)))
for i_hk, hk in enumerate(hk_array):
    for i_l, l in enumerate(l_array):
        K = 2*np.pi/a0*np.array([hk, hk, l])
        K_array[i_hk, i_l] = K 
        q_array[i_hk, i_l] = K-h
        q_array_norm[i_hk, i_l] = np.linalg.norm(K-h)
    print "done %.2f%%" % ((i_hk+1)*1./len(hk_array)*100)

np.savez("qarray_2d.npz", q_array=q_array, K_array=K_array, h=h, \
         hk_array=hk_array, l_array=l_array)

## plot norm of q
import matplotlib.pyplot as plt
plt.imshow(q_array_norm, origin='lower', interpolation='nearest', \
           extent=[l_low, l_high, hk_low, hk_high])
plt.xlabel("l")
plt.ylabel("h or k")
plt.colorbar()
plt.show()
