'''
This code applies only to fcc crystals with {111} dislocation loops
The script is mostly written in loop coordinates
'''
import numpy as np
import sys, os
import matplotlib.pyplot as plt

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

amplitudes = [] 
amplitudes1 = []  # amplitude inside rho<5R, |z|<10R for comparison
for qR in np.linspace(-5., 5., 51):
    q = qR/R * eq
    qloop = np.einsum('ij,j', rot, q)
    K = h+q
    Kloop = np.einsum('ij,j', rot, K)

    amplitude = 0.
    amplitude1 = 0.
    ## first calculate the laue term
    '''
    The dislocation loop is inside the region r <= R, which
        in p coordinates is x^2+y^2+xy <= (R/a1)^2. The constraint on
        x is 3/4*x^2 <= (R/a1)^2
    '''
    _xlim = np.floor(R / a1 * 2./np.sqrt(3))
    for x in np.arange(-_xlim, _xlim+1):
        _ylims = np.roots([1., x, x**2-(R/a1)**2])
        for y in np.arange(np.ceil(np.min(_ylims)), np.floor(np.max(_ylims))+1):
            amplitude += np.cos(qloop.dot(x*ex_p+y*ey_p))
    amplitude1 = amplitude

    ## then calculate the other atoms 
    for x, y, z, sx, sy, sz in sdata:
        qr = qloop.dot([x, y, z])
        Ks = Kloop.dot([-sx, -sy, -sz]) # for interstitial loops
        amplitude += np.cos(qr)*(np.cos(Ks)-1.)-np.sin(qr)*np.sin(Ks)
        if np.sqrt(x**2+y**2)<=5.*R and abs(z)<=10.*R:
            amplitude1 += np.cos(qr)*(np.cos(Ks)-1.)-np.sin(qr)*np.sin(Ks)

    print "qR = %f, amplitude = %f" % (qR, amplitude)
    amplitudes.append([qR, amplitude])
    amplitudes1.append([qR, amplitude1])

# plot and save data
amplitudes = np.array(amplitudes)
amplitudes1 = np.array(amplitudes1)
plt.plot(amplitudes[:,0], amplitudes[:,0]**4*amplitudes[:,1]**2,\
                label='rho<10R, |z|<20R')
plt.plot(amplitudes1[:,0], amplitudes1[:,0]**4*amplitudes1[:,1]**2,\
                label='rho<5R, |z|<10R')
plt.legend()
plt.show()

np.save("data/%s_atoms_amplitudes.npy"%sample, amplitudes)
np.save("data/%s_atoms_amplitudes_1.npy"%sample, amplitudes1)
