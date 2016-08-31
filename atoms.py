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
for qR in np.linspace(-5., 5., 51):
    q = qR/R * eq
    qloop = np.einsum('ij,j', rot, q)
    K = h+q
    Kloop = np.einsum('ij,j', rot, K)

    amplitude = 0.
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

    for x, y, z, sx, sy, sz in sdata:
        qr = qloop.dot([x, y, z])
        Ks = Kloop.dot([sx, sy, sz])
        amplitude += np.cos(qr)*(np.cos(Ks)-1.)-np.sin(qr)*np.sin(Ks)

    print "qR = %f, amplitude = %f" % (qR, amplitude)
    amplitudes.append([qR, amplitude])

# plot and save data
amplitudes = np.array(amplitudes)
plt.plot(amplitudes[:,0], amplitudes[:,0]**4*amplitudes[:,1]**2)
plt.show()

np.save("data/%s_amplitudes.npy"%sample, amplitudes)
