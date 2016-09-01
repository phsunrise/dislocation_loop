'''
This code applies only to fcc crystals with {111} dislocation loops
The script is mostly written in loop coordinates
'''
import numpy as np
import sys, os
import matplotlib.pyplot as plt

sample = 'Cu'
looptype = 'int'
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

## read s data
sdata = np.load("data/%s_atoms_s_%s_R%d_0000.npy"%(sample, looptype, R))
i_file = 1
while os.path.isfile("data/%s_atoms_s_%s_R%d_%04d.npy"%(sample,\
                            looptype, R, i_file)):
    sdata = np.vstack((sdata, \
                np.load("data/%s_atoms_s_%s_R%d_%04d.npy"%(sample, \
                            looptype, R, i_file))))
    i_file += 1
print "read %d files" % i_file

amplitudes = [] 
amplitudes1 = []  # amplitude inside rho<5R, |z|<10R for comparison
for qR in np.linspace(-5., 5., 21):
    q = qR/R * eq
    qloop = np.einsum('ij,j', rot, q)
    K = h+q
    Kloop = np.einsum('ij,j', rot, K)

    amplitude = 0.
    amplitude1 = 0.
    ## first calculate the laue term
    if looptype == 'int':
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
    elif looptype == 'vac':
        rhocut = 10.*R
        _xlim = np.floor(rhocut / a1 * 2./np.sqrt(3))
        for x in np.arange(-_xlim, _xlim+1):
            _ylims = np.roots([1., x, x**2-(rhocut/a1)**2])
            for y in np.arange(np.ceil(np.min(_ylims)), np.floor(np.max(_ylims))+1):
                ''' According to the configuration defined in the 
                        parameter file, the z=0 plane should be type C
                '''
                if np.linalg.norm((x+orig_C[0])*ex_p+(y+orig_C[1])*ey_p) > R:
                    amplitude += np.cos(qloop.dot(
                            (x+orig_C[0])*ex_p+(y+orig_C[1])*ey_p))
                    if np.linalg.norm((x+orig_C[0])*ex_p+(y+orig_C[1])*ey_p) < 5.*R:
                        amplitude1 += np.cos(qloop.dot(
                                (x+orig_C[0])*ex_p+(y+orig_C[1])*ey_p))

    ## then calculate the other atoms 
    for x, y, z, sx, sy, sz in sdata:
        qr = qloop.dot([x, y, z])
        if looptype == 'vac':
            Ks = Kloop.dot([sx, sy, sz])
        elif looptype == 'int':
            Ks = Kloop.dot([-sx, -sy, -sz])

        amplitude += np.cos(qr)*(np.cos(Ks)-1.)-np.sin(qr)*np.sin(Ks)
        if np.sqrt(x**2+y**2)<=5.*R and abs(z)<=10.*R:
            amplitude1 += np.cos(qr)*(np.cos(Ks)-1.)-np.sin(qr)*np.sin(Ks)

    print "q = %f, amplitude = %f" % (qR/R, amplitude)
    amplitudes.append([qR/R, amplitude])
    amplitudes1.append([qR/R, amplitude1])

# plot and save data
amplitudes = np.array(amplitudes)
amplitudes1 = np.array(amplitudes1)
plt.plot(amplitudes[:,0], amplitudes[:,0]**4/R**2*amplitudes[:,1]**2,\
                label='rho<10R, |z|<20R')
plt.plot(amplitudes1[:,0], amplitudes1[:,0]**4/R**2*amplitudes1[:,1]**2,\
                label='rho<5R, |z|<10R')
plt.legend()
plt.show()

#np.save("data/%s_atoms_%s_R%d_amplitudes.npy"%(\
#                    sample, looptype, R), amplitudes)
#np.save("data/%s_atoms_%s_R%d_amplitudes1.npy"%(\
#                    sample, looptype, R), amplitudes1)
