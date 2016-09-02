import matplotlib.pyplot as plt
import numpy as np

sample = 'Cu'
looptypes = ['int', 'vac']
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

fig = plt.figure()
ax = fig.add_subplot(111)

for looptype in looptypes:
    for i_ori, ori in enumerate(orientations):
        ## loading all data for this orientation
        amps = np.load("data/%s_atoms_%s_R%d_amplitudes_ori%d_0000.npy"%(\
                  sample, looptype, R, i_ori))
        amps1 = np.load("data/%s_atoms_%s_R%d_amplitudes1_ori%d_0000.npy"%(\
                  sample, looptype, R, i_ori))
        i_file = 1 
        while True:
            try:
                amps = np.vstack((amps, \
                          np.load("data/%s_atoms_%s_R%d_amplitudes_ori%d_%04d.npy"%(\
                          sample, looptype, R, i_ori, i_file))))
                amps1 = np.vstack((amps1, \
                          np.load("data/%s_atoms_%s_R%d_amplitudes1_ori%d_%04d.npy"%(\
                          sample, looptype, R, i_ori, i_file))))
            except IOError:
                break
            i_file += 1
        print "loaded %d files" % i_file
       
        ## sort the array
        amps = amps[np.argsort(amps[:,0])]
        amps1 = amps1[np.argsort(amps1[:,0])]

        ## calculate the intensities
        try:
            intensities += amps[:,1]**2
            intensities1 += amps1[:,1]**2
        except NameError:
            intensities = amps[:,1]**2
            intensities1 = amps1[:,1]**2

        print "done %s, orientation %d" % (looptype, i_ori)

    ## average over orientations        
    intensities /= len(orientations)
    intensities1 /= len(orientations)

    ax.plot(amps[:,0], amps[:,0]**4/R**2*intensities, label='%s,larger'%looptype)
    ax.plot(amps1[:,0], amps1[:,0]**4/R**2*intensities1, label='%s,smaller'%looptype)

ax.set_yscale("log", nonposy='mask')

plt.legend(loc='lower left')
plt.show()
