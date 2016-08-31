import matplotlib.pyplot as plt
import numpy as np

plotting = 2 

sample = 'Al'
amps = np.load("data/%s_atoms_amplitudes.npy"%sample)
amps1 = np.load("data/%s_atoms_amplitudes_1.npy"%sample)

fig = plt.figure()
ax = fig.add_subplot(111)
if plotting == 2:
    ax.plot(amps[:,0], amps[:,0]**2*amps[:,1]**2, \
            label='larger, pos')
    ax.plot(-amps[:,0], amps[:,0]**2*amps[:,1]**2, \
            label='larger, neg')
    ax.plot(amps1[:,0], amps1[:,0]**2*amps1[:,1]**2, \
            label='smaller, pos')
    ax.plot(-amps1[:,0], amps1[:,0]**2*amps1[:,1]**2, \
            label='smaller, neg')
    ax.set_xlim(0.1, 5)
    ax.set_xscale("log", nonposx='mask')
    ax.set_yscale("log", nonposy='mask')
elif plotting == 4:
    ax.plot(amps[:,0], amps[:,0]**4*amps[:,1]**2, label='larger')
    ax.plot(amps1[:,0], amps1[:,0]**4*amps1[:,1]**2, label='smaller')

plt.legend(loc='lower left')
plt.show()
