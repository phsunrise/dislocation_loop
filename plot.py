import matplotlib.pyplot as plt
import numpy as np

sample = 'Al'
amps = np.load("data/%s_amplitudes.npy"%sample)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(amps[:,0], amps[:,0]**2*amps[:,1]**2)
ax.plot(-amps[:,0], amps[:,0]**2*amps[:,1]**2)
ax.set_xlim(0.1, 5)
ax.set_xscale("log", nonposx='mask')
ax.set_yscale("log", nonposy='mask')
plt.show()
