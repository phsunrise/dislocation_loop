import numpy as np
import matplotlib.pyplot as plt
from settings import basedir

sample = 'W'

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(1,1,1)

for i_R, R in enumerate([20.]):
    datadir = basedir+"%s_R%d/"%(sample, R)
    q_array = np.load(datadir+"q_array.npy")
    for looptype in ['vac', 'int']:
        intensities = np.load(datadir+"%s_atoms_intensity_%s_R%d.npy"%(\
                                sample, looptype, R))
        intensities_poly = np.load(datadir+"%s_atoms_intensity_poly_%s_R%d.npy"%(\
                                sample, looptype, R))*1.e5

        linestyle = '--' if looptype=='vac' else '-'
        ax.plot(q_array, intensities, 'r'+linestyle, \
                label="single, %s"%looptype)
        ax.plot(q_array, intensities_poly, 'b'+linestyle, \
                label="poly*1e5, %s"%looptype)

ax.set_xlim(-0.3, 0.3)
ax.set_ylim(0., 1.5e9)
ax.set_xlabel(r"$q$ (220)")
ax.set_ylabel(r"$I$")
ax.legend(loc='upper right')
plt.tight_layout()
fig.savefig("intensity_compare_%s.pdf"%sample)
plt.show()
