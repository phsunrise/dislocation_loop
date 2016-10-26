import numpy as np
import sys, os
import matplotlib.pyplot as plt
from settings import basedir
from formfactor import formfactor

sample = 'W'
from W_parameters import *
looptypes = ['int', 'vac']
colors = ['r', 'b', 'g', 'c']

fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(1, 1, 1)

for i_R, R in enumerate([20.]):
    print "sample:", sample
    print "R=%.1f"%(R)
    datadir = basedir+"%s_R%d/"%(sample, R)
    q_array = np.load(datadir+"q_array.npy")

    for looptype in looptypes:
        print "starting looptype %s..." % looptype
        intensities = [] 
        for i_ori, ori in enumerate(orientations):
            print "starting orientation %d..." % i_ori

            # get the form factor
            ff_tier2 = []
            ff_tier3 = []
            for i_qval, qval in enumerate(q_array):
                q = qval * np.einsum('ij,j', ori, eq)
                qloop = np.einsum('ij,j', rot, q)
                ff_tier2.append(formfactor(qloop, 2))
                ff_tier3.append(formfactor(qloop, 3))
            ff_tier2 = np.array(ff_tier2)
            ff_tier3 = np.array(ff_tier3)

            amplitudes = np.load(datadir+"%s_atoms_core_%s_R%d_ori%d.npy"%(\
                                    sample, looptype, R, i_ori)) 
            #amplitudes *= 0

            amplitudes[:,1] += np.load(datadir+"%s_atoms_amplitude_%s_T1_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))[:,1]*2.
            amplitudes[:,1] += np.load(datadir+"%s_atoms_amplitude_%s_T2_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))[:,1]*ff_tier2*2.
            amplitudes[:,1] += np.load(datadir+"%s_atoms_amplitude_%s_T3_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))[:,1]*ff_tier3*2.

            intensities.append(amplitudes[:,1]) ## here only appending amplitudes

        # average over all orientations
        intensities = np.array(intensities)
        print intensities[:, :5]
        print intensities[:, -5:]
        intensities = np.mean(intensities, axis=0)**2 ## here turning amplitude into intensity

        # plot q^4/R^2*I 
        if looptype == 'vac':
            linestyle = '--'
        elif looptype == 'int':
            linestyle = '-'
        ax.plot(q_array, intensities*abs(q_array)**4/R**2, \
            color=colors[i_R], ls=linestyle, \
            label=r"$R=%d\mathrm{\AA}$, %s"%(R, looptype))
        #ax.plot(q_array, intensities/R**6, \
        #    color=colors[i_R], ls=linestyle, \
        #    label=r"$R=%d\mathrm{\AA}$, %s"%(R, looptype))
        np.save(datadir+"%s_atoms_intensity_%s_R%d.npy"%(sample, looptype, R), intensities)

ax.set_xlim(-1.0, 1.0)
#ax.set_ylim(0., 5.)
ax.set_xlabel(r"$q$ (220)")
ax.set_ylabel(r"$q^4/R^2 I$")

ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=4, \
          fancybox=True, shadow=True)

plt.tight_layout()
fig.savefig("intensity_%s.pdf"%sample)
plt.show()
