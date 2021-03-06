import numpy as np
import sys, os
import matplotlib.pyplot as plt
from settings import basedir, ff_dir
from formfactor import formfactor
import time

sample = 'W'
from W_parameters import *
looptypes = ['vac', 'int']
colors = ['r', 'b', 'g', 'c']

# get uniform distributed points on sphere
sph_array = np.load("uniformsphere.npy")

fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(1, 1, 1)

for i_R, R in enumerate([10., 20.]):
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
            ff_tier2 = np.load(ff_dir+"ff_tier2_ori%d.npy"%(i_ori)) 
            ff_tier3 = np.load(ff_dir+"ff_tier3_ori%d.npy"%(i_ori)) 

            amplitudes = np.load(datadir+"%s_atoms_core_poly_%s_R%d_ori%d.npy"%(\
                                    sample, looptype, R, i_ori)) 

            amplitudes += np.load(datadir+"%s_atoms_amplitude_poly_%s_T1_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))*2.
            amplitudes += np.load(datadir+"%s_atoms_amplitude_poly_%s_T2_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))*ff_tier2*2.
            amplitudes += np.load(datadir+"%s_atoms_amplitude_poly_%s_T3_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))*ff_tier3*2.

            intensities.append(amplitudes**2)

        # average over all orientations
        intensities = np.mean(np.array(intensities), axis=(0,2))
        #intensities1 = np.mean(np.array(intensities)[:,:,0::4], axis=(0,2))
        #intensities2 = np.mean(np.array(intensities)[:,:,1::4], axis=(0,2))
        #intensities3 = np.mean(np.array(intensities)[:,:,2::4], axis=(0,2))
        #intensities4 = np.mean(np.array(intensities)[:,:,3::4], axis=(0,2))

        # plot q^4/R^2*I 
        if looptype == 'vac':
            linestyle = '--'
        elif looptype == 'int':
            linestyle = '-'
        #ax.plot(q_array, intensities1/R**2, \
        #    color=colors[i_R], ls=linestyle, \
        #    label=r"$R=%d\mathrm{\AA}$, %s"%(R, looptype))
        #ax.plot(q_array, intensities2/R**2, \
        #    color=colors[i_R+1], ls=linestyle, \
        #    label=r"$R=%d\mathrm{\AA}$, %s"%(R, looptype))
        #ax.plot(q_array, intensities3/R**2, \
        #    color=colors[i_R+2], ls=linestyle, \
        #    label=r"$R=%d\mathrm{\AA}$, %s"%(R, looptype))
        #ax.plot(q_array, intensities4/R**2, \
        #    color=colors[i_R+3], ls=linestyle, \
        #    label=r"$R=%d\mathrm{\AA}$, %s"%(R, looptype))
        ax.plot(q_array, intensities/R**3, \
            color=colors[i_R], ls=linestyle, \
            label=r"$R=%d\mathrm{\AA}$, %s"%(R, looptype))
        np.save(datadir+"%s_atoms_intensity_poly_%s_R%d.npy"%(sample, looptype, R), intensities)

ax.set_xlim(-1.0, 1.0)
#ax.set_ylim(0., 1.0)
ax.set_xlabel(r"$q$ (220)")
ax.set_ylabel(r"$|q|^{1.5}/R^2 I$")

ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=4, \
          fancybox=True, shadow=True)

plt.tight_layout()
fig.savefig("intensity_poly_%s.pdf"%sample)
plt.show()
