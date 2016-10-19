import numpy as np
import sys, os
import matplotlib.pyplot as plt
from settings import basedir
from formfactor import formfactor

sample = 'W'
from W_parameters import *
looptypes = ['vac', 'int']
colors = ['r', 'b', 'g']

fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(1, 1, 1)

for i_R, R in enumerate([10., 20., 40.]):
    print "sample:", sample
    print "R=%.1f"%(R)
    datadir = basedir+"%s_R%d/"%(sample, R)
    q_array = np.load(datadir+"q_array.npy")

    for looptype in looptypes:
        print "starting looptype %s..." % looptype
        intensities = [] 
        intensities1 = []
        intensities2 = []
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
            amplitudes1 = np.copy(amplitudes)
            amplitudes2 = np.copy(amplitudes)

            amplitudes[:,1] += np.load(datadir+"%s_atoms_amplitude_%s_T1_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))[:,1]*2.
            amplitudes1[:,1] += np.load(datadir+"%s_atoms_amplitude1_%s_T1_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))[:,1]*2.
            amplitudes2[:,1] += np.load(datadir+"%s_atoms_amplitude2_%s_T1_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))[:,1]*2.
            amplitudes[:,1] += np.load(datadir+"%s_atoms_amplitude_%s_T2_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))[:,1]*ff_tier2*2.
            amplitudes1[:,1] += np.load(datadir+"%s_atoms_amplitude1_%s_T2_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))[:,1]*ff_tier2*2.
            amplitudes2[:,1] += np.load(datadir+"%s_atoms_amplitude2_%s_T2_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))[:,1]*ff_tier2*2.
            amplitudes[:,1] += np.load(datadir+"%s_atoms_amplitude_%s_T3_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))[:,1]*ff_tier3*2.
            amplitudes1[:,1] += np.load(datadir+"%s_atoms_amplitude1_%s_T3_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))[:,1]*ff_tier3*2.
            amplitudes2[:,1] += np.load(datadir+"%s_atoms_amplitude2_%s_T3_R%d_ori%d_combined.npy"%(\
                                sample, looptype, R, i_ori))[:,1]*ff_tier3*2.

            intensities.append(amplitudes[:,1]**2)
            intensities1.append(amplitudes1[:,1]**2)
            intensities2.append(amplitudes2[:,1]**2)

        # average over all orientations
        intensities = np.mean(np.array(intensities), axis=0)
        intensities1 = np.mean(np.array(intensities1), axis=0)
        intensities2 = np.mean(np.array(intensities2), axis=0)

        # plot q^4/R^2*I 
        if looptype == 'vac':
            linestyle = '--'
        elif looptype == 'int':
            linestyle = '-'
        ax.plot(q_array, intensities*abs(q_array)**6/R**0, \
            color=colors[i_R], ls=linestyle, \
            label=r"$R=%d\mathrm{\AA}$, %s"%(R, looptype))
        #ax.plot(q_array, intensities/R**6, \
        #    color=colors[i_R], ls=linestyle, \
        #    label=r"$R=%d\mathrm{\AA}$, %s"%(R, looptype))
        #ax.plot(q_array, intensities1*abs(q_array)**4/R**2, color='b', ls=linestyle)
        #ax.plot(q_array, intensities2*abs(q_array)**4/R**2, color='g', ls=linestyle)
        np.save(datadir+"%s_atoms_intensity_R%d.npy"%(sample, R), intensities)
        np.save(datadir+"%s_atoms_intensity1_R%d.npy"%(sample, R), intensities1)
        np.save(datadir+"%s_atoms_intensity2_R%d.npy"%(sample, R), intensities2)

ax.set_xlim(-0.5, 0.5)
#ax.set_ylim(0., 20.)
ax.set_xlabel(r"$q$ (220)")
ax.set_ylabel(r"$q^4/R^2 I$")
ax.legend(loc='upper right')
plt.tight_layout()
fig.savefig("intensity_%s.pdf"%sample)
plt.show()
