import numpy as np
import sys, os
import matplotlib.pyplot as plt
from settings import basedir
import glob
from formfactor import formfactor

sample = 'W'
from W_parameters import R 
print "sample %s, R=%d" % (sample, R)
i_dir = 0
while os.path.isdir(basedir + "%s_R%d_amp2d_%d/"%(sample, R, i_dir)):
    print "folder %s exists, containing %d files" % (\
            basedir+"%s_R%d_amp2d_%d/"%(sample, R, i_dir), \
            len(glob.glob(basedir+"%s_R%d_amp2d_%d/%s_atoms_amplitude_2d_*_T*_R*_ori*_????.npy"%(\
                sample, R, i_dir, sample))))
    i_dir += 1
val = raw_input("enter folder number: ")
i_dir = int(val)
datadir = basedir+"%s_R%d_amp2d_%d/"%(sample, R, i_dir)

sys.path.insert(0, datadir)
from W_parameters import *
print "D/R=%.1f" % (D/R)
_file = np.load(datadir+"qarray_2d.npz")
q_array = _file['q_array']
tth_array = _file['tth_array']/np.pi*180.
chi_array = _file['chi_array']/np.pi*180.
chi, tth = np.meshgrid(chi_array, tth_array, indexing='xy')

looptypes = ['int', 'vac']
colors = ['r', 'b', 'g', 'c']

for looptype in looptypes:
    print "starting looptype %s..." % looptype
    intensities = [] 
    for i_ori, ori in enumerate(orientations):
        print "starting orientation %d..." % i_ori

        # get the form factor
        ff_tier2 = np.zeros((q_array.shape[0], q_array.shape[1]))
        ff_tier3 = np.zeros((q_array.shape[0], q_array.shape[1]))
        for index, _ in np.ndenumerate(ff_tier2):
            q = q_array[index] 
            qloop = np.einsum('ij,jk,k', ori, rot, q)
            ff_tier2[index] = formfactor(qloop, 2)
            ff_tier3[index] = formfactor(qloop, 3)

        amplitudes = np.load(datadir+"%s_atoms_core_2d_%s_R%d_ori%d.npy"%(\
                                sample, looptype, R, i_ori)) 
        #amplitudes *= 0

        amplitudes += np.load(datadir+"%s_atoms_amplitude_2d_%s_T1_R%d_ori%d_combined.npy"%(\
                            sample, looptype, R, i_ori))*2.
        amplitudes += np.load(datadir+"%s_atoms_amplitude_2d_%s_T2_R%d_ori%d_combined.npy"%(\
                            sample, looptype, R, i_ori))*ff_tier2*2.
        amplitudes += np.load(datadir+"%s_atoms_amplitude_2d_%s_T3_R%d_ori%d_combined.npy"%(\
                            sample, looptype, R, i_ori))*ff_tier3*2.

        intensities.append(amplitudes**2)

        ## plot individual orientations 
        #fig = plt.figure(figsize=(12,6))
        #ax = fig.add_subplot(1, 1, 1)
        #img = ax.imshow(intensities, origin='lower', \
        #                interpolation='nearest', \
        #          extent=[chi_array[0]*1.5-chi_array[1]*0.5, \
        #                  chi_array[-1]*1.5-chi_array[-2]*0.5, \
        #                  tth_array[0]*1.5-tth_array[1]*0.5, \
        #                  tth_array[-1]*1.5-tth_array[-2]*0.5])
        #fig.colorbar(img)

    # average over all orientations
    intensities = np.array(intensities)
    intensities = np.mean(intensities, axis=0)

    # plot diffuse intensity 
    #if looptype == 'vac':
    #    linestyle = '--'
    #elif looptype == 'int':
    #    linestyle = '-'

    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(1, 1, 1)
    img = ax.imshow(intensities, origin='lower', \
                    interpolation='nearest', \
              extent=[chi_array[0]*1.5-chi_array[1]*0.5, \
                      chi_array[-1]*1.5-chi_array[-2]*0.5, \
                      tth_array[0]*1.5-tth_array[1]*0.5, \
                      tth_array[-1]*1.5-tth_array[-2]*0.5])
    fig.colorbar(img)
    ax.set_title(looptype)
    ax.set_xlabel("chi (deg)")
    ax.set_ylabel("tth (deg)")
    fig.savefig(datadir+"intensity_2d_%s_R%d.pdf"%(looptype, R))
    print "saved figure to", datadir+"intensity_2d_%s_R%d.pdf"%(looptype, R)
    plt.show()
    plt.close()
    
    np.save(datadir+"%s_atoms_intensity_2d_%s_R%d.npy"%(sample, looptype, R), intensities)
