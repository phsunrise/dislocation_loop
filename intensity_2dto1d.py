import os, sys
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from settings import plot_dir
from W_parameters import a0

#R_ind_list = [(20., 2, 'r'), (40., 2, 'g'), (60., 0, 'b'), (80., 0, 'm')]
R_ind_list = [(20., 2, 'r'), (40., 2, 'g'), (60., 0, 'b'), (80., 0, 'm'), \
              (20., 3, 'c'), (60., 1, 'k')]

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)

for R, ind, color in R_ind_list:
    for looptype in ['vac', 'int']:
        datadir = "data/W_R%d_amp2d_%d/" % (R, ind)

        intensity = np.load(datadir+"W_atoms_intensity_2d_%s_R%d.npy" % (looptype, R))
        data = np.load(datadir+"qarray_2d.npz")
        q_array = data['q_array']
        q_array_norm = np.sum(q_array**2, axis=2)**0.5
        h = data['h']
        h = h/(2.*np.pi/a0)
        chi_array = data['chi_array']/np.pi*180.
        chi0_ind = np.where(chi_array==0.)[0][0]
        tth_array = data['tth_array']/np.pi*180.

        chi, tth = np.meshgrid(chi_array, tth_array, indexing='xy')

        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #ax.plot_surface(chi, tth, intensity, \
        #        rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0)
        #plt.show()

        ## modify intensity by q^4*I/R^2
        intensity = q_array_norm**4*intensity/R**2
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        img = ax2.imshow(intensity, origin='lower', \
                  extent=[chi_array[0], chi_array[-1], tth_array[0], tth_array[-1]], \
                  interpolation='nearest')
        fig2.colorbar(img)
        ax2.set_title("R%d, %s" % (R, looptype))
        fig2.savefig(plot_dir+"q4I_R%d_h%d%d%d_%s_2d.pdf"%(\
                        R, h[0], h[1], h[2], looptype))
        plt.close(fig2)

        ls = '--' if looptype == 'vac' else '-'
        q_proj_on_h = q_array[:, chi0_ind].dot(h)/np.linalg.norm(h)
        ax.plot(q_proj_on_h, intensity[:,chi0_ind], \
                ls=ls, label="R%d %s"%(R, looptype), color=color)

        ## get maximum
        max_ind = np.unravel_index(intensity.argmax(), intensity.shape)
        print intensity[max_ind]
        print chi[max_ind], tth[max_ind]

ax.legend(bbox_to_anchor=[0., 1.02, 1., 0.102], ncol=4, \
          mode="expand", borderaxespad=0.)
fig.savefig(plot_dir+"q4I_h%d%d%d_all_1d.pdf"%(h[0],h[1],h[2]))
plt.show()
