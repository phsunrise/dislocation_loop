import numpy as np
import sys, os
import matplotlib.pyplot as plt
from settings import basedir, Rlist
import glob

sample = 'W'
from W_parameters import a0, orientations

for R in Rlist: 
    print "sample %s, R=%d" % (sample, R)
    i_dir = 0
    datadir = basedir+"%s_R%d_amp2d_%d/"%(sample, R, i_dir)

    _file = np.load(datadir+"qarray_2d.npz")
    y_array = _file['hk_array']
    x_array = _file['l_array']
    x_mesh, y_mesh = np.meshgrid(x_array, y_array, indexing='xy')
    q4_mesh = ((x_mesh-(x_array[0]+x_array[-1])/2.)**2+\
               (y_mesh-(y_array[0]+y_array[-1])/2.)**2)**2
    xlabel = "l"
    ylabel = "h or k"
    extent=[x_array[0]*1.5-x_array[1]*0.5, x_array[-1]*1.5-x_array[-2]*0.5, \
            y_array[0]*1.5-y_array[1]*0.5, y_array[-1]*1.5-y_array[-2]*0.5]

    for looptype in ["int", "vac"]:
        intensities = np.load(datadir+"%s_atoms_intensity_2d_%s_R%d.npy"%(sample, looptype, R))

#        fig, ax = plt.subplots(1,1, figsize=(12,6))
#        img = ax.imshow(intensities*q4_mesh, origin='lower', \
#                        interpolation='nearest', extent=extent) 
#        fig.colorbar(img)
#        ax.set_title(r"%s, $q^4I$"%(looptype))
#        ax.set_xlabel(xlabel)
#        ax.set_ylabel(ylabel)
#        fig.savefig(datadir+"intensity_2d_%s_R%d_replot.pdf"%(looptype, R))
#        print "saved figure to", datadir+"intensity_2d_%s_R%d_replot.pdf"%(looptype, R)
#        plt.show()
#        plt.close()

        fig, ax = plt.subplots(1,1, figsize=(12,6))
        mode = "verticalslice"
        if mode == "horizontalslice":
            if len(y_array) % 2 == 0:
                intensityslice = 0.5*(intensities[len(y_array)/2-1, :]
                                     +intensities[len(y_array)/2, :])
            else:
                intensityslice = intensities[(len(y_array)-1)/2, :]
            img = ax.plot(x_array, \
                    (x_array-(x_array[0]+x_array[-1])/2.)**4*intensityslice)
            ax.set_title("%s"%(looptype))
            ax.set_xlabel(xlabel)
            ax.set_ylabel(r"$q^4 I$")
        elif mode == "verticalslice":
            from scipy.interpolate import interp2d
            intfunc = interp2d(x_array, y_array, intensities, \
                               kind='cubic', bounds_error=True)
            for xx in [0.00512695]:
                intensityslice = intfunc(xx, y_array).flatten()
                xdata = (y_array-(y_array[0]+y_array[-1])/2.)*2.*np.pi/a0*np.sqrt(2)
                ydata = xdata**4*intensityslice
                ax.plot(xdata, ydata, label="%.3f"%xx)
                np.savez(basedir+"fit/R%d_%s.npz"%(R, looptype), \
                         q=xdata, q4I=ydata)
                        
            ax.set_title("%s"%(looptype))
            ax.set_xlabel(r"$q$ ($\AA^{-1}$)")
            ax.set_ylabel(r"$q^4 I$")
            #ax.set_ylim(0., (70. if looptype=='int' else 60.))
            #ax.set_ylim(0., 60.)
            ax.legend()
        elif mode == "diagonal":
            intensityslice = [intensities[i,i] for i in xrange(len(x_array))]
            img = ax.plot(x_array, \
                    (x_array-(x_array[0]+x_array[-1])/2.)**4*intensityslice)
            ax.set_title("%s"%(looptype))
            ax.set_xlabel(xlabel)
            ax.set_ylabel(r"$q^4 I$")
        #ax.set_xlim(-0.15, 0.15)
        #fig.savefig(datadir+"intensity_1d_%s_R%d_replot.pdf"%(looptype, R))
        #print "saved figure to", datadir+"intensity_1d_%s_R%d_replot.pdf"%(looptype, R)
        plt.show()
        plt.close()

