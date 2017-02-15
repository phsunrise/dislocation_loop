import numpy as np
import sys, os
import matplotlib.pyplot as plt
from settings import basedir
import glob

sample = 'W'
from W_parameters import a0, R 
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

del sys.modules["W_parameters"]
sys.path.insert(0, datadir)
from W_parameters import R, funcform, funcparams, orientations
print funcform, "params =", funcparams

_file = np.load(datadir+"qarray_2d.npz")
q_array = _file['q_array']
if 'tth_array' in _file:
    y_array = _file['tth_array']/np.pi*180.
    x_array = _file['chi_array']/np.pi*180.
    x_mesh, y_mesh = np.meshgrid(x_array, y_array, indexing='xy')
    xlabel = "chi (deg)"
    ylabel = "tth (deg)"
elif 'hk_array' in _file:
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

    fig, ax = plt.subplots(1,1, figsize=(12,6))
    img = ax.imshow(intensities*q4_mesh, origin='lower', \
                    interpolation='nearest', extent=extent) 
    fig.colorbar(img)
    ax.set_title(r"%s, $q^4I$"%(looptype))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.savefig(datadir+"intensity_2d_%s_R%d_replot.pdf"%(looptype, R))
    print "saved figure to", datadir+"intensity_2d_%s_R%d_replot.pdf"%(looptype, R)
    plt.show()
    plt.close()

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
        #for xx in [-0.01, -0.005, 0., 0.005, 0.01]:
        for xx in [0.]:
            intensityslice = intfunc(xx, y_array).flatten()
            ax.plot((y_array-(y_array[0]+y_array[-1])/2.)*2.*np.pi/a0, \
                    np.abs(y_array-(y_array[0]+y_array[-1])/2.)**4\
                        *intensityslice, \
                    label="%.3f"%xx)
                    
        ax.set_title("%s"%(looptype))
        ax.set_xlabel(r"$q$ ($\AA^{-1}$)")
        ax.set_ylabel(r"$q^4 I$")
        ax.legend()
    elif mode == "diagonal":
        intensityslice = [intensities[i,i] for i in xrange(len(x_array))]
        img = ax.plot(x_array, \
                (x_array-(x_array[0]+x_array[-1])/2.)**4*intensityslice)
        ax.set_title("%s"%(looptype))
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$q^4 I$")
    #ax.set_xlim(-0.15, 0.15)
    fig.savefig(datadir+"intensity_1d_%s_R%d_replot.pdf"%(looptype, R))
    print "saved figure to", datadir+"intensity_1d_%s_R%d_replot.pdf"%(looptype, R)
    plt.show()
    plt.close()

