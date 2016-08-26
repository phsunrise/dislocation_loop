import numpy as np
import matplotlib.pyplot as plt
import os, sys

colors = ['r', 'b']
datadir = "/afs/slac.stanford.edu/u/xo/phsun/data1/"
datadir = "data/"
for i_sample, sample in enumerate(['Al', 'Cu']):
    sdata = np.load(datadir+"%s_s_0.npy"%sample)
    rank = 1
    while os.path.isfile(datadir+"%s_s_%d.npy"%(sample, rank)):
        sdata = np.vstack((sdata, \
                 np.load(datadir+"%s_s_%d.npy"%(sample, rank))))
        rank += 1

    sdata_r0 = []
    for i in xrange(len(sdata)):
        if sdata[i, 0] == 0.:
            sdata_r0.append(sdata[i])
    sdata_r0 = np.array(sdata_r0)

    plt.scatter(sdata_r0[:,2], sdata_r0[:,5], c=colors[i_sample], \
            label=sample)

plt.legend()
plt.gca().invert_yaxis()
plt.savefig("uz_r0.pdf")
plt.show()
