import numpy as np
import matplotlib.pyplot as plt
from info import *

sample = 'W'
looptype = 'vac'
R = 20.
datadir = "%s_R%d/" % (sample, R)

lines = []
for i_file in xrange(NFILES):
    try:
        sdata = np.load(datadir+"%s_atoms_s_%s_T1_R%d_%04d.npy"%(\
            sample, looptype, R, i_file))
    except IOError:
        continue

    counter = 0
    for line in sdata:
        x,y,z = line[0:3]
        if (4.*np.pi/40. < 
                np.arccos(z/np.linalg.norm([x,y,z])) < 
                5.*np.pi/40. and 
                3.*np.pi/40. < np.arctan2(y,x) < 4.*np.pi/40.):
            lines.append(line)
            counter += 1
    if counter > 0:
        print "done file %d, added %d points" % (i_file, counter)

lines = np.array(lines).T
np.save("checks_results.npy", lines)
plt.scatter(np.sqrt(lines[0]**2+lines[1]**2+lines[2]**2),\
            np.sqrt(lines[3]**2+lines[4]**2+lines[5]**2))
ax = plt.gca()
ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')

plt.show()
