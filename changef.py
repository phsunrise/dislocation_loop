import numpy as np
import sys, os
from W_parameters import *

print "currently disabled"

'''
for looptype in ['int', 'vac']:
    for i_file in xrange(8000):
        try:
            sdata = np.load("W_R20/W_atoms_s_%s_T1_R20_%04d.npy"%(
                looptype, i_file))
        except IOError:
            break

        sdata1 = []
        for line in sdata:
            if ((looptype == 'int' and abs(line[2]/a2-0.5)<1.e-10) 
                  or (looptype == 'vac' and abs(line[2]/a2-1.)<1.e-10)):
                line[6] += 0.5
            sdata1.append(line)

        np.save("W_R20/W_atoms_s_%s_T1_R20_%04d_changed.npy"%(looptype, i_file),\
            sdata1)
        print "done file %d" % i_file
'''
