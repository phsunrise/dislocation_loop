import numpy as np
import os, sys

sdata_conv = []
rank = 0
while os.path.isfile("data/Al_s_%d.npy"%rank): 
    sdata = np.load("data/Al_s_%d.npy"%rank)
    for i in xrange(len(sdata)):
        if sdata[i,1] == 0.: 
            print sdata[i]
            sdata_conv.append(sdata[i])
    print "done rank %d" % rank
    rank += 1

np.save("data/Al_s_th0.npy", sdata_conv)
