import numpy as np

sdata = np.load("data/Al_s.npy")
sdata_conv = []
for i in xrange(len(sdata)):
    if sdata[i,2] == 0.: 
        sdata_conv.append(sdata[i])

np.save("data/Al_s_th0.npy", sdata_conv)
