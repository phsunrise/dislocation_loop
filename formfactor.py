import numpy as np
from W_parameters import *

xyzf_list = [(0.5, 0.5, -1., 1./4), (-0.5, 0.5, -1., 1./4), (0.5, -0.5, -1., 1./4), (-0.5, -0.5, -1., 1./4),\
             (0., 0., 0., 1.), (1., 0., 0., 1./4), (-1., 0., 0., 1./4), (0., 1., 0., 1./4), (0., -1., 0., 1./4),\
             (0.5, 0.5, 1., 1./4), (-0.5, 0.5, 1., 1./4), (0.5, -0.5, 1., 1./4), (-0.5, -0.5, 1., 1./4)]
xyzf_list = np.array(xyzf_list)

xyzf_list2 = []
for x in [-1., 0., 1.]:
    for y in [-1., 0., 1.]:
        for z in [-1., 0., 1.]:
            f = 1.
            if abs(x) == 1.:
                f /= 2.
            if abs(y) == 1.:
                f /= 2.
            if abs(z) == 1.:
                f /= 2.
            xyzf_list2.append([x,y,z,f])
xyzf_list2 = np.array(xyzf_list2)

def formfactor(qloop, tier):
    ff = 0.
    for xyzf in xyzf_list:
        r = xyzf[0]*ex_p + xyzf[1]*ey_p + xyzf[2]*ez_p
        ff += xyzf[3]*np.cos(qloop.dot(r))
    if tier == 2:
        return ff

    ff2 = 0.
    for xyzf in xyzf_list2:
        r = xyzf[0]*(ex_p+ey_p) + xyzf[1]*(ex_p-ey_p) + xyzf[2]*ez_p*2
        ff2 += xyzf[3]*np.cos(qloop.dot(r))
    if tier == 3:
        return ff*ff2
