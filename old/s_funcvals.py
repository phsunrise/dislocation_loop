import sys, os
import numpy as np
from scipy import interpolate
from W_parameters import *
import time

do_check = True 

npoints = 2000

if not do_check:
    vals = np.zeros((npoints, npoints, 3))
    maxerr = 0.
    maxerr_vals = [] 
    maxerr_ind = []
    for i_u, u in enumerate(\
            np.linspace(-1., 1., npoints, endpoint=False)+1./npoints):
        for i_th, th in enumerate(\
                np.linspace(0., 2.*np.pi, npoints, endpoint=False)):
            x = np.sqrt(1.-u**2)*np.cos(th) 
            y = np.sqrt(1.-u**2)*np.sin(th) 
            z = u
            e = x*ex+y*ey+z*ez
            eloop = np.array([x,y,z])

            g = np.zeros((3,3))
            sum1 = sum(e**2/(C44+d*e**2))
            for i in xrange(3):
                for j in xrange(3):
                    g[i,j] = ((1 if i==j else 0)/(C44+d*e[i]**2) - \
                        (e[i]*e[j]/(C44+d*e[i]**2)/(C44+d*e[j]**2)*\
                         (C44+C12)/(1.+(C44+C12)*sum1)
                        )
                      )
            gloop =  np.einsum('ij,kl,jl', rot, rot, g)
            vals[i_u, i_th] = np.einsum('ij,i,jk', Ploop, eloop, gloop)
            if i_th > 0:
                if maxerr < max(abs(vals[i_u, i_th-1]-vals[i_u, i_th])):
                    maxerr = max(abs(vals[i_u, i_th-1]-vals[i_u, i_th]))
                    maxerr_vals = [vals[i_u, i_th-1], vals[i_u, i_th]]
                    maxerr_ind = [i_u, i_th, [x, y, z]]
            if i_u > 0:
                if maxerr < max(abs(vals[i_u-1, i_th]-vals[i_u, i_th])):
                    maxerr = max(abs(vals[i_u-1, i_th]-vals[i_u, i_th])) 
                    maxerr_vals = [vals[i_u-1, i_th], vals[i_u, i_th]]
                    maxerr_ind = [i_u, i_th, [x, y, z]]
        print "done %d/%d" % (i_u+1, npoints)

    np.save("s_funcvals.npy", vals)
    print "max error =", maxerr
    print maxerr_vals
    print maxerr_ind

else:
    vals = np.load("s_funcvals.npy")
    u_array = np.linspace(-1., 1., npoints, endpoint=False)+1./npoints
    th_array = np.linspace(0., 2.*np.pi, npoints, endpoint=False)
    funcs = []
    for ind in xrange(3):
        funcs.append(interpolate.interp2d(th_array, u_array, \
                vals[:,:,ind], kind='cubic'))
    print "done interpolation"

    errs = np.zeros((npoints-1, npoints-1, 3))
    for i_u in xrange(npoints-1): 
        for i_th in xrange(npoints-1):
            u = np.mean(u_array[i_u: i_u+2])
            th = np.mean(th_array[i_th: i_th+2])
            x = np.sqrt(1.-u**2)*np.cos(th) 
            y = np.sqrt(1.-u**2)*np.sin(th) 
            z = u

            time0 = time.time()
            e = x*ex+y*ey+z*ez
            eloop = np.array([x,y,z])

            g = np.zeros((3,3))
            sum1 = sum(e**2/(C44+d*e**2))
            for i in xrange(3):
                for j in xrange(3):
                    g[i,j] = ((1 if i==j else 0)/(C44+d*e[i]**2) - \
                        (e[i]*e[j]/(C44+d*e[i]**2)/(C44+d*e[j]**2)*\
                         (C44+C12)/(1.+(C44+C12)*sum1)
                        )
                      )
            gloop =  np.einsum('ij,kl,jl', rot, rot, g)
            val = np.einsum('ij,i,jk', Ploop, eloop, gloop)
            time1 = time.time()

            val_interp = np.array([funcs[ind](th, u) 
                            for ind in xrange(3)]).ravel()
            time2 = time.time()
            print time1-time0, time2-time1
            sys.exit(0)
            errs[i_u, i_th] = val_interp - val 
        print "done %d/%d" % (i_u+1, npoints-1)

    np.save("s_funcerrs.npy", errs)
    import matplotlib.pyplot as plt
    for ind in xrange(3):
        plt.figure()
        plt.imshow(errs[:,:,ind])
        plt.colorbar()

        plt.figure()
        plt.imshow(vals[:,:,ind])
        plt.colorbar()
    plt.show()

