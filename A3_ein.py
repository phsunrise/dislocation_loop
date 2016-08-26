import numpy as np
import csv
import os, sys
os.chdir("data")

sample = sys.argv[1]
os.system("rm %s_A3?_ein_ori*.*" % sample)
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

sdata = np.load("%s_s_0.npy"%sample)
rank = 1 
while os.path.isfile("%s_s_%d.npy"%(sample, rank)):
    sdata = np.vstack((sdata, \
         np.load("%s_s_%d.npy"%(sample, rank))))
    rank += 1
print "read %d files in total" % rank

for i_orientation, orientation in enumerate(orientations):
    h1 = np.einsum('ij,j', orientation, h)
    eq1 = np.einsum('ij,j', orientation, eq)
    for qR in np.append(np.linspace(-5., -0.1, 50), np.linspace(0.1, 5., 50)):
        A3s = []
        A3a = []
        q = qR/R*eq1
        qloop = np.einsum('ij,j', rot, q)
        K = h1+q
        Kloop = np.einsum('ij,j', rot, K)
        sintegrand = 0.
        aintegrand = 0.
        for i in xrange(len(sdata)):
            rho = sdata[i,0]*R
            if rho == 0.:
                continue
            th = sdata[i,1]
            z = sdata[i,2]*R
            deltaV = 2* rho*0.1*R * 0.1*R * 2.*np.pi/20 
            Ks = Kloop.dot(sdata[i,3:6]*np.linalg.norm(B))
            qr = qloop.dot([rho*np.cos(th), rho*np.sin(th), z])
            sintegrand += deltaV/Vc*np.cos(qr)*(np.cos(Ks)-1)
            aintegrand += -deltaV/Vc*np.sin(qr)*(np.sin(Ks)-Ks)
        A3s.append([qR, sintegrand])
        A3a.append([qR, aintegrand])

        try:
            A3s = np.vstack((np.load("%s_A3s_ein_ori%d.npy"%(\
                                            sample, i_orientation), A3s))
            A3a = np.vstack((np.load("%s_A3a_ein_ori%d.npy"%(\
                                            sample, i_orientation), A3a))
        except IOError:
            A3s = np.array(A3s)
            A3a = np.array(A3a)
        np.save("%s_A3s_ein_ori%d.npy"%(sample, i_orientation), A3s)
        np.save("%s_A3a_ein_ori%d.npy"%(sample, i_orientation), A3a)
        print "done qR = %f" % qR

    print "done orientation %d" % i_orientation
