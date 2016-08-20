import numpy as np
import csv
from Al_parameters import *

h = 2.*np.pi/a0*np.array([4.,0.,0.])

sdata = np.load("Al_s.npy")

for qR in np.logspace(np.log10(0.1), np.log10(5), num=20):
    A3s = []
    A3a = []
    for qth in np.linspace(0., np.pi/2., 11):
        for qphi in np.linspace(0., 2*np.pi, 20, endpoint=False):
            q = qR/R*np.array([np.sin(qth)*np.cos(qphi), 
                    np.sin(qth)*np.sin(qphi), np.cos(qth)])
            qloop = np.einsum('ij,j', rot, q)
            K = h+q
            Kloop = np.einsum('ij,j', rot, K)
            sintegrand = 0.
            aintegrand = 0.
            for i in xrange(len(sdata)):
                rho = sdata[i,0]
                if rho == 0.:
                    continue
                th = sdata[i,1]
                z = sdata[i,2]
                deltaV = 2* rho*0.1*R * 0.1*R * 2*np.pi/20 
                Ks = Kloop.dot(sdata[i,3:6])
                qr = qloop.dot([rho*np.cos(th), rho*np.sin(th), z])
                sintegrand += deltaV*(1./Vc*np.linalg.norm(B)*\
                        np.cos(qr)*(np.cos(Ks)-1))
                aintegrand += -deltaV*(1./Vc*np.linalg.norm(B)*\
                        np.sin(qr)*(np.sin(Ks)-Ks))
            A3s.append([q[0], q[1], q[2], sintegrand])
            A3a.append([q[0], q[1], q[2], aintegrand])

    try:
        A3s = np.vstack((np.load("Al_A3s.npy"), A3s))
        A3a = np.vstack((np.load("Al_A3a.npy"), A3a))
    except IOError:
        A3s = np.array(A3s)
        A3a = np.array(A3a)
    np.save("Al_A3s.npy", A3s)
    np.save("Al_A3a.npy", A3a)

with open("Al_A3s.csv", 'wb') as f:
    wr = csv.writer(f, delimiter=' ', \
                    quoting=csv.QUOTE_NONNUMERIC)
    A3s = np.load("Al_A3s.npy")
    for i in xrange(len(A3s)):
        wr.writerow(list(A3s[i]))

with open("Al_A3a.csv", 'wb') as f:
    wr = csv.writer(f, delimiter=' ', \
                    quoting=csv.QUOTE_NONNUMERIC)
    A3a = np.load("Al_A3a.npy")
    for i in xrange(len(A3a)):
        wr.writerow(list(A3a[i]))

