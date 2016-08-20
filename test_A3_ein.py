import numpy as np
import csv
import os
os.chdir("data")

sample = 'Al'
os.system("rm %s_A3?_ein.npy"%sample)
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

h = 2.*np.pi/a0*np.array([4.,0.,0.])

sdata = np.load("%s_s_th0.npy"%sample)

for qR in np.append(
        -np.logspace(np.log10(5), np.log10(0.1), num=20),
        np.logspace(np.log10(0.1), np.log10(5), num=20)):
    A3s = []
    A3a = []
    for qth in np.array([np.pi/2.]): 
        for qphi in np.array([0.]): 
            q = qR/R*np.array([np.sin(qth)*np.cos(qphi), 
                    np.sin(qth)*np.sin(qphi), np.cos(qth)])
            qloop = np.einsum('ij,j', rot, q)
            K = h+q
            Kloop = np.einsum('ij,j', rot, K)
            sintegrand = 0.
            aintegrand = 0.
            for i in xrange(len(sdata)):
                rho = sdata[i,0]*R
                z = sdata[i,2]*R
                for th in np.linspace(0., 2*np.pi, 20, endpoint=False):
                    if rho == 0.:
                        continue
                    deltaV = 2* rho*0.1*R * 0.1*R * 2*np.pi/20 
                    Ks = Kloop.dot([sdata[i,3]*np.cos(th), \
                            sdata[i,3]*np.sin(th), sdata[i,5]])
                    qr = qloop.dot([rho*np.cos(th), rho*np.sin(th), z])
                    sintegrand += deltaV*(1./Vc*np.linalg.norm(B)*\
                            np.cos(qr)*(np.cos(Ks)-1))
                    aintegrand += -deltaV*(1./Vc*np.linalg.norm(B)*\
                            np.sin(qr)*(np.sin(Ks)-Ks))
            A3s.append([q[0], q[1], q[2], sintegrand])
            A3a.append([q[0], q[1], q[2], aintegrand])

    try:
        A3s = np.vstack((np.load("%s_A3s_ein.npy"%sample), A3s))
        A3a = np.vstack((np.load("%s_A3a_ein.npy"%sample), A3a))
    except IOError:
        A3s = np.array(A3s)
        A3a = np.array(A3a)
    np.save("%s_A3s_ein.npy"%sample, A3s)
    np.save("%s_A3a_ein.npy"%sample, A3a)
    print "done qR =", qR

with open("%s_A3s_ein.csv"%sample, 'wb') as f:
    wr = csv.writer(f, delimiter=' ', \
                    quoting=csv.QUOTE_NONNUMERIC)
    A3s = np.load("%s_A3s_ein.npy"%sample)
    for i in xrange(len(A3s)):
        wr.writerow(list(A3s[i]))

with open("%s_A3a_ein.csv"%sample, 'wb') as f:
    wr = csv.writer(f, delimiter=' ', \
                    quoting=csv.QUOTE_NONNUMERIC)
    A3a = np.load("%s_A3a_ein.npy"%sample)
    for i in xrange(len(A3a)):
        wr.writerow(list(A3a[i]))

