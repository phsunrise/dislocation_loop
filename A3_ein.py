import numpy as np
import csv
import os, sys
os.chdir("data")

sample = sys.argv[1]
os.system("rm %s_A3?_ein.*" % sample)
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

for qR in np.append(
        -np.logspace(np.log10(5.), np.log10(0.1), num=40),
        np.logspace(np.log10(0.1), np.log10(5.), num=40)):
    A3s = []
    A3a = []
    q = qR/R*eq
    qloop = np.einsum('ij,j', rot, q)
    K = h+q
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

