import numpy as np
import os, sys

sample = 'Cu'
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

sdata = np.load("data/%s_s_0.npy"%sample)
rank = 1
while os.path.isfile("data/%s_s_%d.npy"%(sample,rank)):
    sdata = np.vstack((sdata, 
                    np.load("data/%s_s_%d.npy"%(sample,rank))))
    rank += 1
print "imported %d files" % rank

## parallelism
if len(sys.argv) == 3:
    nproc = int(sys.argv[1])
    rank = int(sys.argv[2])
else:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()

for i_qR, qR in enumerate(np.logspace(np.log10(0.1), np.log10(5), 20)):
    if i_qR % nproc != rank:
        continue
    A3s = []
    A3a = []
    for qth in np.linspace(0., np.pi, 21):
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
        A3s = np.vstack((np.load("data/%s_A3s_%d.npy"%(sample,rank)), A3s))
        A3a = np.vstack((np.load("data/%s_A3a_%d.npy"%(sample,rank)), A3a))
    except IOError:
        A3s = np.array(A3s)
        A3a = np.array(A3a)
    np.save("data/%s_A3s_%d.npy"%(sample,rank), A3s)
    np.save("data/%s_A3a_%d.npy"%(sample,rank), A3a)

#with open("data/%s_A3s.csv"%sample, 'wb') as f:
#    wr = csv.writer(f, delimiter=' ', \
#                    quoting=csv.QUOTE_NONNUMERIC)
#    A3s = np.load("data/%s_A3s.npy"%sample)
#    for i in xrange(len(A3s)):
#        wr.writerow(list(A3s[i]))
#
#with open("data/%s_A3a.csv"%sample, 'wb') as f:
#    wr = csv.writer(f, delimiter=' ', \
#                    quoting=csv.QUOTE_NONNUMERIC)
#    A3a = np.load("data/%s_A3a.npy"%sample)
#    for i in xrange(len(A3a)):
#        wr.writerow(list(A3a[i]))
