import numpy as np
import os, sys
from scipy import special
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

sample = 'Al'
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

P = np.zeros((3, 3))
for i in xrange(3):
    for j in xrange(3):
        P[i,j] = ((C12*F.dot(B) + d*F[i]*B[i]) * \
                    (1 if i==j else 0) + \
                    C44*(F[i]*B[j]+F[j]*B[i])
            )
Ploop =  np.einsum('ij,kl,jl', rot, rot, P)

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

I_data = []
rank = 0
while os.path.isfile("data/%s_A3s_%d.npy"%(sample,rank)):
    A3s = np.load("data/%s_A3s_%d.npy"%(sample,rank))
    A3a = np.load("data/%s_A3a_%d.npy"%(sample,rank))
    for i in xrange(len(A3s)):
        q = A3s[i, 0:3] 
        if q.dot(q) <= (q.dot(ez))**2:
            Q = 1.e-6
        else:
            Q = np.sqrt(q.dot(q) - (q.dot(ez))**2)
        G = g / (q.dot(q))
        
        A = np.linalg.norm(B)*np.pi*R**2/Vc*\
                2.*special.jv(1, Q*R)/(Q*R)) # A1s
        A += -np.einsum("i,ij,jk,k",h,G,P,q)/Vc*\
                    2.*special.jv(1, Q*R)/(Q*R)) # A2a
        A += -np.einsum("i,ij,jk,k",h,G,P,q)/Vc*\
                    2.*special.jv(1, Q*R)/(Q*R)) # A2s
        A += A3s[i,3] + A3a[i,3]
        I_data.append(q[0], q[1], q[2], A**2)

    print "done rank %d" % rank
    rank += 1

np.save("data/%s_intensities.npy"%sample, I_data)
