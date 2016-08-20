import numpy as np
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

qR_array = []
A1s = []
A2a = []
A2s = []
A3s = []
A3a = []
for qR in np.logspace(np.log10(0.1), np.log10(5.), num=50):
    qR_array.append(qR)
    q = qR / R * e
    if q.dot(q) <= (q.dot(ez))**2:
        Q = 1.e-6
    else:
        Q = np.sqrt(q.dot(q) - (q.dot(ez))**2)
    G = g / (q.dot(q))
    A1s.append(np.linalg.norm(B)*np.pi*R**2/Vc*\
                    2.*special.jv(1, Q*R)/(Q*R))
    A2a.append(-np.einsum("i,ij,jk,k",h,G,P,q)/Vc*\
                    2.*special.jv(1, Q*R)/(Q*R))
    A2s.append(-np.einsum("i,ij,jk,k",q,G,P,q)/Vc*\
                    2.*special.jv(1, Q*R)/(Q*R))
qR_array = np.array(qR_array)
A1s = np.array(A1s)
A2a = np.array(A2a)
A2s = np.array(A2s)

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax1.plot(qR_array, A1s, 'r--')
ax1.plot(qR_array, -A1s, 'r-', label='A1s')
ax1.plot(qR_array, A2a, 'g--')
ax1.plot(qR_array, -A2a, 'g-', label='A2a')
ax1.plot(qR_array, A2s, 'b--')
ax1.plot(qR_array, -A2s, 'b-', label='A2s')

## now get and plot A3 data
A3s = np.load("data/%s_A3s_ein.npy"%sample)
A3a = np.load("data/%s_A3a_ein.npy"%sample)
ax1.plot(A3s[:,0]*R, A3s[:,3], 'c--')
ax1.plot(A3s[:,0]*R, -A3s[:,3], 'c-', label='A3s')
ax1.plot(A3a[:,0]*R, A3a[:,3], 'm--')
ax1.plot(A3a[:,0]*R, -A3a[:,3], 'm-', label='A3a')

ax1.set_xscale("log", nonposx='clip')
ax1.set_yscale("log", nonposy='clip')
ax1.set_xlim(0.1, 5)
ax1.set_ylim(1.e1, 1.e5)
ax1.set_xlabel(r"$qR$")
ax1.set_ylabel(r"$A(\vec{K})$")
ax1.legend(loc='upper right')

## plot intensity
ax2 = fig.add_subplot(1,2,2)
f_A1s = interp1d(qR_array, A1s, kind='cubic')
f_A2a = interp1d(qR_array, A2a, kind='cubic')
f_A2s = interp1d(qR_array, A2s, kind='cubic')
f_A3a = interp1d(A3a[:,0]*R, A3a[:,3], kind='cubic')
f_A3s = interp1d(A3s[:,0]*R, A3s[:,3], kind='cubic')
ax2.plot(qR_array, \
         (qR_array/R)**4*(f_A1s(qR_array)+
            f_A2a(qR_array)+f_A2s(qR_array)+
            f_A3a(qR_array)+f_A3s(qR_array))**2,
         'b-', label="pos")
ax2.plot(qR_array, \
         (qR_array/R)**4*(f_A1s(qR_array)-
            f_A2a(qR_array)+f_A2s(qR_array)+
            f_A3a(-qR_array)+f_A3s(-qR_array))**2,
         'r-', label="neg")

#ax2.set_xscale("log", nonposx='clip')
#ax2.set_yscale("log", nonposy='clip')
ax2.set_xlim(0.1, 5)
ax2.set_ylim(1.e1, 10000)
ax2.set_xlabel(r"$qR$")
ax2.set_ylabel(r"$q^2|A(\vec{K})|^2$ $\mathrm{\AA^{-2}}$")
ax2.legend(loc='lower left')

fig.tight_layout()
fig.savefig("plots/%s_amplitudes.pdf"%sample)
plt.show()
