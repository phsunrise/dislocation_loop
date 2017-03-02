import numpy as np
import scipy.optimize 
from settings import basedir
fitdir = basedir+"fit/"
import matplotlib.pyplot as plt

sample = "0.2dpa"
fdata = np.load(fitdir+"%s_q4I_2.npz"%sample)
q = fdata['q']
q4I = fdata['q4I']

## plot every 8 points
n = 8
length = len(q)
q = q.reshape(length/n, n)
q = np.mean(q, axis=1)
q4I = q4I.reshape(length/n, n)
q4I_err = np.std(q4I, axis=1)
q4I = np.mean(q4I, axis=1)

fig, ax = plt.subplots(1,1, figsize=(10,5))
ax.errorbar(q, q4I, yerr=q4I_err, \
            ls='', color='r', marker='o')
ax.set_xlim(-0.5, 0.5)
ax.set_ylim(0., 2.e-6)

#print q-(np.linspace(-0.25,0.25,128, endpoint=False)+0.5/128/2.)*\
#        (np.sqrt(2)*2.*np.pi/3.1648)

#f = np.load(fitdir+"R5_int.npz")
#q = f['q']
#print q-(np.linspace(-0.25,0.25,128, endpoint=False)+0.5/128/2.)*\
#        (np.sqrt(2)*2.*np.pi/3.16)
q4I_th = []
c0 = []
for looptype in ['int', 'vac']:
    for R in np.arange(5, 55, 5):
        _data = np.load(fitdir+"R%d_%s.npz"%(R, looptype))['q4I']
        q4I_th.append(_data)
        c0.append(2.e-7/np.max(_data))
q4I_th = np.array(q4I_th)
c0 = np.array(c0)
## the curve fitting procedure is as follows:
## for an array c = [c5_int, c10_int, ..., c50_int, c5_vac, ..., c50_vac],
## calculate the theoretical q4I values, then calculate the diffrence from
## actual data, sum the squares weighted by 1/err^2
fitrg = (-0.4, 0.4)
indmin, indmax = np.searchsorted(q, fitrg)
q = q[indmin:indmax]
q4I = q4I[indmin:indmax]
q4I_err = q4I_err[indmin:indmax]
q4I_th = q4I_th[:, indmin:indmax]
def difsq(c):
    q4I_th_sum = c.dot(q4I_th)
    diff = q4I-q4I_th_sum
    return np.sum(diff**2/q4I_err**2)

res = scipy.optimize.minimize(difsq, c0)
print res    

ax.plot(q, res.x.dot(q4I_th), 'b-')
ax.plot(q, c0.dot(q4I_th), 'g-')

fig, ax = plt.subplots(1,1)
ax.plot(q, res.x.dot(q4I_th), 'b-')
plt.show()
