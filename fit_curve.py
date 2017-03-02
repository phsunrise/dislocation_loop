import numpy as np
import scipy.optimize 
from settings import basedir, Rlist
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

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1,1,1)
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
    for R in Rlist: 
        _data = np.load(fitdir+"R%d_%s.npz"%(R, looptype))['q4I']
        q4I_th.append(_data)
        c0.append(2.e-7/np.max(_data))
q4I_th = np.array(q4I_th)
c0 = np.array(c0)
#c0 = np.array(c0)**0.5
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
def difsq(ind, *c):
    return np.abs(c).dot(q4I_th[:,ind])

bounds = (1.e-11, 1.e-6)
res = scipy.optimize.curve_fit(difsq, np.arange(len(q4I)), q4I, \
                p0=c0, sigma=q4I_err)#, bounds=bounds)
print res[0]

ax.plot(q, np.abs(res[0]).dot(q4I_th), 'b-')
ax.plot(q, np.abs(c0).dot(q4I_th), 'g-')

plt.show()
