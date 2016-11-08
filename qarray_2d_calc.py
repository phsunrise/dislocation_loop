import numpy as np

sample = 'W'
if sample == 'W':
    from W_parameters import *

E = 15.5 # keV
wavelength = 12.39842 / E # Angstrom
k = 2*np.pi/wavelength

h = 2*np.pi/a0*np.array([2.,0.,0.])
ehp = np.array([0.,0.,1.]) # this is a direction perpendicular to h
ki = -0.5*h + ehp*np.sqrt(k**2-h.dot(h)/4.)
kf0 = 0.5*h + ehp*np.sqrt(k**2-h.dot(h)/4.) # kf for which q=0

# xray coordinates
ez_b = ki / np.linalg.norm(ki)
ey_b = np.cross(ez_b, h)
ey_b = ey_b / np.linalg.norm(ey_b)
ex_b = np.cross(ey_b, ez_b)

tth0 = 2.*np.arcsin(np.linalg.norm(h)/2./k) # 2theta for Bragg peak
tth_array = np.linspace(-5., 5., 21)/180.*np.pi + tth0
tth_low = tth_array[0]*1.5-tth_array[1]*0.5
tth_high = tth_array[-1]*1.5-tth_array[-2]*0.5
chi_array = np.linspace(-5., 5., 21)/180.*np.pi
chi_low = chi_array[0]*1.5-chi_array[1]*0.5
chi_high = chi_array[-1]*1.5-chi_array[-2]*0.5
#chi, tth = np.meshgrid(chi_array, tth_array, indexing='xy')
#print "chi=", chi
#print "tth=", tth

q_array = np.zeros((len(tth_array), len(chi_array), 3))
K_array = np.zeros((len(tth_array), len(chi_array), 3))
q_array_norm = np.zeros((len(tth_array), len(chi_array)))
for i_tth, tth in enumerate(tth_array):
    for i_chi, chi in enumerate(chi_array):
        kf = k*(ex_b*np.sin(tth)*np.cos(chi)+\
                ey_b*np.sin(tth)*np.sin(chi)+\
                ez_b*np.cos(tth))
        K_array[i_tth, i_chi] = kf - ki
        q_array[i_tth, i_chi] = kf - ki - h 
        q_array_norm[i_tth, i_chi] = np.linalg.norm(kf - ki - h)
    print "done %.2f%%" % ((i_tth+1)*1./len(tth_array)*100)

np.savez("qarray_2d.npz", q_array=q_array, K_array=K_array, h=h)

import matplotlib.pyplot as plt
plt.imshow(q_array_norm, origin='lower', interpolation='nearest', \
           extent=(chi_low, chi_high, tth_low, tth_high), aspect='auto')
plt.xlabel("chi")
plt.ylabel("tth")
plt.colorbar()
plt.show()
