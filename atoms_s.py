'''
This code applies only to fcc crystals with {111} dislocation loops
The script is mostly written in loop coordinates
'''
import numpy as np
import sys, os

#sample = sys.argv[1]
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


'''
In loop coordinates, there should be 3 different planes,
    each of which is a triangular lattice. The ex_p and ey_p vectors
    indicate the orientation of the two sides of the triangle
The lattice orientation is the same for all 3 planes; they only
    differ by a constant translation
'''
## side length of triangle
a1 = a0*np.sqrt(2)/2.
## distance between planes
a2 = a0*np.sqrt(3)/3.
ex_p_th = np.pi/6.
ey_p_th = np.pi/2.
ex_p = a1*np.array([np.cos(ex_p_th), np.sin(ex_p_th), 0.])
ey_p = a1*np.array([np.cos(ey_p_th), np.sin(ey_p_th), 0.])
ez_p = a2*np.array([0., 0., 1.])
## From now on, the coordinates will be in units of ex_p, ey_p, ez_p
'''
Here we assume that the dislocation loop has the configuration of plane A
In this case, starting from the origin, the planes will be (A), C, A, B, C, A, ... 
    with z coordinates (0), 0.5, 1.5, 2.5, 3.5, 4.5, ... in units of a2
The following is the coordinates for one atom on each plane 
'''
orig_A = np.array([0., 0., 1.5])
orig_B = np.array([1./3, 1./3, 2.5])
orig_C = np.array([-1./3, 2./3, 0.5])

## the summation goes to rho = 5R in xy plane, and z from -10R to 10R
##     (both roughly, since the origins will give some offset)
rhomax = 5.*R
zmax = 10.*R
for zind in np.arange(np.ceil(-zmax/(3.*a2)), np.floor(zmax/(3.*a2))):
    if zind == 0.:
        continue
    _xlim = np.floor(rhomax / a1 * 2./np.sqrt(3))
    for xind in np.arange(-_xlim, _xlim+1):
        _ylims = np.roots([1., x, x**2-(rhomax/a1)**2])
        for yind in np.arange(np.ceil(np.min(_ylims)), np.floor(np.max(_ylims))+1):
            x, y, z = (xind*ex_p + yind*ey_p + zind*ez_p)/R

            n = 30
            n1 = n
            integrand = 0.
            for jSum1, th1 in enumerate(np.linspace(0., 2.*np.pi, n1, endpoint=False)):
                coeff = (2 if jSum1%2==0 else 4)*2*np.pi/n1/3
                        # coefficient for th1 integration
                a = -(x*np.cos(th1)+y*np.sin(th1)+1.)/z
                b = -(x*np.cos(th1)+y*np.sin(th1)-1.)/z

                for jSum in np.arange(1, n+1):
                    t = (b+a)/2.+(b-a)*np.cos((2.*jSum-1)*
                                            np.pi/2./n)/2.
                    e = (np.cos(th1)*ex+np.sin(th1)*ey+t*ez)/(1.+t**2)**0.5
                    eloop = np.array([np.cos(th1),np.sin(th1),t])/(1.+t**2)**0.5
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
                    integrand += (coeff*eloop.dot(np.array([x,y,z]))*\
                        np.einsum('ij,i,jk', Ploop, eloop, gloop))
            s = -1./(4.*np.pi**2*R**2*n*abs(z))*integrand
     
