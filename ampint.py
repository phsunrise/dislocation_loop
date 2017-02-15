import numpy as np

## this function calculates the amplitude integral from r=r0 to r=inf
## under the assumption of isotropic continuum. 
## TODO: currently only works for qvec = (q,q,0), Kvec = (k,k,0), and bvec = (b,b,b)
def ampint(q, K, r0, b, R, Vc, C12, C44):
    nu = C12/2./(C12+C44)
    A = np.pi*b*R**2*K/(2.*(1.-nu)*Vc)
    B = (2.*nu-1)*np.sin(np.sqrt(2.)*r0*q)/(np.sqrt(2.)*r0*q**2)
    C = 1.5/(np.sqrt(2.))*(np.sin(np.sqrt(2.)*r0*q)- np.sqrt(2.)*r0*q*np.cos(np.sqrt(2.)*r0*q))/(q**4*r0**3)
    return A*(B+C)
