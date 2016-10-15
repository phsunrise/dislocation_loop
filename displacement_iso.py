import numpy as np
from scipy import special
import sys, os
from W_parameters import a0, B, C12, C44

def solang(x, y, z, R):
    ''' 
    this function returns the solid angle of a disk of radius R,
    in +z direction, as seen from (x,y,z)
    Reference: Paxton, Solid angle calculation for a solid disk, 
               Rev. Sci. Ins. 30, 254 (1959)
    '''
    r0 = np.sqrt(x*x+y*y)
    Rmax = np.sqrt(z**2 + (r0+R)**2)
    R1 = np.sqrt(z**2 + (r0-R)**2)
    k = 1. - (R1/Rmax)**2
    kp = (R1/Rmax)**2
    xi = np.arcsin(z/R1)
    gamma0 = 2./np.pi*(
            special.ellipe(k)*special.ellipkinc(xi, kp)
           +special.ellipk(k)*special.ellipeinc(xi, kp)
           -special.ellipk(k)*special.ellipkinc(xi, kp))

    if r0 < R:
        return 2.*np.pi-2.*z/Rmax*special.ellipk(k) - np.pi*gamma0
    elif r0 == R:
        return np.pi - 2.*z/Rmax*special.ellipk(k)
    else: # r0 > R
        return -2.*z/Rmax*special.ellipk(k) + np.pi*gamma0

def disp(b, x, y, z, R, C12, C44):
    nu = C12/2./(C12+C44)

    r = np.sqrt(x*x+y*y)
    R1 = np.sqrt(z**2+(r-R)**2)
    E = special.ellipe(-4.*R*r/R1**2)
    K = special.ellipk(-4.*R*r/R1**2)
    #print "E=%f, K=%f, nu=%f" % (E, K, nu)
    # first term, containing the solid angle
    term1 = -b/(4.*np.pi)*solang(x,y,z,R)

    # second term
    term2 = 1./(4.*np.pi)*(\
              np.cross(b, [-y, x, 0])*2./r**2*\
              (-R1*E+(r**2+R**2+z**2)/R1*K))

    # third term
    bx, by, bz = b
    term3 = np.array([0.,0.,0.])
    term3[0] = (2*(-(((R**2*np.sqrt(x**2 + y**2) - 2*R*(x**2 + y**2) + \
                    np.sqrt(x**2 + y**2)*(x**2 + y**2 + z**2))*\
                    ((R**2 - x**2 - y**2)*(-(by*x*y*\
                    (-2*R**2 + x**2 + y**2)) + bx*(R**2*(x - y)*(x + y)\
                     + y**2*(x**2 + y**2)))*z + (bx*(2*R**2*x**2 + \
                     x**4 - (2*R**2 + x**2)*y**2 - 2*y**4) + \
                     by*x*y*(4*R**2 + 3*(x**2 + y**2)))*z**3 + \
                     (2*by*x*y + bx*(x - y)*(x + y))*z**5 + \
                     bz*x*(x**2 + y**2)*((-R**2 + x**2 + y**2)**2 + \
                     (R**2 + x**2 + y**2)*z**2))*E)/(x**2 + y**2)**2.5)\
                      + (((-R**2 + x**2 + y**2)**2 + 2*(R**2 + x**2 + y**2)\
                      *z**2 + z**4)*(bz*x*(x**2 + y**2)*(R**2 + x**2 + \
                      y**2) + z*(by*x*y*(2*R**2 + x**2 + y**2 + 2*z**2) + \
                      bx*(R**2*(x - y)*(x + y) - y**2*(x**2 + y**2) + \
                      (x - y)*(x + y)*z**2)))*K)/(x**2 + y**2)**2))\
                      /((R**2 + x**2 + y**2 - 2*R*np.sqrt(x**2 + y**2) + \
                      z**2)**1.5*(R**2 + x**2 + y**2 + 2*R*np.sqrt(x**2 + \
                      y**2) + z**2))
    term3[1] = (2*(-(((R**2*np.sqrt(x**2 + y**2) - 2*R*(x**2 + y**2) + \
                    np.sqrt(x**2 + y**2)*(x**2 + y**2 + z**2))*((-R**2 + \
                    x**2 + y**2)*(by*(R - x)*x**2*(R + x) - by*(R**2 + x**2)*y**2 + \
                    bx*x*y*(-2*R**2 + x**2 + y**2))*z + (by*(-2*x**2*(R**2 + x**2) + \
                    (2*R**2 - x**2)*y**2 + y**4) + bx*x*y*(4*R**2 + \
                    3*(x**2 + y**2)))*z**3 + (2*bx*x*y + by*(-x**2 + y**2))*z**5 + \
                    bz*y*(x**2 + y**2)*((-R**2 + x**2 + y**2)**2 + \
                    (R**2 + x**2 + y**2)*z**2))*E)/(x**2 + y**2)**2.5) + (((-R**2 + \
                    x**2 + y**2)**2 + 2*(R**2 + x**2 + y**2)*z**2 + \
                    z**4)*(bz*y*(x**2 + y**2)*(R**2 + x**2 + y**2) + \
                    z*(bx*x*y*(2*R**2 + x**2 + y**2 + 2*z**2) - by*(x**4 + \
                    R**2*(x - y)*(x + y) - y**2*z**2 + x**2*(y**2 + \
                    z**2))))*K)/(x**2 + y**2)**2))/((R**2 + x**2 + \
                    y**2 - 2*R*np.sqrt(x**2 + y**2) + z**2)**1.5*(R**2 + x**2 + \
                    y**2 + 2*R*np.sqrt(x**2 + y**2) + z**2))
    term3[2] = (2*(((bx*x + by*y)*(-R**2 + x**2 + y**2)**2 - bz*(x**2 + \
                    y**2)*(-R**2 + x**2 + y**2)*z + 3*(bx*x + by*y)*(R**2 + \
                    x**2 + y**2)*z**2 - bz*(x**2 + y**2)*z**3 + 2*(bx*x + \
                    by*y)*z**4)*E - (R**2 + x**2 + y**2 + 2*R*np.sqrt(x**2 + \
                    y**2) + z**2)*((bx*x + by*y)*(R**2 + x**2 + y**2) - \
                    bz*(x**2 + y**2)*z + 2*(bx*x + by*y)*z**2)*K))/((x**2 + \
                    y**2)*np.sqrt(R**2 + x**2 + y**2 - 2*R*np.sqrt(x**2 + \
                    y**2) + z**2)*(R**2 + x**2 + y**2 + 2*R*np.sqrt(x**2 + \
                    y**2) + z**2))
    term3 /= (8.*np.pi*(1.-nu))

    #print "B =", b
    #print "term1 =", term1
    #print "term2 =", term2
    #print "term3 =", term3
    #print "term3*(8pi(1-nu)) =", term3*(8.*np.pi*(1.-nu))
    #print (term1+term2+term3)/a0
    return term1-term2-term3
    
if __name__=='__main__':
   x, y, z, R = sys.argv[1:] 
   x = float(x)*a0
   y = float(y)*a0
   z = float(z)*a0
   R = float(R)*a0
   disp(B, x, y, z, R, C12, C44)
