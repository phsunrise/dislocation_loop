import numpy as np

## this function calculates the integral of e^(iq.r) over the Wigner-Seitz
## cell of BCC crystal. Note that q should be in crystal coordinates
## TODO: still haven't considered the case qx^2=qy^2!=0
def WSint(q, a0):
    qx, qy, qz = q
    if np.abs(qx)+np.abs(qy) < 1.e-16:
        return WSint_z(qz, a0)
    elif np.abs(qy)+np.abs(qz) < 1.e-16:
        return WSint_z(qx, a0)
    elif np.abs(qz)+np.abs(qx) < 1.e-16:
        return WSint_z(qy, a0)
    else: ## at most one component is zero 
        XY = qx**2-qy**2
        YZ = qy**2-qz**2
        ZX = qz**2-qx**2
        Sx2 = np.sin(qx*a0/2.)
        Cx2 = np.cos(qx*a0/2.)
        Sx4 = np.sin(qx*a0/4.)
        Cx4 = np.cos(qx*a0/4.)
        Sy2 = np.sin(qy*a0/2.)
        Cy2 = np.cos(qy*a0/2.)
        Sy4 = np.sin(qy*a0/4.)
        Cy4 = np.cos(qy*a0/4.)
        Sz2 = np.sin(qz*a0/2.)
        Cz2 = np.cos(qz*a0/2.)
        Sz4 = np.sin(qz*a0/4.)
        Cz4 = np.cos(qz*a0/4.)

        A1 = Cz2*(qx*Sx4*YZ+qy*Sy4*ZX) + qz*Sz2*(YZ*Cx4+ZX*Cy4)
        A2 = Cx2*(qy*Sy4*ZX+qz*Sz4*XY) + qx*Sx2*(ZX*Cy4+XY*Cz4)
        A3 = Cy2*(qz*Sz4*XY+qx*Sx4*YZ) + qy*Sy2*(XY*Cz4+YZ*Cx4)
        print XY, YZ, ZX, A1, A2, A3
        return -8./(XY*YZ*ZX)*(A1+A2+A3)
    

## special case where at least two components of q is nonzero
## The parameter qz is the third component which may or may not be zero
def WSint_z(qz, a0):
    if np.abs(qz) < 1.e-16: # all components of q are 0
        return a0**3/2.
    else:
        return (8./qz**4)*(qz**2*a0/4.*(1.-np.cos(qz*a0/2.))-\
                    qz*np.sin(qz*a0/2.)*(1.-(qz*a0)**2/32.)+\
                    2.*qz*np.sin(qz*a0/4.))

## use a sphere to approximate the cell
def sphint(qnorm, a0):
    if qnorm < 1.e-16:
        return a0**3/2.
    else:
        r = (3./8./np.pi)**(1./3.)*a0
        return 4.*np.pi/qnorm**3*(-qnorm*r*np.cos(qnorm*r)+np.sin(qnorm*r))

## this function calculates the amplitude integral from r=r0 to r=inf
## under the isotropic continuum approximation
def ampint(qvec, Kvec, r0, bvec, R, C12, C44, a0):
    qvec = np.array(qvec)
    Kvec = np.array(Kvec)
    bvec = np.array(bvec)
    
    nu = C12/2./(C12+C44)
    Vc = sphint(np.linalg.norm(qvec), a0) 
    A = np.pi*R**2/(4.*(1.-nu))/Vc
    ## first need to find the new coordinates where qvec'=(0,0,q)
    if np.linalg.norm(qvec[0:2]) < 1.e-16:
            ## if qvec is already in z direction
        q = qvec[2]
        qr = q*r0 
        Sin = np.sin(qr)
        Cos = np.cos(qr)
        return -A/q*(2.*(2.*nu-1.)*(bvec.dot(Kvec))*Sin/qr\
          -3.*((bvec[0]*Kvec[0]+bvec[1]*Kvec[1])*(2.*Sin-2.*qr*Cos)/(qr**3)\
               +2.*bvec[2]*Kvec[2]*(2.*qr*Cos+(qr**2-2.)*Sin)/(qr**3))) # minus sign comes from integration to inf

    else:
        ## variables ending in "p" indicates the new (primed) coordinate system
        q = np.linalg.norm(qvec)
        ezp = qvec/q
        eyp = np.cross(ezp, [0.,0.,1.])
        eyp = eyp/np.linalg.norm(eyp)
        exp = np.cross(eyp, ezp)
        rot = np.array([exp, eyp, ezp]) # rotation matrix from unprimed to primed system 
        Sina, _, Cosa = rot.dot([0.,0.,1.])
        qr = q*r0 
        Sin = np.sin(qr)
        Cos = np.cos(qr)
        bvecp = rot.dot(bvec)
        Kvecp = rot.dot(Kvec)
        
        B = 2.*(2.*nu-1.)*(bvec.dot(Kvec)*Cosa+bvecp[2]*Kvec[2]-bvec[2]*Kvecp[2])*Sin/qr
        C = 3.*Cosa*((bvecp[0]*Kvecp[0]+bvecp[1]*Kvecp[1])*(2.*Sin-2.*qr*Cos)/(qr**3)+2.*bvecp[2]*Kvecp[2]*(2.*qr*Cos+(qr**2-2.)*Sin)/(qr**3))
        D = 3.*(bvecp[0]*Kvecp[2]+bvecp[2]*Kvecp[0])*Sina*(2.*Sin-2.*qr*Cos)/(qr**3)
        return -A/q*(B-C-D) # minus sign comes from integration to inf

if __name__=="__main__":
    from W_parameters import R, a0, B, C12, C44
    ## check: Kvec=(K,K,0), qvec=(q,q,0), bvec=(b,b,b)
    q = 0.3
    K = 4.*np.pi/a0+q
    b = B[0]
    r0 = 400.
    nu = C12/2./(C12+C44)
    Vc = sphint(np.linalg.norm([q,q,0.]), a0) 
    print ampint(qvec=[q,q,0.], Kvec=[K,K,0.], \
                 bvec=[b,b,b], R=R, r0=r0, C12=C12, C44=C44, a0=a0)

    qr = np.sqrt(2.)*q*r0
    Sin = np.sin(qr)
    Cos = np.cos(qr)
    print np.pi*b*R**2*K/(2.*(1.-nu))*((2.*nu-1.)/q*Sin/qr\
            +3./q*(Sin-qr*Cos)/(qr**3))/Vc

    import matplotlib.pyplot as plt
    qz_list = np.linspace(-0.5, 0.5, 101)
    Vc = np.zeros_like(qz_list)
    for i_qz, qz in enumerate(qz_list):
        Vc[i_qz] = sphint(np.abs(qz), a0)

    plt.plot(qz_list, Vc/(a0**3/2.), 'ro')
    plt.show()
