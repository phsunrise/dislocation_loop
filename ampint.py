import numpy as np

## this function calculates the amplitude integral from r=r0 to r=inf
## under the isotropic continuum approximation
def ampint(qvec, Kvec, r0, bvec, R, C12, C44, a0):
    qvec = np.array(qvec)
    Kvec = np.array(Kvec)
    bvec = np.array(bvec)
    
    nu = C12/2./(C12+C44)
    A = np.pi*R**2/(4.*(1.-nu))
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
    print ampint(qvec=[q,q,0.], Kvec=[K,K,0.], \
                 bvec=[b,b,b], R=R, r0=r0, C12=C12, C44=C44, a0=a0)

    qr = np.sqrt(2.)*q*r0
    Sin = np.sin(qr)
    Cos = np.cos(qr)
    print np.pi*b*R**2*K/(2.*(1.-nu))*((2.*nu-1.)/q*Sin/qr\
            +3./q*(Sin-qr*Cos)/(qr**3))
