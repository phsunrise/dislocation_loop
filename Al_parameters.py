import numpy as np

crystaltype = 'FCC'

a0 = 4.040 # Angstrom
R = 10.
Vc = 16.61

ex = np.array([1.,1.,-2.])
ey = np.array([-1.,1.,0.])
ez = np.array([1.,1.,1.])
ex = ex / np.linalg.norm(ex)
ey = ey / np.linalg.norm(ey)
ez = ez / np.linalg.norm(ez)
rot = np.array([ex, ey, ez])
    # rotation matrix from crystal coordinates to 
    # dislocation loop coordinates

F = np.pi*R**2 * ez
B = a0/3.*np.array([1.,1.,1.])

#Al elastic constants, taken from https://arxiv.org/pdf/1605.09237.pdf
C11 = 107.3 # GPa
C12 = 60.08
C44 = 28.3
#given in Ben's email
C11 = 11.63
C12 = 6.48
C44 = 3.09

d = C11 - C12 - 2.*C44

h = 2.*np.pi/a0*np.array([4.,0.,0.])
eq = np.array([1.,0.,0.])
eq = eq / np.linalg.norm(eq)

## all four orientations are equivalent for h and q
## parallel to <100>
orientations = [np.diag([1.,1.,1.])]

## P tensor
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
Here we assume that the interstitial loop has the configuration of plane A
In this case, starting from the origin, the planes will be (A), C, A, B, C, A, ... 
    with z coordinates (0), 0.5, 1.5, 2.5, 3.5, 4.5, ... in units of a2
The following is the coordinates for one atom on each plane 
'''
orig_A = np.array([0., 0., 1.5])
orig_B = np.array([-1./3, -1./3, 2.5])
orig_C = np.array([1./3, 1./3, 0.5])
