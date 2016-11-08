import numpy as np

crystaltype = 'BCC'

a0 = 3.160 # Angstrom
R = 20. 
D = 4.*R
Vc = 15.78

ex = np.array([-1.,1.,0.])
ey = np.array([0.,0.,1.])
ez = np.array([1.,1.,0.])
ex = ex / np.linalg.norm(ex)
ey = ey / np.linalg.norm(ey)
ez = ez / np.linalg.norm(ez)
rot = np.array([ex, ey, ez])
    # rotation matrix from crystal coordinates to 
    # dislocation loop coordinates

F = np.pi*R**2 * ez
B = a0/2.*np.array([1.,1.,1.])
Bloop = np.einsum('ij,j', rot, B)

## W elastic constants at 300K, taken from Featherston (1963)
C11 = 5.2327 
C12 = 2.0453
C44 = 1.6072
d = C11 - C12 - 2.*C44

h = 2.*np.pi/a0 * np.array([2.,0.,0.])
eq = np.array([1.,0.,0.])
eq = eq / np.linalg.norm(eq)

## four possibilities for loop orientations
## in order to accelerate the calculation, we rotate h and eq,
## but not ex, ey, ez, so that the s field results can still be used

#orientations_pm = [np.diag([1.,1.,1.]),
#                   np.diag([1.,-1.,-1.]),
#                   np.diag([-1.,1.,-1.]),
#                   np.diag([-1.,-1.,1.])]
orientations_pm = [np.diag([1.,1.,1.])]
orientations_perm = [np.diag([1.,1.,1.]),
                     np.array([[0.,1.,0.],[0.,0.,1.],[1.,0.,0.]]),
                     np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,0.]])]
orientations = [np.dot(mat1, mat2) for mat1 in orientations_pm \
                        for mat2 in orientations_perm]

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
In loop coordinates, there should be 2 different planes,
    each of which is a triangular lattice. The ex_p and ey_p vectors
    indicate the orientation of the two sides of the triangle
The lattice orientation is the same for all 2 planes; they only
    differ by a constant translation
'''
## side length on x, y direction of triangle
a1 = a0*np.sqrt(3)/2.
## distance between planes
a2 = a0*np.sqrt(2)/2.
ex_p_th = -np.arcsin(1./np.sqrt(3.))
ey_p_th = np.arcsin(1./np.sqrt(3.))
ex_p = a1*np.array([np.cos(ex_p_th), np.sin(ex_p_th), 0.])
ey_p = a1*np.array([np.cos(ey_p_th), np.sin(ey_p_th), 0.])
ez_p = a2*np.array([0., 0., 1.])
## From now on, the coordinates will be in units of ex_p, ey_p, ez_p
'''
Here we assume that the dislocation loop has the configuration of plane A
In this case, starting from the origin, the planes will be (A), B, A, B, ...
    with z coordinates (0), 0.5, 1.5, 2.5, ... in units of a2
The following is the coordinates for one atom on each plane 
'''
orig_A = np.array([0., 0., 1.5])
orig_B = np.array([0.5, 0.5, 0.5])
