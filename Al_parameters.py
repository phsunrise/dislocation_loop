import numpy as np

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

orientations = [np.diag([1.,1.,1.]),
                np.diag([1.,-1.,-1.]),
                np.diag([-1.,1.,-1.]),
                np.diag([-1.,-1.,1.])]
