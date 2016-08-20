import numpy as np

a0 = 3.6147 # Angstrom
R = 10.
Vc = 11.82

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

## Cu elastic constants, taken from Ohr (1974)
C11 = 1.69
C12 = 1.22
C44 = 0.755
d = C11 - C12 - 2.*C44

h = 2.*np.pi/a0 * np.array([2.,2.,2.])
e = np.array([1.,1.,1.])
