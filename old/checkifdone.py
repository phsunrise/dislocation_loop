import numpy as np

nproc = 16

rho_array = np.linspace(0., 5., 51)
th_array = np.linspace(0., 2.*np.pi, 20, endpoint=False)
z_array = np.linspace(0.1, 10., 100)
for rank in xrange(nproc):
    data = np.load("data/Al_s_%d.npy"%rank)
    temp = 50/nproc*nproc+rank
    if temp <= 50:
        rho = rho_array[temp]
    else:
        rho = rho_array[temp-nproc]

    if np.linalg.norm(data[-1, 0:3] - np.array(
            [rho, th_array[-1], z_array[-1]])) > 1.e-6:
        print "rank %d not done!, last entry:\n\t"%rank,
        print data[-1]
        print "should be:\n\t",
        print np.array([rho, th_array[-1], z_array[-1]])
    else:
        print "rank %d done!"%rank
