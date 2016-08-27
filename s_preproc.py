import numpy as np
import os

# first remove existing files
os.system("rm preproc/s_????.npy")

r_array = np.linspace(0., 10., 101)
th_array = np.linspace(0., 2*np.pi, 20, endpoint=False)
z_array = np.linspace(0.1, 20., 200)
th_grid, r_grid, z_grid = np.meshgrid(r_array, th_array, z_array, indexing='ij')
thrz_list = np.array([th_grid, r_grid, z_grid]).ravel()
thrz_list = thrz_list.reshape(len(thrz_list)/3, 3)
thrz_list_proc = []
for th, r, z in thrz_list:
    if r == 0. and th != 0.:
        continue
    thrz_list_proc.append([th, r, z])
nproc = raw_input("Total length %d, number of processors? "%len(thrz_list_proc))
nproc = int(nproc)

sublistlen = int(np.ceil(len(thrz_list_proc)*1./nproc))
for i in xrange(nproc):
    np.save("preproc/s_%04d.npy"%i, \
        np.array(thrz_list_proc[i*sublistlen:(i+1)*sublistlen]))
    print "file %d saved"%i
