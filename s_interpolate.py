import sys, os
import numpy as np
from scipy.interpolate import griddata

sample = sys.argv[1]
if sample == 'Al':
    from Al_parameters import *
elif sample == 'Cu':
    from Cu_parameters import *

sdata = np.load("data/%s_s_0.npy"%sample)
rank = 1 
while os.path.isfile("data/%s_s_%d.npy"%(sample, rank)):
    sdata = np.vstack((sdata, \
         np.load("data/%s_s_%d.npy"%(sample, rank))))
    rank += 1
print "read %d files in total" % rank

## change to cartesian coordinates
sdata_xyz = np.copy(sdata)
sdata_xyz[:,0] = R*sdata[:,0]*np.cos(sdata[:,1])
sdata_xyz[:,1] = R*sdata[:,0]*np.sin(sdata[:,1])
sdata_xyz[:,2] = R*sdata[:,2]

xmin, xmax = np.min(sdata_xyz[:,0]), np.max(sdata_xyz[:,0])
ymin, ymax = np.min(sdata_xyz[:,1]), np.max(sdata_xyz[:,1])
zmin, zmax = np.min(sdata_xyz[:,2]), np.max(sdata_xyz[:,2])
xx = np.arange(xmin, xmax, a0/3.)
yy = np.arange(ymin, ymax, a0/3.)
zz = np.arange(zmin, zmax, a0/3.)
x_grid, y_grid, z_grid = np.meshgrid(xx,yy,zz,indexing='ij') 

print "begin interpolation..."
sx_interp = griddata(sdata_xyz[:,0:3], sdata_xyz[:,3],\
       np.array([x_grid.ravel(), y_grid.ravel(), 
           z_grid.ravel()]).T).reshape((len(xx), len(yy), len(zz)))
print "sx interpolated"
sy_interp = griddata(sdata_xyz[:,0:3], sdata_xyz[:,4],\
       np.array([x_grid.ravel(), y_grid.ravel(), 
           z_grid.ravel()]).T).reshape((len(xx), len(yy), len(zz)))
print "sy interpolated"
sz_interp = griddata(sdata_xyz[:,0:3], sdata_xyz[:,5],\
       np.array([x_grid.ravel(), y_grid.ravel(), 
           z_grid.ravel()]).T).reshape((len(xx), len(yy), len(zz)))
print "sz interpolated"

np.savez("data/%s_s_interp.npy"%sample, 
        xx=xx, yy=yy, zz=zz,
        sx=sx_interp, sy=sy_interp, sz=sz_interp)
