import sys
import numpy as np

model = np.loadtxt(sys.argv[1])

depth = model[:,0]
vpMod = model[:,1]

nx = int(sys.argv[2])
ny = int(sys.argv[3])
nz = int(sys.argv[4])

vp = np.zeros((nz,nx,ny))
vs = np.zeros((nz,nx,ny))
rho = np.zeros((nz,nx,ny))

vp[:int(depth[0]),:,:] = np.ones((int(depth[0]),nx,ny)) * vpMod[0]
vs[:int(depth[0]),:,:] = np.zeros((int(depth[0]),nx,ny))
rho[:int(depth[0]),:,:] = np.ones((int(depth[0]),nx,ny)) * 1000.0

for i in range(1,len(depth)):
    var = slice(int(depth[i-1]),int(depth[i]))
    vp[var,:,:] = np.ones((int(depth[i]-depth[i-1]),nx,ny)) * vpMod[i]
    vs[var,:,:] = vp[var,:,:] / np.sqrt(3)
    rho[var,:,:] = 310 * np.power(vp[var,:,:],0.25)

vp.astype("float32",order="C").tofile(sys.argv[5])
vs.astype("float32",order="C").tofile(sys.argv[6])
rho.astype("float32",order="C").tofile(sys.argv[7])