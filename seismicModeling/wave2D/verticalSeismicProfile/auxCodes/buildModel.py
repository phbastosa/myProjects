import sys
import numpy as np
import matplotlib.pyplot as plt

depth = sys.argv[1]
vpMod = sys.argv[2]

depth = np.array(depth[1:-1].split(","),dtype=int) 
vpMod = np.array(vpMod[1:-1].split(","),dtype=int)

vsMod = vpMod / 1.7
rhoMod = 310 * vpMod ** 0.25

vsMod[0] = 0.0
rhoMod[0] = 1000.0

nx = int(sys.argv[3])
nz = int(sys.argv[4])

vp = np.zeros((nz,nx))
vs = np.zeros((nz,nx))
rho = np.zeros((nz,nx))

vp[:int(depth[0]),:]  = np.ones((int(depth[0]),nx)) * vpMod[0]
vs[:int(depth[0]),:]  = np.ones((int(depth[0]),nx)) * vsMod[0]
rho[:int(depth[0]),:] = np.ones((int(depth[0]),nx)) * rhoMod[0]

for i in range(1,len(depth)):
    var = slice(int(depth[i-1]),int(depth[i]))
    vp[var,:]  = np.ones((int(depth[i]-depth[i-1]),nx)) * vpMod[i]
    vs[var,:]  = np.ones((int(depth[i]-depth[i-1]),nx)) * vsMod[i]
    rho[var,:] = np.ones((int(depth[i]-depth[i-1]),nx)) * rhoMod[i]

vp.astype("float32",order="C").tofile(sys.argv[5])
vs.astype("float32",order="C").tofile(sys.argv[6])
rho.astype("float32",order="C").tofile(sys.argv[7])