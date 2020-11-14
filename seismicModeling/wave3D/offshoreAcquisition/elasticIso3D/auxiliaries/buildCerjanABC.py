import sys
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits import mplot3d

nx = 300
ny = 300
nz = 300
abcLayer = 50
parameter = 0.005

nxx = nx + 2*abcLayer
nyy = ny + 2*abcLayer
nzz = nz + 2*abcLayer

dampXZ = np.ones((nzz,nxx))
dampYZ = np.ones((nzz,nyy))
dampXY = np.ones((nyy,nxx))

damp3D = np.ones((nzz,nxx,nyy))

factor = np.zeros(abcLayer)

for k in range(abcLayer):   
    factor[k] = np.exp(-pow(parameter*(abcLayer-k),2))

# Left - Right damping
for i in range(abcLayer,nzz-abcLayer):
    dampXZ[i,:abcLayer] = factor
    dampXZ[i,nxx-abcLayer:nxx] = factor[::-1]
    dampYZ[i,:abcLayer] = factor
    dampYZ[i,nyy-abcLayer:nyy] = factor[::-1]

for i in range(abcLayer,nyy-abcLayer):
    dampXY[i,:abcLayer] = factor
    dampXY[i,nxx-abcLayer:nxx] = factor[::-1]

# Up - Bottom damping
for j in range(abcLayer,nxx-abcLayer):
    dampXZ[:abcLayer,j] = factor
    dampXZ[nzz-abcLayer:nzz,j] = factor[::-1]
    
    dampXY[:abcLayer,j] = factor
    dampXY[nyy-abcLayer:nyy,j] = factor[::-1]

for k in range(abcLayer,nyy-abcLayer):    
    dampYZ[:abcLayer,k] = factor
    dampYZ[nzz-abcLayer:nzz,k] = factor[::-1]

for i in range(abcLayer):
    # Up left corner
    dampXZ[i:abcLayer,i] = factor[i]
    dampYZ[i:abcLayer,i] = factor[i]
    dampXY[i:abcLayer,i] = factor[i]
    dampXZ[i,i:abcLayer] = factor[i]    
    dampYZ[i,i:abcLayer] = factor[i]
    dampXY[i,i:abcLayer] = factor[i]
    
    # Up right corner
    dampXZ[i:abcLayer,nxx-i-1] = factor[i]
    dampYZ[i:abcLayer,nyy-i-1] = factor[i]
    dampXY[i:abcLayer,nxx-i-1] = factor[i]
    dampXZ[i,nxx-abcLayer:nxx-i-1] = factor[i]
    dampYZ[i,nyy-abcLayer:nyy-i-1] = factor[i]
    dampXY[i,nxx-abcLayer:nxx-i-1] = factor[i]

    # Bottom left
    dampXZ[nzz-i-1,i:abcLayer] = factor[i]    
    dampYZ[nzz-i-1,i:abcLayer] = factor[i]
    dampXY[nyy-i-1,i:abcLayer] = factor[i]
    dampXZ[nzz-abcLayer-1:nzz-i-1,i] = factor[i]
    dampYZ[nzz-abcLayer-1:nzz-i-1,i] = factor[i]
    dampXY[nyy-abcLayer-1:nyy-i-1,i] = factor[i]
    

    # Bottom right
    dampXZ[nzz-abcLayer-1:nzz-i,nxx-i-1] = factor[i]
    dampXZ[nzz-i-1,nxx-abcLayer-1:nxx-i] = factor[i]   
    dampYZ[nzz-abcLayer-1:nzz-i,nyy-i-1] = factor[i]
    dampYZ[nzz-i-1,nyy-abcLayer-1:nyy-i] = factor[i]
    dampXY[nyy-abcLayer-1:nyy-i,nxx-i-1] = factor[i]
    dampXY[nyy-i-1,nyy-abcLayer-1:nxx-i] = factor[i]

    # 3D corners
    damp3D[i:abcLayer,i:abcLayer,i] = factor[i]
    damp3D[i:abcLayer,i,i:abcLayer] = factor[i]
    damp3D[i,i:abcLayer,i:abcLayer] = factor[i]
    damp3D[i:abcLayer,i:abcLayer,nyy-i-1] = factor[i]
    damp3D[i:abcLayer,nxx-i-1,i:abcLayer] = factor[i]
    damp3D[nzz-i-1,i:abcLayer,i:abcLayer] = factor[i]
    damp3D[i:abcLayer,nxx-abcLayer-1:nxx-i,i] = factor[i]
    damp3D[nzz-abcLayer-1:nzz-i,i:abcLayer,i] = factor[i]
    damp3D[i:abcLayer,i,nyy-abcLayer-1:nyy-i] = factor[i]
    damp3D[nzz-abcLayer-1:nzz-i,i,i:abcLayer] = factor[i]
    damp3D[i,i:abcLayer,nyy-abcLayer-1:nyy-i] = factor[i]
    damp3D[i,nxx-abcLayer-1:nxx-i,i:abcLayer] = factor[i]
    damp3D[i:abcLayer,nxx-abcLayer-1:nxx-i,nyy-i-1] = factor[i]
    damp3D[nzz-abcLayer-1:nzz-i,i:abcLayer,nyy-i-1] = factor[i]    
    damp3D[i:abcLayer,nxx-i-1,nyy-abcLayer-1:nyy-i] = factor[i]
    damp3D[nzz-abcLayer-1:nzz-i,nxx-i-1,i:abcLayer] = factor[i]
    damp3D[nzz-i-1,i:abcLayer,nyy-abcLayer-1:nyy-i] = factor[i]
    damp3D[nzz-i-1,nxx-abcLayer-1:nxx-i,i:abcLayer] = factor[i]
    damp3D[nzz-abcLayer-1:nzz-i,nxx-abcLayer-1:nxx-i,i] = factor[i]
    damp3D[nzz-abcLayer-1:nzz-i,i,nyy-abcLayer-1:nyy-i] = factor[i]
    damp3D[i,nxx-abcLayer-1:nxx-i,nyy-abcLayer-1:nyy-i] = factor[i]
    damp3D[nzz-abcLayer-1:nzz-i,nxx-abcLayer-1:nxx-i,nyy-i-1] = factor[i]
    damp3D[nzz-abcLayer-1:nzz-i,nxx-i-1,nyy-abcLayer-1:nyy-i] = factor[i]
    damp3D[nzz-i-1,nxx-abcLayer-1:nxx-i,nyy-abcLayer-1:nyy-i] = factor[i]

# Prisms XZ
for i in range(abcLayer,nyy-abcLayer):
    damp3D[:,:,i] = dampXZ

# Prisms YZ
for i in range(abcLayer,nxx-abcLayer):
    damp3D[:,i,:] = dampYZ

# Prisms XY
for i in range(abcLayer,nzz-abcLayer):
    damp3D[i,:,:] = dampXY.T

plt.figure(1,figsize=(10,10))
plt.imshow(damp3D[310,:,:])
plt.colorbar()
plt.show()

damp3D.astype("float32",order="C").tofile("testeBorda3D.bin")