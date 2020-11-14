import sys
import numpy as np
import matplotlib.pyplot as plt

def readbinaryfile(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='C')
    return matrix.T

nx = 300
ny = 300
nz = 300
borda = 50

# filename = 

nxx = nx + 2*borda
nyy = ny + 2*borda
nzz = nz + 2*borda

# model = readbinaryfile(nx,nz,filename)
model = np.random.rand(nz,nx,ny)

modelExpanded = np.zeros((nzz,nxx,nyy))

modelExpanded[borda:nzz-borda,borda:nxx-borda,borda:nyy-borda] = model[:,:,:]

for i in range(borda):
    modelExpanded[i,borda:nxx-borda,borda:nyy-borda] = model[0,:,:]
    modelExpanded[nzz-i-1,borda:nxx-borda,borda:nyy-borda] = model[-1,:,:]

for i in range(borda):
    modelExpanded[:,i,borda:nyy-borda] = modelExpanded[:,borda,borda:nyy-borda]
    modelExpanded[:,nxx-i-1,borda:nyy-borda] = modelExpanded[:,nxx-borda-1,borda:nyy-borda]

for i in range(borda):
    modelExpanded[:,:,i] = modelExpanded[:,:,borda]
    modelExpanded[:,:,nyy-i-1] = modelExpanded[:,:,nyy-borda-1]

# modelName = 

plt.show()
# modelExpanded.astype("float32",order="C").tofile(modelName)