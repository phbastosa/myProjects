import sys
import numpy as np

def readbinaryfile(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='C')
    return matrix.T

nx = int(sys.argv[1])
nz = int(sys.argv[2])
borda = int(sys.argv[3])

filename = sys.argv[4]

nxx = nx + 2*borda
nzz = nz + 2*borda

model = readbinaryfile(nx,nz,filename)

modelExpanded = np.zeros((nzz,nxx))

modelExpanded[borda:nzz-borda,borda:nxx-borda] = model[:,:]

for i in range(borda):
    modelExpanded[i,borda:nxx-borda] = model[0,:]
    modelExpanded[nzz-i-1,borda:nxx-borda] = model[-1,:]

for i in range(borda):
    modelExpanded[:,i] = modelExpanded[:,borda]
    modelExpanded[:,nxx-i-1] = modelExpanded[:,nxx-borda-1]

modelExpanded.astype("float32",order="C").tofile(sys.argv[5])