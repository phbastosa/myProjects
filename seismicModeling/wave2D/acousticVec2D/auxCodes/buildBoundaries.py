import sys
import numpy as np

def readBinaryFloatMatrix(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='C')
    return matrix.T
    
nx = int(sys.argv[1])
nz = int(sys.argv[2])
nabc = int(sys.argv[3])

vpPath = sys.argv[4]
rhoPath = sys.argv[5]

nxx = nx + 2*nabc
nzz = nz + 2*nabc

vp = readBinaryFloatMatrix(nx,nz,vpPath)
rho = readBinaryFloatMatrix(nx,nz,rhoPath)

vpExpanded = np.zeros((nzz,nxx))
rhoExpanded = np.zeros((nzz,nxx))

vpExpanded[nabc:nzz-nabc,nabc:nxx-nabc] = vp[:,:]
rhoExpanded[nabc:nzz-nabc,nabc:nxx-nabc] = rho[:,:]

for i in range(nabc):
    vpExpanded[i,nabc:nxx-nabc] = vp[0,:]
    vpExpanded[nzz-i-1,nabc:nxx-nabc] = vp[-1,:]

    rhoExpanded[i,nabc:nxx-nabc] = rho[0,:]
    rhoExpanded[nzz-i-1,nabc:nxx-nabc] = rho[-1,:]

for i in range(nabc):
    vpExpanded[:,i] = vpExpanded[:,nabc]
    vpExpanded[:,nxx-i-1] = vpExpanded[:,nxx-nabc-1]

    rhoExpanded[:,i] = rhoExpanded[:,nabc]
    rhoExpanded[:,nxx-i-1] = rhoExpanded[:,nxx-nabc-1]

vpExpanded.astype("float32",order="C").tofile(sys.argv[6])
rhoExpanded.astype("float32",order="C").tofile(sys.argv[7])
