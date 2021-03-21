import sys
import numpy as np

def readBinaryFloatMatrix(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='C')
    return matrix
    
nx = int(sys.argv[1])
nz = int(sys.argv[2])
nabc = int(sys.argv[3])

vpPath = sys.argv[4]
vsPath = sys.argv[5]
rhoPath = sys.argv[6]

nxx = nx + 2*nabc
nzz = nz + 2*nabc

vp = readBinaryFloatMatrix(nz,nx,vpPath)
vs = readBinaryFloatMatrix(nz,nx,vsPath)
rho = readBinaryFloatMatrix(nz,nx,rhoPath)

vpExpanded = np.zeros((nzz,nxx))
vsExpanded = np.zeros((nzz,nxx))
rhoExpanded = np.zeros((nzz,nxx))

vpExpanded[nabc:nzz-nabc,nabc:nxx-nabc] = vp[:,:]
vsExpanded[nabc:nzz-nabc,nabc:nxx-nabc] = vs[:,:]
rhoExpanded[nabc:nzz-nabc,nabc:nxx-nabc] = rho[:,:]

for i in range(nabc):
    vpExpanded[i,nabc:nxx-nabc] = vp[0,:]
    vpExpanded[nzz-i-1,nabc:nxx-nabc] = vp[-1,:]

    vsExpanded[i,nabc:nxx-nabc] = vs[0,:]
    vsExpanded[nzz-i-1,nabc:nxx-nabc] = vs[-1,:]

    rhoExpanded[i,nabc:nxx-nabc] = rho[0,:]
    rhoExpanded[nzz-i-1,nabc:nxx-nabc] = rho[-1,:]

for i in range(nabc):
    vpExpanded[:,i] = vpExpanded[:,nabc]
    vpExpanded[:,nxx-i-1] = vpExpanded[:,nxx-nabc-1]

    vsExpanded[:,i] = vsExpanded[:,nabc]
    vsExpanded[:,nxx-i-1] = vsExpanded[:,nxx-nabc-1]

    rhoExpanded[:,i] = rhoExpanded[:,nabc]
    rhoExpanded[:,nxx-i-1] = rhoExpanded[:,nxx-nabc-1]

vpInput = vpExpanded[:,:]
vsInput = vsExpanded[:,:]
rhoInput = rhoExpanded[:,:]

vpInput.astype("float32",order="C").tofile(sys.argv[7])
vsInput.astype("float32",order="C").tofile(sys.argv[8])
rhoInput.astype("float32",order="C").tofile(sys.argv[9])
