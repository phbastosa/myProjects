import sys
import numpy as np

def readbinaryfile(dim1,dim2,dim3,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2*dim3)
        volume = np.reshape(data, [dim1,dim2,dim3], order='C')
    return volume

nx = int(sys.argv[1])
ny = int(sys.argv[2])
nz = int(sys.argv[3])
borda = int(sys.argv[4])

vpFile = sys.argv[5]
vsFile = sys.argv[6]
rhoFile = sys.argv[7]

nxx = nx + 2*borda
nyy = ny + 2*borda
nzz = nz + 2*borda

vpModel = readbinaryfile(nz,nx,ny,vpFile)
vsModel = readbinaryfile(nz,nx,ny,vsFile)
rhoModel = readbinaryfile(nz,nx,ny,rhoFile)

vpModelExpanded = np.zeros((nzz,nxx,nyy))
vsModelExpanded = np.zeros((nzz,nxx,nyy))
rhoModelExpanded = np.zeros((nzz,nxx,nyy))

vpModelExpanded[borda:nzz-borda,borda:nxx-borda,borda:nyy-borda] = vpModel[:,:,:]
vsModelExpanded[borda:nzz-borda,borda:nxx-borda,borda:nyy-borda] = vsModel[:,:,:]
rhoModelExpanded[borda:nzz-borda,borda:nxx-borda,borda:nyy-borda] = rhoModel[:,:,:]

for i in range(borda):
    vpModelExpanded[i,borda:nxx-borda,borda:nyy-borda] = vpModel[0,:,:]
    vpModelExpanded[nzz-i-1,borda:nxx-borda,borda:nyy-borda] = vpModel[-1,:,:]

    vsModelExpanded[i,borda:nxx-borda,borda:nyy-borda] = vsModel[0,:,:]
    vsModelExpanded[nzz-i-1,borda:nxx-borda,borda:nyy-borda] = vsModel[-1,:,:]

    rhoModelExpanded[i,borda:nxx-borda,borda:nyy-borda] = rhoModel[0,:,:]
    rhoModelExpanded[nzz-i-1,borda:nxx-borda,borda:nyy-borda] = rhoModel[-1,:,:]

for i in range(borda):
    vpModelExpanded[:,i,borda:nyy-borda] = vpModelExpanded[:,borda,borda:nyy-borda]
    vpModelExpanded[:,nxx-i-1,borda:nyy-borda] = vpModelExpanded[:,nxx-borda-1,borda:nyy-borda]

    vsModelExpanded[:,i,borda:nyy-borda] = vsModelExpanded[:,borda,borda:nyy-borda]
    vsModelExpanded[:,nxx-i-1,borda:nyy-borda] = vsModelExpanded[:,nxx-borda-1,borda:nyy-borda]

    rhoModelExpanded[:,i,borda:nyy-borda] = rhoModelExpanded[:,borda,borda:nyy-borda]
    rhoModelExpanded[:,nxx-i-1,borda:nyy-borda] = rhoModelExpanded[:,nxx-borda-1,borda:nyy-borda]

for i in range(borda):
    vpModelExpanded[:,:,i] = vpModelExpanded[:,:,borda]
    vpModelExpanded[:,:,nyy-i-1] = vpModelExpanded[:,:,nyy-borda-1]

    vsModelExpanded[:,:,i] = vsModelExpanded[:,:,borda]
    vsModelExpanded[:,:,nyy-i-1] = vsModelExpanded[:,:,nyy-borda-1]

    rhoModelExpanded[:,:,i] = rhoModelExpanded[:,:,borda]
    rhoModelExpanded[:,:,nyy-i-1] = rhoModelExpanded[:,:,nyy-borda-1]

vpModelExpanded.astype("float32",order="C").tofile(sys.argv[8])
vsModelExpanded.astype("float32",order="C").tofile(sys.argv[9])
rhoModelExpanded.astype("float32",order="C").tofile(sys.argv[10])