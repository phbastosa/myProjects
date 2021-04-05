import sys
import numpy as np

nx = int(sys.argv[1])
nz = int(sys.argv[2])
abcLayer = int(sys.argv[3])
parameter = float(sys.argv[4])

nxx = nx + 2*abcLayer
nzz = nz + 2*abcLayer

damp = np.ones((nzz,nxx))
factor = np.zeros(abcLayer)

for k in range(abcLayer):   
    factor[k] = np.exp(-pow(parameter*(abcLayer-k),2))

for i in range(abcLayer,nzz-abcLayer):
    damp[i,:abcLayer] = factor
    damp[i,nxx-abcLayer:nxx] = factor[::-1]

for j in range(abcLayer,nxx-abcLayer):
    damp[:abcLayer,j] = factor
    damp[nzz-abcLayer:nzz,j] = factor[::-1]

for i in range(abcLayer):
    damp[i:abcLayer,i] = factor[i]
    damp[i,i:abcLayer] = factor[i]

    damp[i:abcLayer,nxx-i-1] = factor[i]
    damp[i,nxx-abcLayer:nxx-i-1] = factor[i]

    damp[nzz-abcLayer-1:nzz-i-1,i] = factor[i]
    damp[nzz-i-1,i:abcLayer] = factor[i]

    damp[nzz-abcLayer-1:nzz-i,nxx-i-1] = factor[i]
    damp[nzz-i-1,nxx-abcLayer-1:nxx-i] = factor[i]

damp.astype("float32",order="C").tofile(sys.argv[5])
