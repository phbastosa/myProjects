import sys
import numpy as np

model = np.loadtxt(sys.argv[1],dtype=int)

nz = int(sys.argv[2])
nodes = model[:,0]
veloc = model[:,1]

nabc = int(sys.argv[3])
parb = float(sys.argv[4])
nodes[1:] += nabc

nzz = nz + 2*nabc
vp = np.zeros(nzz)

vp[nodes[-1]:] = veloc[-1]
for i in range(1,len(nodes)):
    vp[nodes[i-1]:nodes[i]] = veloc[i-1]

factor = np.zeros(nabc)
for k in range(nabc):   
    factor[k] = np.exp(-pow(parb*(nabc-k),2.0))

damp = np.ones(nzz)
damp[:nabc] = factor
damp[nzz-nabc:] = factor[::-1]

vp.astype("float32",order="C").tofile(sys.argv[5])
damp.astype("float32",order="C").tofile(sys.argv[6])
