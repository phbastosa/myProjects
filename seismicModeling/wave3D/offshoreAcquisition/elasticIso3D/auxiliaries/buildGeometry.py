import sys
import numpy as np

dx = float(sys.argv[1])

nshot = int(sys.argv[2])
ishot = int(float(sys.argv[3]) / dx)
dshot = int(float(sys.argv[4]) / dx)

spread = int(sys.argv[5])
irecep = int(float(sys.argv[6]) / dx)
drecep = int(float(sys.argv[7]) / dx)

xsrc = np.zeros(nshot,dtype=int)
xrec = np.zeros((nshot,spread),dtype=int)

for i in range(nshot):
    xsrc[i] = ishot + dshot*i

for i in range(spread):        
    xrec[0,i] = irecep + drecep*i

for i in range(1,nshot):
    xrec[i,:] = xrec[i-1,:] + dshot

xsrc += int(sys.argv[8])
xrec += int(sys.argv[8])

zsrc = np.ones(nshot,dtype=int) * int(sys.argv[8]) + 1
zrec = np.ones(nshot*spread,dtype=int) * int(sys.argv[8]) 

xsrc.astype("int32",order="C").tofile(sys.argv[9])
zsrc.astype("int32",order="C").tofile(sys.argv[10])

xrec = np.reshape(xrec,[nshot*spread])
xrec.astype("int32",order="C").tofile(sys.argv[11])
zrec.astype("int32",order="C").tofile(sys.argv[12])
