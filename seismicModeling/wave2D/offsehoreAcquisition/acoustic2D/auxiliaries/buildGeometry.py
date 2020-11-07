import sys
import numpy as np

dx = float(sys.argv[1])

nshot = int(sys.argv[2])
ishot = int(float(sys.argv[3]) / dx)
dshot = int(float(sys.argv[4]) / dx)

spread = int(sys.argv[5])
irecep = int(float(sys.argv[6]) / dx)
drecep = int(float(sys.argv[7]) / dx)

shots = np.zeros(nshot,dtype=int)
recps = np.zeros((nshot,spread),dtype=int)

for i in range(nshot):
    shots[i] = ishot + dshot*i

for i in range(spread):        
    recps[0,i] = irecep + drecep*i

for i in range(1,nshot):
    recps[i,:] = recps[i-1,:] + dshot

shots += int(sys.argv[10])
recps += int(sys.argv[10])

shots.astype("float32",order="C").tofile(sys.argv[8])

recps = np.reshape(recps,[nshot*spread])
recps.astype("float32",order="C").tofile(sys.argv[9])