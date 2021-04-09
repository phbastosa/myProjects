import numpy as np
import matplotlib.pyplot as plt

def readBinaryMatrix(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype=np.float32, count=dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='C')
    return matrix

def maxCatch(vector):
    m = 0
    k = 0
    for i in range(len(vector)):
        if vector[i] > m:
            m = vector[i]
            k = i
    return k

nt = 10000
dt = 0.0005
spread = 101
nshots = 201

AS = readBinaryMatrix(nshots*nt,spread,"../acoustic_marmousi2_dh5_dataset.bin")
AV = readBinaryMatrix(nshots*nt,spread,"../acousticVec_marmousi2_dh5_dataset.bin")
EI = readBinaryMatrix(nshots*nt,spread,"../elasticIsotropic_marmousi2_dh5_dataset.bin")

seismicAS = np.zeros((nt,spread,nshots))
seismicAV = np.zeros((nt,spread,nshots))
seismicEI = np.zeros((nt,spread,nshots))

for i in range(nshots):
    seismicAS[:,:,i] = AS[i*nt:i*nt+nt,:]
    seismicAV[:,:,i] = AV[i*nt:i*nt+nt,:]
    seismicEI[:,:,i] = EI[i*nt:i*nt+nt,:]

lag = 300
shot = 100
trace = 50

time = np.arange(nt) * dt

traceAS = seismicAS[:,trace,shot]
traceAV = seismicAV[:,trace,shot]
traceEI = seismicEI[:,trace,shot]

tAS = (maxCatch(traceAS) - lag) * dt
tAV = (maxCatch(traceAV) - lag) * dt
tEI = (maxCatch(traceEI) - lag) * dt 

plt.figure(1,figsize=(4,10))

plt.plot(traceAS,time)
plt.plot(traceAV,time)
plt.plot(traceEI,time)

plt.gca().invert_yaxis() 
plt.show()