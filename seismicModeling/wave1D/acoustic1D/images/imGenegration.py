import numpy as np
import matplotlib.pyplot as plt
import imageio

def readBinaryMatrix(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='F')
    return matrix.T

def readBinaryVector(dim,filename):
    with open(filename, 'rb') as f:    
        data = np.fromfile(filename, dtype=np.float32, count= dim)
    return data

nz = 400
dz = 10.0

nt = 7001
dt = 0.0005

nSnap = 1000

snap = readBinaryMatrix(nz,nSnap,"../results/snapshots.bin")
seis = readBinaryVector(nt,"../results/seismogram.bin")
vp   = readBinaryVector(nz,"../model/velocities.bin")

depth = np.arange(nz) * dz
times = np.arange(nt) * dt

plt.figure(1,figsize=(15,7))
plt.subplot(131)
plt.plot(vp,depth)
plt.gca().invert_yaxis()
plt.title("Velocity model",fontsize=18)
plt.xlabel("Velocities [m/s]",fontsize=12)
plt.ylabel("Depth [m]",fontsize=12)

i = 480

plt.subplot(132)
plt.plot(snap[i,:],depth)
plt.xlim([-100,100])
plt.gca().set_yticklabels([])
plt.gca().invert_yaxis()
plt.title(f"Wave field at {i*dt:.2f} sec",fontsize=18)
plt.xlabel("Amplitude",fontsize=12)

plt.subplot(133)
plt.plot(seis,times)
plt.xlim([-100,100])
plt.gca().invert_yaxis()
plt.title("Surface seismogram",fontsize=18)
plt.xlabel("Amplitude",fontsize=12)
plt.ylabel("Times [s]",fontsize=12)
plt.savefig(f"imageGithub{i}.png",bbox_inches="tight")

plt.show()