import sys
import numpy as np
import matplotlib.pyplot as plt

def readBinaryVolume(dim1,dim2,dim3,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2*dim3)
        volume = np.reshape(data, [dim1,dim2,dim3], order='C')
    
    return volume

def readBinaryArray(dim,filename):
    with open(filename, 'rb') as f:    
        data = np.fromfile(filename, dtype= np.int32, count= dim)
    
    return data

nrecx = int(sys.argv[1])
nrecy = int(sys.argv[2])

nt = int(sys.argv[3])
dt = float(sys.argv[6])

xplain = int(sys.argv[4])
yplain = int(sys.argv[5])

time = np.arange(nt) * dt

nrecs = nrecx*nrecy
nshot = nrecx*nrecy

seismVx = readBinaryVolume(nt,nrecy,nrecx,"results/seismVx.bin")
seismVy = readBinaryVolume(nt,nrecy,nrecx,"results/seismVy.bin")
seismVz = readBinaryVolume(nt,nrecy,nrecx,"results/seismVz.bin")
seismPs = readBinaryVolume(nt,nrecy,nrecx,"results/seismPs.bin")

xrec = readBinaryArray(nrecs,"parameters/xrecPositions.bin")
yrec = readBinaryArray(nrecs,"parameters/yrecPositions.bin")

xsrc = readBinaryArray(nshot,"parameters/xsrcPositions.bin")
ysrc = readBinaryArray(nshot,"parameters/ysrcPositions.bin")

recx = np.reshape(xrec,[nrecy,nrecx])
recy = np.reshape(yrec,[nrecy,nrecx])

data = np.array([seismVx,seismVy,seismVz,seismPs])
titles = ["Vx","Vy","Vz","Pressure"]

plt.figure(1,figsize=(15,10))
for i in range(len(data)):
    plt.subplot(int("22"+str(i+1)))
    plt.imshow(data[i,:,yplain,:],aspect="auto",cmap="Greys")
    plt.title(titles[i]+" - Y slice")
plt.savefig("results/seismogramsSliceY.png",dpi=200,bbox_inches="tight")

plt.figure(2,figsize=(15,10))
for i in range(len(data)):
    plt.subplot(int("22"+str(i+1)))
    plt.imshow(data[i,:,:,xplain],aspect="auto",cmap="Greys")
    plt.title(titles[i]+" - X slice")
plt.savefig("results/seismogramsSliceX.png",dpi=200,bbox_inches="tight")

plt.figure(3,figsize=(15,10))
plt.scatter(recx[yplain,:],recy[yplain,:])
plt.scatter(recx[:,xplain],recy[:,xplain])
plt.scatter(recx[yplain,xplain],recy[yplain,xplain])
plt.scatter(xsrc[int(nshot/2)],ysrc[int(nshot/2)])
plt.title("Geometria de aquisição",fontsize=20)
plt.xlabel("Direção de nrecx",fontsize=15)
plt.ylabel("Direção de nrecy",fontsize=15)
plt.legend(["Receptores","Receptores","Traço coletado","Ponto de tiro"],fontsize=15)
plt.xlim([50,250])
plt.ylim([50,250])
plt.savefig("results/geometria.png",dpi=200,bbox_inches="tight")

plt.figure(4,figsize=(15,10))

plt.subplot(131)
plt.plot(seismPs[:,yplain,xplain],time)
plt.gca().invert_yaxis()
plt.xlabel("Amplitude",fontsize=15)
plt.ylabel("Tempo [s]",fontsize=15)

plt.subplot(132)
plt.plot(seismPs[500:,yplain,xplain],time[500:])
plt.gca().invert_yaxis()
plt.xlabel("Amplitude",fontsize=15)
plt.ylabel("Tempo [s]",fontsize=15)

plt.subplot(133)
plt.plot(seismPs[900:,yplain,xplain],time[900:])
plt.gca().invert_yaxis()
plt.xlabel("Amplitude",fontsize=15)
plt.ylabel("Tempo [s]",fontsize=15)

plt.savefig("results/nearOffsetTrace.png",dpi=200,bbox_inches="tight")

plt.show()