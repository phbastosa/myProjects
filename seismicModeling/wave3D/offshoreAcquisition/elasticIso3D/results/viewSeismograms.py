import sys
import numpy as np
import matplotlib.pyplot as plt

def readbinaryfile(dim1,dim2,dim3,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2*dim3)
        volume = np.reshape(data, [dim1,dim2,dim3], order='C')
    return volume

nrecx = int(sys.argv[1])
nrecy = int(sys.argv[2])
nt = int(sys.argv[3])

seismVx = readbinaryfile(nt,nrecy,nrecx,"results/seismVx.bin")
seismVy = readbinaryfile(nt,nrecy,nrecx,"results/seismVy.bin")
seismVz = readbinaryfile(nt,nrecy,nrecx,"results/seismVz.bin")
seismPs = readbinaryfile(nt,nrecy,nrecx,"results/seismPs.bin")

seismP = readbinaryfile(nt,nrecy,nrecx,"results/seismP.bin")
seismSv = readbinaryfile(nt,nrecy,nrecx,"results/seismSv.bin")
seismShx = readbinaryfile(nt,nrecy,nrecx,"results/seismShx.bin")
seismShy = readbinaryfile(nt,nrecy,nrecx,"results/seismShy.bin")

time = np.arange(nt)

data = np.array([seismVx,seismVy,seismVz,seismPs,seismShx,seismShy,seismSv,seismP])
titles = ["Vx","Vy","Vz","Pressure","Shx","Shy","Sv","P"]

xplain = int(sys.argv[4])
yplain = int(sys.argv[5])

plt.figure(1,figsize=(15,10))
for i in range(len(data)):
    plt.subplot(int("24"+str(i+1)))
    plt.imshow(data[i,:,yplain,:],aspect="auto",cmap="Greys")
    plt.title(titles[i]+" X slice")

plt.figure(2,figsize=(15,10))
for i in range(len(data)):
    plt.subplot(int("24"+str(i+1)))
    plt.imshow(data[i,:,:,xplain],aspect="auto",cmap="Greys")
    plt.title(titles[i]+" Y slice")

plt.show()