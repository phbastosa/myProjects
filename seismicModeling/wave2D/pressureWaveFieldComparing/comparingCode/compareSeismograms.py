import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from matplotlib.patches import Rectangle

def readbinaryfile(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='F')
    return matrix.T

nx = int(sys.argv[1])
dx = float(sys.argv[2])

nt = int(sys.argv[3])
dt = float(sys.argv[4])

trace = int(nx/2)
timeCut = int(sys.argv[5])
lagSource = int(int(sys.argv[6])/2-1)

seismSA = readbinaryfile(nx,nt,"results/seismograms/pressureScalarAcoustic2D.bin")
seismVA = readbinaryfile(nx,nt,"results/seismograms/pressureVectorialAcoustic2D.bin")
seismEI = readbinaryfile(nx,nt,"results/seismograms/pressureElasticIsotropic2D.bin")

seismSA = seismSA[lagSource:timeCut,:]
seismVA = seismVA[lagSource:timeCut,:]
seismEI = seismEI[lagSource:timeCut,:]

traceSA = seismSA[:,trace]
traceVA = seismVA[:,trace]
traceEI = seismEI[:,trace]

specSA = np.fft.fft(traceSA)
specVA = np.fft.fft(traceVA)
specEI = np.fft.fft(traceEI)

freq = np.fft.fftfreq(timeCut-lagSource,dt)
mask = freq >= 0.0

t = np.arange(timeCut-lagSource) * dt

tloc = np.linspace(0,timeCut-lagSource,5)
tlab = np.around(tloc * dt,decimals=2)

xloc = np.linspace(0,nx,5,dtype=int)
xlab = np.linspace(0,nx*dx,5,dtype=int)

plt.figure(1,figsize=(15, 8))

G = gridspec.GridSpec(3, 4)

axes_1 = plt.subplot(G[0,:2])
plt.imshow(seismSA,aspect="auto",cmap="Greys")
rect1 = Rectangle((250,0),0.001,2000,linewidth=1,edgecolor='r',facecolor='none')
plt.gca().add_patch(rect1)
plt.yticks(tloc,tlab)
plt.xticks(xloc,xlab)
plt.title("Scalar Acoustic")
plt.ylabel("Time [s]")

axes_2 = plt.subplot(G[1,:2])
plt.imshow(seismVA,aspect="auto",cmap="Greys")
rect2 = Rectangle((250,0),0.001,2000,linewidth=1,edgecolor='r',facecolor='none')
plt.gca().add_patch(rect2)
plt.yticks(tloc,tlab)
plt.xticks(xloc,xlab)
plt.title("Vectorial Acoustic")
plt.ylabel("Time [s]")

axes_3 = plt.subplot(G[2,:2])
plt.imshow(seismEI,aspect="auto",cmap="Greys")
rect3 = Rectangle((250,0),0.001,2000,linewidth=1,edgecolor='r',facecolor='none')
plt.gca().add_patch(rect3)
plt.yticks(tloc,tlab)
plt.xticks(xloc,xlab)
plt.title("Elastic Isotropic")
plt.ylabel("Time [s]")
plt.xlabel("Displacement [m]")

axes_4 = plt.subplot(G[0:,2])
plt.plot(traceSA,t)
plt.plot(traceVA,t)
plt.plot(traceEI,t)
plt.gca().invert_yaxis()
plt.title("Traces")
plt.xlabel("Amplitude")
plt.ylabel("Time [s]")
plt.legend(["Scalar acoustic","Vectorial acoustic","Elastic isotropic"],loc="lower right")

axes_5 = plt.subplot(G[0:,3])
plt.plot(np.abs(specSA[mask]),freq[mask])
plt.plot(np.abs(specVA[mask]),freq[mask])
plt.plot(np.abs(specEI[mask]),freq[mask])
plt.ylim([0,30])
plt.gca().invert_yaxis()
plt.title("Spectrum")
plt.ylabel("Frequency [Hz]")
plt.xlabel("Amplitude")
plt.legend(["Scalar acoustic","Vectorial acoustic","Elastic isotropic"],loc="lower right")

plt.tight_layout()
plt.savefig("results/modelingAnalysis.png",dpi=200,bbox_inches="tight")
plt.show()