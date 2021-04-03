import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle

def readbinaryfile(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='F')
    return matrix.T

nx = 500
nt = 2000
dt = 0.001
dx = 5.0

timeCut = 750
trace = 250

seismRicker = readbinaryfile(nx,nt,"seismPressure_ricker.bin")
seismIntRicker = readbinaryfile(nx,nt,"seismPressure_intRicker.bin")
seismHankelIntRicker = readbinaryfile(nx,nt,"seismPressure_hankelIntRicker.bin")

seismRicker = seismRicker[:timeCut,:]
seismIntRicker = seismIntRicker[:timeCut,:]
seismHankelIntRicker = seismHankelIntRicker[:timeCut,:]

traceRicker = seismRicker[:,trace]
traceIntRicker = seismIntRicker[:,trace]
traceHankelIntRicker = seismHankelIntRicker[:,trace]
traceHankelIntRicker[:164] = 0
traceHankelIntRicker[480:] = 0

t = np.arange(timeCut) * dt

tloc = np.linspace(0,timeCut,5)
tlab = np.around(tloc * dt,decimals=2)

xloc = np.linspace(0,nx,5,dtype=int)
xlab = np.linspace(0,nx*dx,5,dtype=int)
plt.figure(2,figsize=(15, 8))
G = gridspec.GridSpec(3, 5)

axes_1 = plt.subplot(G[0,:2])
plt.imshow(seismRicker,aspect="auto",cmap="Greys")
rect1 = Rectangle((250,0),0.001,2000,linewidth=1,edgecolor='r',facecolor='none')
plt.gca().add_patch(rect1)
plt.yticks(tloc,tlab)
plt.xticks(xloc,xlab)
plt.title("Usando a Ricker")
plt.ylabel("Tempo [s]")

axes_2 = plt.subplot(G[1,:2])
plt.imshow(seismIntRicker,aspect="auto",cmap="Greys")
rect2 = Rectangle((250,0),0.001,2000,linewidth=1,edgecolor='r',facecolor='none')
plt.gca().add_patch(rect2)
plt.yticks(tloc,tlab)
plt.xticks(xloc,xlab)
plt.title("Usando a integral da Ricker")
plt.ylabel("Tempo [s]")

axes_3 = plt.subplot(G[2,:2])
plt.imshow(seismHankelIntRicker,aspect="auto",cmap="Greys")
rect3 = Rectangle((250,0),0.001,2000,linewidth=1,edgecolor='r',facecolor='none')
plt.gca().add_patch(rect3)
plt.yticks(tloc,tlab)
plt.xticks(xloc,xlab)
plt.title("Usando a integral da Ricker com atraso")
plt.ylabel("Tempo [s]")
plt.xlabel("Distância [m]")

axes_4 = plt.subplot(G[0:,2])
plt.plot(traceRicker,t)
plt.gca().invert_yaxis()
plt.title("Ricker")
plt.xlabel("Amplitude")
plt.ylabel("Tempo [s]")

axes_5 = plt.subplot(G[0:,3])
plt.plot(traceIntRicker,t)
plt.gca().invert_yaxis()
plt.title("Integral da Ricker")
plt.xlabel("Amplitude")
plt.ylabel("Tempo [s]")

axes_6 = plt.subplot(G[0:,4])
plt.plot(traceHankelIntRicker,t)
plt.gca().invert_yaxis()
plt.title("Resultado")
plt.xlabel("Amplitude")
plt.ylabel("Tempo [s]")

plt.tight_layout()
plt.savefig("analiseFontes.png",dpi=200,bbox_inches="tight")
plt.show()