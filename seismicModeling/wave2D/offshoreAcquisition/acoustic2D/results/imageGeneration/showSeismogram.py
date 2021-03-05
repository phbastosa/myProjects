import sys
import numpy as np
import matplotlib.pyplot as plt

from skimage import exposure 

def perc(matrix,value):
    p = np.percentile(matrix,[.5, value])
    image = exposure.rescale_intensity(matrix,in_range=(p[0],p[1]),out_range=(0,255))    

    return image 

def readbinaryfile(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='F')
    return matrix.T

nx = 400
nt = 5000

dx = 10
dt = 0.001

seism = readbinaryfile(nx,nt,"../seismograms/seism_shot_1.bin")[:,::-1]

pseism = perc(seism,99.5)

plt.figure(1,figsize=(10,7))
cbar = plt.colorbar(plt.imshow(seism,aspect="auto",cmap="Greys"))
cbar.set_label("Amplitude",fontsize=15)
plt.imshow(pseism,aspect="auto",cmap="Greys")

tlocks = np.linspace(0,nt,11)
tlabel = np.around(tlocks * dt,decimals=2)

xlocks = np.linspace(0,nx,5)
xlabel = xlocks * dx

plt.xticks(xlocks,xlabel)
plt.yticks(tlocks,tlabel)

plt.title("Sismograma de pressão hidrostática",fontsize=20)
plt.xlabel("Distância [m]",fontsize=15)
plt.ylabel("Tempo [s]",fontsize=15)
plt.tight_layout()
plt.savefig("../../sismograma.png",dpi=200,bbox_inches="tight")
