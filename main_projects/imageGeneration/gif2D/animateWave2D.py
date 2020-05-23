import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def readbinaryfile(dim1,dim2,dim3,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2*dim3)
        matrix = np.reshape(data, [dim1,dim2,dim3], order='F')
    return matrix

# Lendo os snapshots brutos
nx = 583
nz = 341

dx = 30
dz = 20

borda = 100
snap = 67

snapshots = readbinaryfile(nx,nz,snap,"snapshots.bin")

# Coletando cada snapshot e o transpondo
SemBorda = snapshots[borda:nx-borda,borda:nz-borda,:]

campoLimpo = np.zeros((nz - 2*borda, nx - 2*borda, snap))
for i in range(snap):
    campoLimpo[:,:,i] = SemBorda[:,:,i].T

fig = plt.figure()

plt.title("Propagação de ondas no modelo Marmousi acústico")
plt.xlabel("Distância [m]")
plt.ylabel("Profundidade [m]")

qLabelsX = 13
qLabelsZ = 17

locks_x = np.linspace(0,nx,qLabelsX)
labels_x = []
for j in locks_x:
    labels_x.append(str(int(j*dx)))

locks_z = np.linspace(0,nz,qLabelsZ)
labels_z = []
for i in locks_z:
    labels_z.append(str(int(i*dz)))

plt.xticks(locks_x,labels_x)
plt.yticks(locks_z,labels_z)

ims = []
for i in range(snap):
    im = plt.imshow(campoLimpo[:,:,i],cmap="Greys",animated=True)
    ims.append([im])

ani = animation.ArtistAnimation(fig, ims, interval=40, blit=True, repeat=True)

writer = animation.FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)

ani.save("wave2D.mp4", writer=writer)

plt.show()