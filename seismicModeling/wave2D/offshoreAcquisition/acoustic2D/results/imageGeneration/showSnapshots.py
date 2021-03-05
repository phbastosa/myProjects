import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def readbinaryfile(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='F')
    return matrix

nx = int(sys.argv[1])
nz = int(sys.argv[2])
dh = float(sys.argv[3])
nsnaps = int(sys.argv[4])

snapshots = readbinaryfile(nx,nz*nsnaps,sys.argv[5])

snaps = np.zeros((nz,nx,nsnaps))

for i in range(nsnaps):
    snaps[:,:,i] = snapshots[:,i*nz:i*nz + nz].T 

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (15,5), dpi = 96)

xloc = np.linspace(0,nx,11,dtype=int)
xlab = np.linspace(0,nx*dh,11,dtype=int)

zloc = np.linspace(0,nz,5)
zlab = np.linspace(0,nz*dh,5,dtype=int)

ims = []
for i in range(4,nsnaps):
   
    im = ax.imshow(snaps[:,:,i], animated=True,aspect="auto",cmap="Greys")

    if i == 4:
        ax.imshow(snaps[:,:,i],aspect="auto",cmap="Greys")  # show an initial one first
    
    ims.append([im])

ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True, repeat_delay=1000)

ax.set_title("Snapshots",fontsize=20)
ax.set_ylabel("Profundidade [m]",fontsize=15)
ax.set_xlabel("Distância [m]",fontsize=15)

plt.xticks(xloc,xlab)
plt.yticks(zloc,zlab)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Paulo'), bitrate=1000)

ani.save('snapshots.mp4', writer=writer)
