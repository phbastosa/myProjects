import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def readbinaryfile(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='F')
    return matrix

nx = 1700
nz = 351
dh = 10
nsnaps = 100

snapshots = readbinaryfile(nx,nz*nsnaps,"../snapshots/snaps_shot_1.bin")

snaps = np.zeros((nz,nx,nsnaps))

for i in range(nsnaps):
    snaps[:,:,i] = snapshots[:,i*nz:i*nz + nz].T 

fig, ax = plt.subplots()

ims = []
for i in range(4,nsnaps):
    
    im = ax.imshow(snaps[:,:,i], animated=True,aspect="auto",cmap="Greys")

    if i == 4:
        ax.imshow(snaps[:,:,i],aspect="auto",cmap="Greys")  # show an initial one first
    
    ims.append([im])

ani = animation.ArtistAnimation(fig, ims, interval=80, blit=True, repeat_delay=1000)

# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Paulo'), bitrate=1000)

# ani.save('../../snapshots.mp4', writer=writer)

plt.show()
