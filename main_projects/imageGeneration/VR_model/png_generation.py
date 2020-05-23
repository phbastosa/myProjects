import numpy as np
import matplotlib.pyplot as plt

def readbinaryfile(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='F')
    return matrix

nx = 2500
nz = 200

vp = readbinaryfile(nz,nx,"modelo_vp_bruna1_topo_2500x200.bin")
vp_s = readbinaryfile(nz,nx,"modelo_vp_smooth_topo_2500x200.bin")

cont = 0
images = np.array([vp,vp_s])

fig, axes = plt.subplots(nrows = 2, ncols = 1, figsize = (10,5), dpi = 200)

axes[0].set_ylabel("Depth [m]")
axes[0].set_xticklabels([])

axes[1].set_ylabel("Depth [m]")

for ax in axes.flat:
    im = ax.imshow(images[cont],aspect="auto",vmin=0,vmax=3000)
    cont += 1

plt.xlabel("Distance [m]",fontsize=15)
cbar = fig.colorbar(im,ax=axes.flat,fraction=0.02)
cbar.set_label("Velocities [m/s]",fontsize=15,rotation=90) 
plt.savefig("smooth_model.png")
plt.show()
 