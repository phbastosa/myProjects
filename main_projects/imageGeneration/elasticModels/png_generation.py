import numpy as np
import matplotlib.pyplot as plt

def readbinaryfile(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='F')
    return matrix

nx = 2500
nz = 200

vp = readbinaryfile(nz,nx,"Modelo_vp_200x2500.bin")
vs = readbinaryfile(nz,nx,"Modelo_vs_200x2500.bin")
rho = readbinaryfile(nz,nx,"Modelo_rho_200x2500.bin")

cont = 0
images = np.array([vp,vs,rho])

fig, axes = plt.subplots(nrows = 3, ncols = 1, figsize = (10,5), dpi = 150)

axes[0].set_ylabel("Depth [m]")
axes[0].set_xticklabels([])

axes[1].set_ylabel("Depth [m]")
axes[1].set_xticklabels([])

axes[2].set_ylabel("Depth [m]")

for ax in axes.flat:
    im = ax.imshow(images[cont],aspect="auto",vmin=0,vmax=3000)
    cont += 1

plt.xlabel("Distance [m]",fontsize=15)
cbar = fig.colorbar(im,ax=axes.flat,fraction=0.025)
cbar.set_label("Velocities $[m/s]$ and Densities $[kg/m^3]$",fontsize = 15, rotation = 90)
plt.savefig("VpVsRho_models.png")
plt.show()
 