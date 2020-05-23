from skimage import exposure
import numpy as np
import matplotlib.pyplot as plt

def readbinaryfile(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='C')
    return matrix
dt = 0.00012
spread = 72
nt = 5000
t = np.array([0.0, 0.0, 0.12, 0.24, 0.36, 0.48, 0.60])

raw_shot = readbinaryfile(nt,spread,"Raw_shots_50Hz/Tiro_100.bin")
hom_shot = readbinaryfile(nt,spread,"Homo_shots_50Hz/Shot_100.bin")
mig_shot = readbinaryfile(nt,spread,"Mig_shots_50Hz/VR_mig_100.bin")

perc1 = np.percentile(raw_shot,[.5, 99.6])
perc2 = np.percentile(hom_shot,[.5, 99.6])
perc3 = np.percentile(mig_shot,[.5, 99.6])

raw_shot = exposure.rescale_intensity(raw_shot,in_range=(perc1[0],perc1[1]),out_range=(0,255))
hom_shot = exposure.rescale_intensity(hom_shot,in_range=(perc1[0],perc1[1]),out_range=(0,255))
mig_shot = exposure.rescale_intensity(mig_shot,in_range=(perc1[0],perc1[1]),out_range=(0,255))

cont = 0
images = np.array([raw_shot,hom_shot,mig_shot])

fig, axes = plt.subplots(ncols=3,nrows=1,dpi=150)
fig.set_size_inches(10,8)

axes[0].imshow(images[0],aspect="auto",cmap="Greys")
axes[0].set_ylabel("Time [s]",fontsize=15)
axes[0].set_xlabel("Traces",fontsize=15)
axes[0].set_yticklabels(t)    

axes[1].imshow(images[1],aspect="auto",cmap="Greys")
axes[1].set_yticklabels([])
axes[1].set_xlabel("Traces",fontsize=15)

axes[2].imshow(images[2],aspect="auto",cmap="Greys")
axes[2].set_yticklabels([])
axes[2].set_xlabel("Traces",fontsize=15)

for ax in axes.flat:
    im = ax.imshow(images[cont],aspect="auto",cmap="Greys")
    cont += 1
cbar = fig.colorbar(im,ax=axes.flat,fraction=0.025)
cbar.set_label("Amplitude",fontsize = 15, rotation = 90)
plt.savefig("elastic_data.png")
plt.show()

 
