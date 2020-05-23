from skimage import exposure
import numpy as np
import matplotlib.pyplot as plt

def readbinaryfile(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='F')
    return matrix

nx = 2500
nz = 200

image_30Hz = readbinaryfile(nz,nx,"results_30Hz/Image_without_Laplaciano_200x2500.bin")
image_45Hz = readbinaryfile(nz,nx,"results_45Hz/Image_without_Laplaciano_200x2500.bin")
image_60Hz = readbinaryfile(nz,nx,"results_60Hz/Image_without_Laplaciano_200x2500.bin")
image_75Hz = readbinaryfile(nz,nx,"results_75Hz/Image_without_Laplaciano_200x2500.bin")

# cont = 0
# images = np.array([image_30Hz,image_45Hz,image_60Hz,image_75Hz])

# fig, axes = plt.subplots(ncols = 2, nrows = 2, dpi = 150)
# fig.set_size_inches(20,6.5)

# axes[0][0].imshow(images[0],aspect="auto",cmap="Greys")
# axes[0][0].set_xticklabels([])
# axes[0][0].set_ylabel("Depht [m]",fontsize=15)
# axes[0][0].set_xlabel("a)",fontsize=15)

# axes[0][1].imshow(images[1],aspect="auto",cmap="Greys")
# axes[0][1].set_xticklabels([])
# axes[0][1].set_yticklabels([])
# axes[0][1].set_xlabel("b)",fontsize=15)

# axes[1][0].imshow(images[2],aspect="auto",cmap="Greys")
# axes[1][0].set_ylabel("Depht [m]",fontsize=15)
# axes[1][0].set_xlabel("Distance [m] \n c)",fontsize=15)

# axes[1][1].imshow(images[3],aspect="auto",cmap="Greys")
# axes[1][1].set_yticklabels([])
# axes[1][1].set_xlabel("Distance [m] \n d)",fontsize=15)

# for ax in axes.flat:
#     im = ax.imshow(images[cont],aspect="auto",cmap="Greys")
#     cont += 1

# cbar = fig.colorbar(im,ax=axes.flat,fraction=0.02)
# cbar.set_label("Amplitude square",fontsize = 15, rotation = 90)
# plt.savefig("Raw_images.png")
# plt.show()
 
image_30Hz = readbinaryfile(nz,nx,"results_30Hz/Compens_src_field_200x2500.bin")
image_45Hz = readbinaryfile(nz,nx,"results_45Hz/Compens_src_field_200x2500.bin")
image_60Hz = readbinaryfile(nz,nx,"results_60Hz/Compens_src_field_200x2500.bin")
image_75Hz = readbinaryfile(nz,nx,"results_75Hz/Compens_src_field_200x2500.bin")

perc1 = np.percentile(image_30Hz,[.5, 99.0])
perc2 = np.percentile(image_45Hz,[.5, 99.0])
perc3 = np.percentile(image_60Hz,[.5, 99.0])
perc4 = np.percentile(image_75Hz,[.5, 99.0])

image_30Hz = exposure.rescale_intensity(image_30Hz,in_range=(perc1[0],perc1[1]),out_range=(0,255))
image_45Hz = exposure.rescale_intensity(image_45Hz,in_range=(perc2[0],perc2[1]),out_range=(0,255))
image_60Hz = exposure.rescale_intensity(image_60Hz,in_range=(perc3[0],perc3[1]),out_range=(0,255))
image_75Hz = exposure.rescale_intensity(image_75Hz,in_range=(perc4[0],perc4[1]),out_range=(0,255))

plt.figure(figsize=(20,6.2),dpi=200)
plt.imshow(image_30Hz,aspect="auto",cmap="Greys")
plt.xlabel("Distance [m]",fontsize=30)
plt.ylabel("Depth [m]",fontsize=30)
cbar = plt.colorbar()
plt.savefig("image_norm_filt_30Hz.png")

plt.figure(figsize=(20,6.2),dpi=200)
plt.imshow(image_45Hz,aspect="auto",cmap="Greys")
plt.xlabel("Distance [m]",fontsize=30)
plt.ylabel("Depth [m]",fontsize=30)
cbar = plt.colorbar()
plt.savefig("image_norm_filt_45Hz.png")

plt.figure(figsize=(20,6.2),dpi=200)
plt.imshow(image_60Hz,aspect="auto",cmap="Greys")
plt.xlabel("Distance [m]",fontsize=30)
plt.ylabel("Depth [m]",fontsize=30)
cbar = plt.colorbar()
plt.savefig("image_norm_filt_60Hz.png")

plt.figure(figsize=(20,6.2),dpi=200)
plt.imshow(image_75Hz,aspect="auto",cmap="Greys")
plt.xlabel("Distance [m]",fontsize=30)
plt.ylabel("Depth [m]",fontsize=30)
cbar = plt.colorbar()
plt.savefig("image_norm_filt_75Hz.png")
