import numpy as np
import matplotlib.pyplot as plt

def readbinaryfile(dim1,dim2,dim3,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2*dim3)
        volume = np.reshape(data, [dim1,dim2,dim3], order='C')
    return volume

nrecx = 31
nrecy = 51
nt = 2000

seismVx = readbinaryfile(nt,nrecy,nrecx,"seismVx.bin")
seismVy = readbinaryfile(nt,nrecy,nrecx,"seismVy.bin")
seismVz = readbinaryfile(nt,nrecy,nrecx,"seismVz.bin")
seismPs = readbinaryfile(nt,nrecy,nrecx,"seismPs.bin")

seismP = readbinaryfile(nt,nrecy,nrecx,"seismP.bin")
seismSv = readbinaryfile(nt,nrecy,nrecx,"seismSv.bin")
seismShx = readbinaryfile(nt,nrecy,nrecx,"seismShx.bin")
seismShy = readbinaryfile(nt,nrecy,nrecx,"seismShy.bin")

time = np.arange(nt)

plt.subplot(241)
plt.imshow(seismVx[:,0,:],aspect="auto",cmap="Greys")
# plt.plot(seismVx[:,10,0],time)
# plt.gca().invert_yaxis()
plt.title("Vx")

plt.subplot(242)
plt.imshow(seismVy[:,0,:],aspect="auto",cmap="Greys")
# plt.plot(seismVy[:,10,0],time)
# plt.gca().invert_yaxis()
plt.title("Vy")

plt.subplot(243)
plt.imshow(seismVz[:,0,:],aspect="auto",cmap="Greys")
# plt.plot(seismVz[:,10,0],time)
# plt.gca().invert_yaxis()
plt.title("Vz")

plt.subplot(244)
plt.imshow(seismPs[:,0,:],aspect="auto",cmap="Greys")
# plt.plot(seismPs[:,10,0],time)
# plt.gca().invert_yaxis()
plt.title("Pressure")

plt.subplot(245)
plt.imshow(seismShx[:,0,:],aspect="auto",cmap="Greys")
# plt.plot(seismVx[:,10,0],time)
# plt.gca().invert_yaxis()
plt.title("Shx")

plt.subplot(246)
plt.imshow(seismShy[:,0,:],aspect="auto",cmap="Greys")
# plt.plot(seismVy[:,10,0],time)
# plt.gca().invert_yaxis()
plt.title("Shy")

plt.subplot(247)
plt.imshow(seismSv[:,0,:],aspect="auto",cmap="Greys")
# plt.plot(seismVz[:,10,0],time)
# plt.gca().invert_yaxis()
plt.title("Sv")

plt.subplot(248)
plt.imshow(seismP[:,0,:],aspect="auto",cmap="Greys")
# plt.plot(seismPs[:,10,0],time)
# plt.gca().invert_yaxis()
plt.title("P")

plt.show()