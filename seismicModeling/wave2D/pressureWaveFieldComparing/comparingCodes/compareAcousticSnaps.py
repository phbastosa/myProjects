import numpy as np
import matplotlib.pyplot as plt

def readbinaryfile(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='F')
    return matrix

nx = 500
nSnaps = 11

nxx = nx*nSnaps

acoustic = readbinaryfile(nx,nxx,"pressure.bin")
txx = readbinaryfile(nx,nxx,"snapTxx.bin")
tzz = readbinaryfile(nx,nxx,"snapTzz.bin")
elastic = (txx + tzz) / 2

snapAcoustic = acoustic[:,3000:3500]
snapElastic  = elastic[:,3000:3500]

plt.figure(1)
plt.imshow(snapElastic - snapAcoustic)
plt.colorbar()

plt.figure(2)
plt.plot(snapAcoustic[:,250])
plt.plot(snapElastic[:,250])

acoustic = readbinaryfile(nx,nxx,"pressure_corr.bin")
txx = readbinaryfile(nx,nxx,"snapTxx_corr.bin")
tzz = readbinaryfile(nx,nxx,"snapTzz_corr.bin")
elastic = (txx + tzz) / 2

snapAcoustic = acoustic[:,3000:3500]
snapElastic  = elastic[:,3000:3500]

plt.figure(3)
plt.imshow(snapElastic - snapAcoustic)
plt.colorbar()

plt.figure(4)
plt.plot(snapAcoustic[:,250])
plt.plot(snapElastic[:,250])
plt.show()
