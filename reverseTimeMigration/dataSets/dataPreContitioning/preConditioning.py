import numpy as np
import matplotlib.pyplot as plt

def readBinaryMatrix(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype=np.float32, count=dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='C')
    return matrix

def mute(matrix,p):
    y = np.arange(p[0][1],p[1][1],dtype=int)
    x = np.array(((p[0][0] - p[1][0]) * y + p[1][0]*p[0][1] - p[0][0]*p[1][1]) / (p[0][1] - p[1][1]),dtype=int)
    for i in range(len(x)):
        matrix[y[i],x[i]:] = 0.0 

nt = 10000
spread = 320
nshots = 357

# seismic = readBinaryMatrix(nshots*nt,spread,"../acoustic_marmousi2_dh5_dataset.bin")
seismic = readBinaryMatrix(nshots*nt,spread,"../acousticVec_marmousi2_dh5_dataset.bin")
# seismic = readBinaryMatrix(nshots*nt,spread,"../acoustic_marmousi2_dh5_dataset.bin")

points = [[0,700],[270,nt]]
for i in range(nshots):
    mute(seismic[i*nt:i*nt + nt,::-1],points)
    seismic[i*nt:i*nt + 700,::-1] = 0.0
    if i % 50 == 0:
        print(f"Onda direta removida até o sismograma {i}...")

# seismic.astype("float32",order="C").tofile("../acousticDataSetInput.bin")
seismic.astype("float32",order="C").tofile("../acousticVecDataSetInput.bin")
# seismic.astype("float32",order="C").tofile("../acousticVecDataSetInput.bin")

plt.show()