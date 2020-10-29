import numpy as np
import matplotlib.pyplot as plt

def readbinaryfile(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='C')
    return matrix

model = readbinaryfile(1180,3201,'sigsbee2a_stratigraphy_1180x3201.bin')

locks_x = np.linspace(0,3201,11)
labels_x = []
for j in locks_x:
    labels_x.append(str(int(j*10)))

locks_z = np.linspace(0,1180,11)
labels_z = []
for i in locks_z:
    labels_z.append(str(int(i*10)))

plt.figure(1,figsize=(10,5))
plt.imshow(model,aspect='auto',cmap='Greys')
plt.title("Modelo de velocidades Sigsbee")
plt.xlabel("Distância [amostras]")
plt.ylabel("Profundidade [amostras]")
cbar = plt.colorbar()
cbar.set_label("Velocidades [m/s]")

plt.figure(2,figsize=(10,5))
plt.imshow(model,aspect='auto',cmap='Greys')
plt.title("Modelo de velocidades Sigsbee")
plt.xlabel("Distância [m]")
plt.ylabel("Profundidade [m]")
cbar = plt.colorbar()
cbar.set_label("Velocidades [m/s]")
plt.xticks(locks_x,labels_x)
plt.yticks(locks_z,labels_z)
plt.show()