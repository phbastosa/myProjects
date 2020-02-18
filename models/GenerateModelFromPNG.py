import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

img = mpimg.imread('model.jpg')

<<<<<<< HEAD
maximo = 1490
minimo = 2800
=======
minimo = 1000
maximo = 4000
>>>>>>> ad1f0e4c82591805051ef40604a954b4fe930e1d

Lx = img.shape[0]
Ly = img.shape[1]

N = img[0:Lx,0:Ly,2].copy()

cor = np.arange(0,256)
vel = cor/256 * (maximo - minimo) + minimo
V = np.zeros((Lx,Ly))

for kk in range(len(cor)):
    i,j = np.where(cor[kk] == N)
    V[i,j] = vel[kk]

plt.figure(figsize=(10,5))
plt.subplot(121)
plt.imshow(img)
plt.colorbar()

plt.subplot(122)
plt.imshow(V)
plt.colorbar()
plt.show()

np.savetxt('model.txt',V)
