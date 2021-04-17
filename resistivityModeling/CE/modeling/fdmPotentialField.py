import numpy as np
import matplotlib.pyplot as plt

nx = 200
nz = 100

hE = 0.5
hW = 0.5
hN = 0.5
hS = 0.5

resistivity = np.ones((nz,nx)) * 10
condutivity = 1 / resistivity

qx = 1
qz = 1

q = np.zeros((nz,nx))

v_pas = np.zeros((nz,nx))
v_fut = np.zeros((nz,nx))

q[qz,qx] = 10

tol = 1e5

while True:    
    for i in range(1,nz-1):
        for j in range(1,nx-1):
        
            E = 2 * (condutivity[i,j] + condutivity[i,j+1]) * (hE*(hE + hW))**(-1)
            N = 2 * (condutivity[i,j] + condutivity[i-1,j]) * (hN*(hN + hS))**(-1)
            W = 2 * (condutivity[i,j] + condutivity[i,j-1]) * (hW*(hW + hE))**(-1)
            S = 2 * (condutivity[i,j] + condutivity[i+1,j]) * (hS*(hS + hN))**(-1)
            
            P = E + N + W + S 

            v_fut[i,j] += (1.0 / P) * (E*v_pas[i,j+1] + N*v_fut[i-1,j] + W*v_fut[i,j-1] + S*v_pas[i+1,j] + q[i,j] - P*v_pas[i,j]) 

    if np.max(np.abs(v_fut - v_pas)) > tol:
        break

xloc = np.linspace(0,nx,11)
xlab = np.around(xloc * hE, decimals=1)

zloc = np.linspace(0,nz,9)
zlab = np.around(zloc * hS, decimals=1)

plt.imshow(v_fut,aspect="auto")   
plt.xticks(xloc,xlab)
plt.yticks(zloc,zlab)
plt.show()
