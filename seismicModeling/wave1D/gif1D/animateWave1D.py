import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def readbinaryfile(dim1,dim2,filename):
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='F')
    return matrix

nx = 300
nt = 3000

snapshots = readbinaryfile(nt,nx,"seism.bin")

# plt.imshow(snapshots,aspect="auto",cmap="Greys")
# plt.show()

fig, ax = plt.subplots()

ax.set_xlim(0, nx)
ax.set_ylim(-2,2)

x = np.arange(0,nx)                  # Dominio
line, = ax.plot(x,snapshots[0,:])    # Condição inicial

def init():  # Começando a imagem sempre limpa
    line.set_ydata([np.nan] * len(x))
    return line,

def animate(i): # Atualizo o dado com o passo seguinte
    line.set_ydata(snapshots[i,:])  
    return line,

ani = animation.FuncAnimation(fig,animate,init_func=init,interval=0.5,frames=range(0,nt,10),blit=True,repeat=True)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Paulo'), bitrate=1000)

ani.save('wave1D.mp4', writer=writer)

plt.show()