import numpy as np
import matplotlib.pyplot as plt

nx = 300                # Definindo o tamanho da corda
dx = 5.0                # Definindo o tamanho do espaço entru um ponto e outro
v = np.ones(nx) * 2500  # Definindo a velocidade, equivalente à tensão na corda

nt = 3000               # Total de amostras no tempo
dt = 0.0005             # Discretização temporal 

# Definindo a fonte sísmica
nsrc = 1200             # Quantidade de amostras na fonte
freq = 20.0             # Frequencia de pico da fonte 
source = np.zeros(nsrc) # Array que armazenará os valores da fonte

s = int(nsrc/2)         # Parâmetro para montar wavelet fase zero

for ii in range(-s,s):  # Calculo da wavelet Ricker
    aux1 = np.exp(-np.pi**2 * (ii*dt)**2 * freq**2)
    aux2 = 1.0 - (2.0 * np.pi**2 * (ii*dt)**2 * freq**2)
    source[ii + s] = aux2 * aux1

# Visualizando a wavelet 
# plt.plot(np.arange(nsrc)*dt,source)
# plt.show()

# Definindo campos de onda
u_pas = np.zeros(nx)    # Campo de onda do tempo passado    
u_pre = np.zeros(nx)    # Campo de onda no tempo presente
u_fut = np.zeros(nx)    # Campo de onda no tempo futuro 

seismogram = np.zeros((nt,nx)) # Sismograma registra os campos de onda em todos os instantes de tempo

for tt in range(nt):

    if tt % 500 == 0:   # Flag para 
        print(f"Tempo de propagação: {tt*dt:0.2f}")

    if tt < nsrc:       # Disparando a fonte no modelo
        u_pre[0] += source[tt]

    for ii in range(1,nx-1): # Calculo do campo de onda de tempo futuro com diferenças finitas
        aux1 = (u_pre[ii + 1] - 2.0*u_pre[ii] + u_pre[ii - 1]) / dx**2.0
        aux2 = dt**2.0 * v[ii]**2.0
        aux3 = 2.0*u_pre[ii] - u_pas[ii]

        u_fut[ii] = aux1 * aux2 + aux3

    u_pas[:] = u_pre[:]  # Atualizando os campos de onda para o calculo no tempo seguinte 
    u_pre[:] = u_fut[:]

    seismogram[tt,:] = u_fut # Atribuindo o campo de onda calculado ao sismograma do tempo correspondente

inst = 500       # Instante de tempo capturado no snapshot

# Visualização da propagação da onda
plt.figure(1)
plt.title(f"Propagação da onda na corda num instante t = {inst*dt} s")
plt.xlabel("Corda [m]")
plt.ylabel("Amplitude da onda")
plt.plot(seismogram[inst,:])

plt.figure(2)
plt.title("Propagação da onda em todos os instantes de tempo")
plt.xlabel("Corda [m]")
plt.ylabel("Tempo de propagação [s]")

# Alterando o eixo de tamanho de amostra para tempo
locks = np.linspace(0,nt,9)
labels = []
for tt in locks:
    labels.append(str(format(tt*dt)))

# Aplicando novas labels na imagem
plt.yticks(locks,labels)

plt.imshow(seismogram,aspect="auto",cmap="Greys")
cbar = plt.colorbar()
cbar.set_label("Amplitude da onda")
plt.show()

# print(np.shape(seismogram))
# np.savetxt("snapshots.txt",np.reshape(seismogram,[nx*nt]))
