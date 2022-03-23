import sys
import numpy as np
import matplotlib.pyplot as plt

def rickerGenerator(fcut,nsrc,dt):
    
    ricker = np.zeros(nsrc)
    fc = fcut/(3*np.sqrt(np.pi))
    
    s = int(nsrc/2)
    for i in range(-s,s):
        aux1 = 1 - 2 * np.pi * pow(i*dt,2)*pow(fc,2)*pow(np.pi,2)
        aux2 = np.exp(-np.pi * pow(i*dt,2)*pow(fc,2)*pow(np.pi,2))
        ricker[i + s] = aux1 * aux2 

    return ricker

def halfDerivative(wavelet, dt):

    i = complex(0,1)
    f = np.fft.fftfreq(len(wavelet),dt)
    w = f
    
    w = (w / np.max(w)) * np.pi
    
    # fwavelet = np.fft.fft(wavelet) * np.sqrt(2.0 * np.pi * i * w)
    fwavelet = np.fft.fft(wavelet) * np.sqrt(i * w)

    parameter = 1.0 / np.max(np.real(np.fft.ifft(fwavelet)))
    parameter = 1.0 

    return parameter * np.real(np.fft.ifft(fwavelet)) 

dt = 0.001              # Parâmetro de discretização espacial
nsrc = 300              # Quantidade de amostras na fonte
fcut = 30            # Frequencia de pico da fonte 

ricker = rickerGenerator(fcut,nsrc,dt)
# source = halfDerivative(rickerGenerator(fcut,nsrc,dt),dt)
# source.astype("float32",order="C").tofile(sys.argv[4])

plt.plot(ricker)
plt.show()

# times = np.arange(nsrc) * dt
# freqs = np.fft.fftfreq(nsrc,dt)
# fftricker = np.fft.fft(ricker)
# fftsource = np.fft.fft(source) 

# plt.figure(1,figsize=(15,10))
# plt.subplot(221)
# plt.plot(times,ricker)
# plt.title("Ricker")
# plt.ylabel("Amplittude")

# plt.subplot(222)
# plt.plot(freqs,np.abs(fftricker))
# plt.xlim([0,fcut])
# plt.title("Espectro da Ricker")
# plt.ylabel("Amplittude")

# plt.subplot(223)
# plt.plot(times,source)
# plt.title("Fonte aplicada")
# plt.xlabel("Tempo [s]")
# plt.ylabel("Amplittude")

# plt.subplot(224)
# plt.plot(freqs,np.abs(fftsource))
# plt.xlim([0,fcut])
# plt.title("Espectro da Fonte aplicada")
# plt.xlabel("Frequência [Hz]")
# plt.ylabel("Amplittude")

# plt.savefig("fontes.png",dpi=200,bbox_inches="tight")
# plt.show()