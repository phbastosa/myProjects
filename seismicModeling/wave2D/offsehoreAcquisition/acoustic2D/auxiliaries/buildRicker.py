import sys
import numpy as np

def rickerGenerator(fcut,nsrc,dt):
    
    ricker = np.zeros(nsrc)
    fc = fcut/(3*np.sqrt(np.pi))
    
    s = int(nsrc/2)
    for i in range(-s,s):
        aux1 = 1 - 2 * pow(i*dt,2)*pow(fc,2)*pow(np.pi,2)
        aux2 = np.exp(-pow(i*dt,2)*pow(fc,2)*pow(np.pi,2))
        ricker[i + s] = aux1 * aux2 

    return ricker

def halfDerivative(wavelet, dt):
    
    fwavelet = np.fft.fft(wavelet)
    w = np.fft.fftfreq(len(wavelet),dt)
    w = (w / np.max(w)) * np.pi
    fwavelet = fwavelet * np.sqrt(complex(0,1) * w)

    return np.real(np.fft.ifft(fwavelet))

dt = float(sys.argv[1])              # Parâmetro de discretização espacial
nsrc = int(sys.argv[2])              # Quantidade de amostras na fonte
fcut = float(sys.argv[3])            # Frequencia de pico da fonte 

source = halfDerivative(rickerGenerator(fcut,nsrc,dt),dt)
np.savetxt(str(sys.argv[4]),source,fmt="%.5f")