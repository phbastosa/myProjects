import sys
import numpy as np
import matplotlib.pyplot as plt

def rickerGenerator(fcut,nsrc,dt):
    
    ricker = np.zeros(nsrc)
    fc = fcut/(3*np.sqrt(np.pi))
    
    s = int(nsrc/2)
    for i in range(-s,s):
        aux1 = 1 - 2 * np.pi*pow(i*dt,2)*pow(fc,2)*pow(np.pi,2)
        aux2 = np.exp(-np.pi*pow(i*dt,2)*pow(fc,2)*pow(np.pi,2))
        ricker[i + s] = aux1 * aux2 

    return ricker

def intRickerGenerator(ricker):    
    
    intRicker = np.zeros(len(ricker))
     
    for i in range(len(ricker)):
        for j in range(i+1):
            intRicker[i] += ricker[j]

    return intRicker

def halfDerivative(wavelet, dt):
    
    fwavelet = np.fft.fft(wavelet)
    
    w = np.fft.fftfreq(len(wavelet),dt)
    w = (w / np.max(w)) * np.pi

    fwavelet *= np.sqrt(complex(0,1) * w)
    
    return np.real(np.fft.ifft(fwavelet))

dt = float(sys.argv[1])
nsrc = int(sys.argv[2])
fcut = float(sys.argv[3])

simpGridRicker = halfDerivative(rickerGenerator(fcut,nsrc,dt),dt)
stagGridRicker = intRickerGenerator(simpGridRicker)

simpGridRicker.astype("float32",order="C").tofile(sys.argv[4])
stagGridRicker.astype("float32",order="C").tofile(sys.argv[5])

