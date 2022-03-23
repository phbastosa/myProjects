import sys
import numpy as np

def rickerGenerator(fcut,nsrc,dt):
    
    ricker = np.zeros(nsrc)
    fc = fcut/(3*np.sqrt(np.pi))
    
    s = int(nsrc/2)
    for i in range(-s,s):
        aux1 = 1.0-2.0*np.pi*pow(i*dt,2.0)*pow(fc,2.0)*pow(np.pi,2.0)
        aux2 = np.exp(-np.pi*pow(i*dt,2.0)*pow(fc,2.0)*pow(np.pi,2.0))
        ricker[i + s] = aux1 * aux2 

    return ricker

def intRickerGenerator(ricker):    
    
    intRicker = np.zeros(len(ricker))
     
    for i in range(len(ricker)):
        for j in range(i+1):
            intRicker[i] += ricker[j]

    return intRicker

dt = float(sys.argv[1])              # Parâmetro de discretização espacial
nsrc = int(sys.argv[2])              # Quantidade de amostras na fonte
fcut = float(sys.argv[3])            # Frequência de pico da fonte 

source = intRickerGenerator(rickerGenerator(fcut,nsrc,dt))

source.astype("float32",order="C").tofile(sys.argv[4])