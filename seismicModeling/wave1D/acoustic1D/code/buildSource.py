import sys
import numpy as np

dt = float(sys.argv[1])              
nsrc = int(sys.argv[2])              
fcut = float(sys.argv[3])             
source = np.zeros(nsrc)              

fc = fcut / (3 * np.sqrt(np.pi))     
s = int(nsrc/2)                     
for i in range(-s,s):              
    aux1 = 1 - 2 * pow(i*dt,2)*pow(fc,2)*pow(np.pi,2)
    aux2 = np.exp(-pow(i*dt,2)*pow(fc,2)*pow(np.pi,2))    
    source[i + s] = aux2 * aux1

source.astype("float32",order="C").tofile(sys.argv[4])
