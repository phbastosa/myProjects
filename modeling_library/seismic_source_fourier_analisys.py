import numpy as np
import matplotlib.pyplot as plt

nsrc = 1000
dt = 0.0002
f_cut = 100

ricker = np.zeros(nsrc)
t = np.arange(0,nsrc)*dt

fc = f_cut/(3*np.sqrt(np.pi))
for i in range(-int(nsrc/2),int(nsrc/2)):
    aux1 = np.e**(-np.pi*pow(i*dt,2.0)*pow(fc,2.0)*pow(np.pi,2.0));
    aux2 = 1 - 2*np.pi*pow(i*dt,2.0)*pow(fc,2.0)*pow(np.pi,2.0);
    ricker[i + int(nsrc/2)] = aux1 * aux2

plt.plot(t,ricker)
plt.show()

freq = np.fft.fftfreq(nsrc,dt)
ricker_fft = np.abs(np.fft.fft(ricker))

plt.plot(freq,ricker_fft)
plt.xlim([0,100])
plt.show()
