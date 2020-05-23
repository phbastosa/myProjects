import numpy as np
import matplotlib.pyplot as plt

ns = 1000
dt = 0.00012
t = np.arange(ns) * dt

def readbinaryfile(dim1,filename):
    with open(filename, 'rb') as f:    
        data = np.fromfile(filename, dtype = np.float32, count = dim1)
    return data
    
fgd = readbinaryfile(ns,"first_gaussian_derivative.bin")    

plt.figure(figsize=(10,5),dpi=150)
plt.plot(t,fgd)
plt.xlabel("Time in seconds",fontsize=15)
plt.ylabel("Amplitude",fontsize=15)
plt.savefig("source_data.png")
plt.show()

freqs = np.fft.fftfreq(ns,dt)
wfgd = np.abs(np.fft.fft(fgd))

plt.figure(figsize=(10,5),dpi=150)
plt.stem(freqs,wfgd,use_line_collection=True)
plt.xlabel("Frequency in Hertz",fontsize=15)
plt.ylabel("Amplitude",fontsize=15)
plt.xlim([0,210])
plt.savefig("spectrum_data.png")
plt.show()

rickers = []
for ii in range(30,76,15):
    rickers.append("ricker_" + str(ii) + "Hz.bin")   
    
r30 = readbinaryfile(ns,rickers[0])
r45 = readbinaryfile(ns,rickers[1])
r60 = readbinaryfile(ns,rickers[2])
r75 = readbinaryfile(ns,rickers[3])

plt.figure(figsize=(10,5),dpi=150)
plt.plot(t,r30,t,r45,t,r60,t,r75)
plt.xlabel("Time in seconds",fontsize=15)
plt.ylabel("Amplitude",fontsize=15)
plt.legend(["Ricker 30Hz","Ricker 45Hz","Ricker 60Hz","Ricker 75Hz"],fontsize=15)
plt.savefig("sources_mig.png")
plt.show()

wr30 = np.abs(np.fft.fft(r30))
wr45 = np.abs(np.fft.fft(r45))
wr60 = np.abs(np.fft.fft(r60))
wr75 = np.abs(np.fft.fft(r75))

plt.figure(figsize=(10,5),dpi=150)
plt.plot(freqs,wr30,freqs,wr45,freqs,wr60,freqs,wr75)
plt.xlabel("Frequency in Hertz",fontsize=15)
plt.ylabel("Amplitude",fontsize=15)
plt.legend(["Ricker 30Hz","Ricker 45Hz","Ricker 60Hz","Ricker 75Hz"],fontsize=15)
plt.xlim([0,210])
plt.savefig("spectrum_mig.png")
plt.show()
    
