import numpy as np
import matplotlib.pyplot as plt

def theoreticalReflections(R):
    
    Rs = np.zeros(len(R)) # Reflections that reach the surface
    Tl = np.zeros(len(R)) # Transmission in each layer 

    Rs[0] = R[0]          # First reflection
    Tl[0] = 1 - R[0]      # First transmission

    # Solving for each refletion after first one
    for i in range(1,len(R)):   
        Tl[i] = Tl[i-1] * (1 - R[i]) # For other layers
        Rs[i] = Tl[i-1] * R[i] # Reflections are dependent  

        for j in range(i,0,-1):  
            Rs[i] *= (1 - R[i-j]) # Solving dependencies

    return Rs # Returning theoretical reflections     

# Model dimensions

nz = 350
dz = 10.0
prof = np.arange(nz) * dz

# Model proprieties

vels = np.array([1500,1700,1900,2300,2500],dtype=float)
rhos = np.array([1000,1990,2046,2146,2192],dtype=float)
depth = np.array([150,200,250,300])

vel = np.zeros(nz)
rho = np.zeros(nz) 

# Putting proprieties in each location using corretly dimensions
 
for i in range(len(vels)):
    if i == 0:
        vel[:150] = vels[i]
        rho[:150] = rhos[i]
    else:
        vel[int(i*50 + 100):int(i*50 + 150)] = vels[i]
        rho[int(i*50 + 100):int(i*50 + 150)] = rhos[i]

# Calculation each reflectivity 

R = np.zeros(len(vels)-1)    # P reflectivity
for i in range(len(vels)-1):
    R[i] = (rhos[i+1] * vels[i+1] - rhos[i] * vels[i]) / (rhos[i+1] * vels[i+1] + rhos[i] * vels[i])

# Putting reflectivity in correctly location

Ref = np.zeros(nz)
for i in range(len(R)):
    Ref[depth[i]] = R[i]

# Calculating theoretical reflections arrived in receptors

Rs = theoreticalReflections(Ref)
Rw = theoreticalReflections(R)

# Building wavelet 

fcut = 30
nsrc = 300
dt = 0.002

t = np.arange(nsrc) * dt
ricker = np.zeros(nsrc)

s = int(nsrc/2)

fc = fcut / (3 * np.sqrt(np.pi))

for i in range(-s,s):
    aux1 = 1.0 - 2.0*np.pi*pow(i*dt,2.0)*pow(fc,2.0)*pow(np.pi,2.0)
    aux2 = np.exp(-np.pi*pow(i*dt,2.0)*pow(fc,2.0)*pow(np.pi,2.0))
    ricker[i + s] = aux1 * aux2 

# Convolving wavelet with teoretical reflectivity

nt = 1501
sz = 100
ts = np.zeros(len(depth))

time = np.arange(nt) * dt
reflectivity = np.zeros(nt)

soma = 2 * (depth[0] - sz - 1) * dz / vels[0]
reflectivity[int(soma/dt)] = Rw[0]

for i in range(1,len(ts)):
    soma += 2 * (depth[i] - depth[i-1] - 1) * dz / vels[i] 
    
    reflectivity[int(soma/dt)] = Rw[i]

trace = np.convolve(reflectivity,ricker,"same")

# Generating images in order

plt.figure(1,figsize=(25,10))

plt.subplot(151)
plt.plot(vel,prof)
plt.plot(rho,prof)
plt.gca().invert_yaxis()
plt.title("Properties Model",fontsize=15)
plt.xlabel("Values",fontsize=15)
plt.ylabel("Depth [m]",fontsize=15)
plt.legend(["Velocity [m/s]","Density [kg/m³]"],fontsize=10)

##############################

plt.subplot(152)
plt.plot(Ref,prof)
plt.plot(1000,"*")
plt.gca().invert_yaxis()
plt.title("Raw reflectivity",fontsize=15)
plt.xlabel("Coefficients",fontsize=15)
plt.legend(["$R = \dfrac{v_{i+1}\u03C1_{i+1} - v_{i}\u03C1_{i}}{v_{i+1}\u03C1_{i+1} + v_{i}\u03C1_{i}}$","Source point"],fontsize=10)
plt.xlim([-0.1,0.5])

##############################

plt.subplot(153)
plt.plot(Rs,prof)
plt.plot(1000,"*")
plt.gca().invert_yaxis()

plt.title("Ajusted reflectivity",fontsize=15)
plt.xlabel("Theoretical coefficients",fontsize=15)
plt.legend(["Interfaces","Source point"],fontsize=10)
plt.xlim([-0.1,0.5])

##############################
plt.subplot(154)
plt.plot(reflectivity,time)
plt.plot(0,"*")
plt.gca().invert_yaxis()

plt.title("Ajusted reflectivity",fontsize=15)
plt.xlabel("Theoretical coefficients",fontsize=15)
plt.ylabel("Time [s]",fontsize=15)
plt.legend(["interfaces","Source point"],fontsize=10)
plt.xlim([-0.1,0.5])

#############################

# plt.figure(5,figsize=(6,10))
# plt.plot(ricker,t)
# plt.gca().invert_yaxis()
# plt.title("Fonte Ricker - $f_{corte} = 30 Hz$",fontsize=15)
# plt.xlabel("Amplitude",fontsize=15)
# plt.ylabel("Tempo [s]",fontsize=15)
# plt.savefig("fonteRescalada.png",dpi=200,bbox_inches="tight")

#############################

plt.subplot(155)
plt.plot(trace,time)
plt.gca().invert_yaxis()
plt.title("Synthetic trace",fontsize=15)
plt.xlabel("Amplitude",fontsize=15)
plt.ylabel("Time [s]",fontsize=15)
plt.xlim([-0.3,0.5])
plt.savefig("result.png",bbox_inches="tight")
plt.show()
