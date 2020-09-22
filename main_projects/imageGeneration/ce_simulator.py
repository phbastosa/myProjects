import numpy as np
import matplotlib.pyplot as plt

################# Data input ##########################
print(35*"-=")
print("Dipole-dipole survey simulation for eletrical resistivity tomography")
print(35*"-=")

lineLength = float(input("Type the line length: "))

electSpace = float(input("Type the electrode spacing unit: "))

invesLevel = int(input("Type the investigation levels: "))

################## Processing ################### 

nElectrodes = int(lineLength / electSpace + 1)

diags = (nElectrodes - (invesLevel + 2)) * invesLevel
others = sum(range(1,invesLevel))
nPoints = diags + others 

electrodesPosition = np.linspace(0,lineLength,nElectrodes)

displPointsPosition = np.array([])
depthPointsPosition = np.array([])

for i in range(invesLevel):
    displInc = 3*electSpace/2 + i*electSpace/2
    levelPoints = np.arange(displInc,lineLength - displInc + 1,electSpace) 
    displPointsPosition = np.append(displPointsPosition,levelPoints)

    depthInc = i*electSpace/2
    depth = np.ones(len(levelPoints)) * electSpace 
    depthPointsPosition = np.append(depthPointsPosition,depth+depthInc)

flag = 1
plt.figure(f"Simulation {flag}",figsize=(12,5))
plt.scatter(electrodesPosition,np.zeros(nElectrodes))
plt.scatter(displPointsPosition,depthPointsPosition)
plt.gca().invert_yaxis()

plt.title(f"Acquisition geometry - dipole-dipole survey with {int(nPoints)} points and {int(nElectrodes)} electrodes")
plt.xlabel("Displacement [m]")
plt.ylabel("Depth [m]")
plt.legend(["Electrodes position","Data position"],fontsize=12,loc="lower left")
plt.savefig("acquisition.png")
plt.show()