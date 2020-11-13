import numpy as np
import matplotlib.pyplot as plt

z = 10                                    # Profundidade da interface
rho1 = 5000                               # Resistividade da primeira camada
rhos = np.array([2500,1000,500,100,50])   # Resistividade variada para a segunda camada
a = np.linspace(1,500,1000)               # Espaçamento dos eletrodos para os arranjos
nmax = 1000                               # Simulando o infinito  

for i in range(len(rhos)): 

    # Coeficiente de reflexão
    k = (rhos[i] - rho1) / (rhos[i] + rho1) 

    # Cálculo das resistividades aparentes usando arranjo Wenner
    # soma = 0
    # for n in range(1,nmax):
    #     soma += (k**n * (1/np.sqrt(1 + (2 * n * z/a)**2) - 1/np.sqrt(4 + (2 * n * z/a)**2)))

    # rhoa = rho1*(1 + 4 * soma)

    # Cálculo das resistividades aparentes usando arranjo Schlumberger
    soma = 0
    for n in range(1,nmax):
        soma += (k**n / (1 + (2 * n * (z/a))**2)**1.5)
    
    rhoa = rho1*(1 + 2 * soma)    
    
    plt.loglog(a,rhoa)

labels = ['2500 $\Omega m$','1000 $\Omega m$','500 $\Omega m$','100 $\Omega m$','50 $\Omega m$']

plt.ylim([10,10000])
# plt.title("Curva de resistividade aparente arranjo Wenner")
plt.title("Curva de resistividade aparente arranjo Schlumberger")
# plt.xlabel("Abertura fixa dos eletrodos $m$")
plt.xlabel("Abertura dos eletrodos de corrente [$m$]")
plt.ylabel("Resistividade aparente [$\Omega m$]")
plt.figlegend(labels,loc='lower right', bbox_to_anchor=(0.35, 0.15))
plt.show()