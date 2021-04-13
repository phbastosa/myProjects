import numpy as np
import matplotlib.pyplot as plt

ng = 240        # Quantidade de sensores ativos por tiro
dg = 40         # Espaçamento entre os sensores

ns = 1000       # Quantidade de fontes na aquisição
# offset = 0    # Distância entre a fonte e o sensor mais próximo
ds = 40         # Espaçamento entre disparos  

sx = np.zeros(ng*ns)
gx = np.zeros(ng*ns)
id = np.zeros(ng*ns)

gx[:ng] = np.arange(ng) * dg
# sx[:ng] = np.ones(ng) * (ng-1) * dg + offset # Formulação para aquisição End-On 
sx[:ng] = np.ones(ng) * (ng/2) * dg  # Formulação para aquisição Split-Spread
id[:ng] = np.ones(ng)

for i in range(1,ns):
    gx[i*ng:i*ng + ng] = gx[:ng] + i*ds
    sx[i*ng:i*ng + ng] = sx[:ng] + i*ds
    id[i*ng:i*ng + ng] = i+1

cmpx = np.array([])
cmpc = np.array([])

cmps = sx - (sx - gx) / 2
for cmp in cmps:
    if cmp not in cmpx:
        cmpx = np.append(cmpx,cmp)
        cmpc = np.append(cmpc,len(np.where(cmp == cmps)[0]))

plt.figure(1,figsize=(10,7))
plt.subplot(211)
plt.scatter(gx,id)
plt.scatter(sx,id)
plt.scatter(cmpx,np.ones(len(cmpx))*ns+1)
plt.gca().invert_yaxis()
plt.xlim([-5,sx[-1]+5])
plt.gca().set_xticklabels([])
plt.title(f"Geometria com {ng*ns} traços no total")
plt.ylabel("Identificador de fontes")
plt.legend(["Receptores","Fontes","Pontos médios"],loc="lower left")

plt.subplot(212)
plt.stem(cmpx,cmpc,use_line_collection=True)
plt.xlim([-5,sx[-1]+5])
plt.grid(axis="y")
plt.xlabel("Distância [m]")
plt.ylabel("Traços por CMP")
# plt.savefig("figures/cmpTraceCount.png",dpi=200,bbox_inches="tight")
plt.show()
