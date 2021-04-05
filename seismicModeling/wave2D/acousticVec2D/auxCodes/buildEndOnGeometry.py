import sys
import numpy as np

ng = int(sys.argv[1])      # Quantidade de sensores ativos por tiro
dg = int(sys.argv[2])      # Espaçamento entre os sensores

ns = int(sys.argv[3])      # Quantidade de fontes na aquisição
offset = int(sys.argv[4])  # Distância entre a fonte e o sensor mais próximo
ds = int(sys.argv[5])      # Espaçamento entre disparos  

dh = float(sys.argv[6])    # Parâmetro de discretização espacial da malha

sx = np.zeros(ng*ns)
gx = np.zeros(ng*ns)
id = np.zeros(ng*ns)

gx[:ng] = np.arange(ng) * dg
sx[:ng] = np.ones(ng) * (ng-1) * dg + offset
id[:ng] = np.ones(ng)

for i in range(1,ns):
    gx[i*ng:i*ng + ng] = gx[:ng] + i*ds
    sx[i*ng:i*ng + ng] = sx[:ng] + i*ds
    id[i*ng:i*ng + ng] = i+1

xsrc = np.array([ ],dtype=int)
xrec = np.array([ ],dtype=int)

for i in range(len(sx)):
    if sx[i] not in xsrc:
	    xsrc = np.append(xsrc,sx[i])
    if gx[i] not in xrec:
        xrec = np.append(xrec,gx[i])

xsrc /= dh
xrec /= dh 

xsrc.astype("int32",order="C").tofile(sys.argv[7])
xrec.astype("int32",order="C").tofile(sys.argv[8])