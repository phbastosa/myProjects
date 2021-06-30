import sys
import numpy as np

ng = int(sys.argv[1])      # Quantidade de sensores ativos por tiro
dg = int(sys.argv[2])      # Espaçamento entre os sensores

ns = int(sys.argv[3])      # Quantidade de fontes na aquisição
ds = int(sys.argv[4])      # Espaçamento entre disparos  

dh = float(sys.argv[5])    # Parâmetro de discretização espacial da malha

sx = np.arange(0,ns*ds,ds)
gx = np.arange(0,ng*dg,dg) 

xsrc = sx / dh
xrec = gx / dh 

xsrc.astype("int32",order="C").tofile(sys.argv[6])
xrec.astype("int32",order="C").tofile(sys.argv[7])