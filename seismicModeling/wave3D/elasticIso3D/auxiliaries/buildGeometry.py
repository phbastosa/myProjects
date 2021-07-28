import sys
import numpy as np

nx = int(sys.argv[1])
ny = int(sys.argv[2])

nrecx = int(sys.argv[3])
nrecy = int(sys.argv[4])

abc = int(sys.argv[5])

dsx = nx / (nrecx - 1)
dsy = ny / (nrecy - 1)

recx, recy = np.meshgrid(np.linspace(abc,nx+abc,nrecx,dtype=int),np.linspace(abc,ny+abc,nrecy,dtype=int))
srcx, srcy = np.meshgrid(np.linspace(abc,nx+abc,nrecx,dtype=int),np.linspace(abc,ny+abc,nrecy,dtype=int))

xrec = np.reshape(recx,[nrecx*nrecy],order="C")
yrec = np.reshape(recy,[nrecx*nrecy],order="C")
xsrc = np.reshape(srcx,[nrecx*nrecy],order="C")
ysrc = np.reshape(srcy,[nrecx*nrecy],order="C")

xrec.astype("int32",order="C").tofile(sys.argv[6])
yrec.astype("int32",order="C").tofile(sys.argv[7])
xsrc.astype("int32",order="C").tofile(sys.argv[8])
ysrc.astype("int32",order="C").tofile(sys.argv[9])
