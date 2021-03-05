#!/bin/bash

##########################################################################
# User area - Fill the parameters
##########################################################################

# Model parameters 
nx=1700         # horizontal samples in models 
nz=351          # vertical samples in models
dx=10           # horizontal discretization parameter 
dz=10           # vertical discretization parameter

# Time parameters
nt=8000         # total samples in time modeling 
dt=0.001        # temporal discretization parameter

# Cerjan damping parameters 
nabc=50         # samples in Cerjan absorbing boundary condition 
parb=0.0055     # parameter to use in exponential function of damp

# Acquisition Geometry parameters
ns=1            # number of shots in modeling 
ds=20           # sources spacing
nr=400          # number of receivers in modeling
dr=20           # receivers spacing
spread=400      # active receivers per shot
mOffset=100     # minimum offset in acquisition geometry

nsnaps=100      # total snapshots to visualize

# Source parameters
fcut=30         # maximum frequency of Ricker source in Hz
nsrc=600        # total samples of source 

# Models filename
vpPath="model/marmousi2_vp_351x1700_dh10.bin"

# Acquisition geometry
xrecPath="parameters/xrec.bin"
xsrcPath="parameters/xsrc.bin"

####################################################################### 
# Processing - running auxiliary codes to build parameters 
#######################################################################
echo "Pre-contitioning parameters:"

xrecPath="parameters/xrec.bin"; xsrcPath="parameters/xsrc.bin" 
python3 auxiliaries/buildEndOnGeometry.py $spread $dr $ns $mOffset $ds $dx $xsrcPath $xrecPath  
echo -e "\nAcquisition was built..."

inputModel="model/vpInput.bin"
python3 auxiliaries/buildBoundaries.py $nx $nz $nabc $vpPath $inputModel
echo -e "Model was built..."

inputDamp="model/damp.bin"
python3 auxiliaries/buildCerjanABC.py $nx $nz $nabc $parb $inputDamp
echo -e "Cerjan condition was built..."

sourceFile="parameters/wavelet.bin"
python3 auxiliaries/buildSource.py $dt $nsrc $fcut $sourceFile
echo "Wavelet was built..."

parFileName="parameters/modelingParameters.txt"
echo -e "$nx\n$nz\n$nt\n$dx\n$dz\n$dt\n$nabc\n$nr\n$ns\n$nsrc\n$spread\n$nsnaps" > $parFileName
echo "Modeling parameters was built..."

gcc acoustic2D.c -lm -O3 -o run.exe
./run.exe $parFileName $inputModel $inputDamp $sourceFile $xsrcPath $xrecPath 

rm run.exe
rm model/damp.bin model/vpInput.bin 
rm parameters/* 

