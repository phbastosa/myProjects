#!/bin/bash

##########################################################################
# User area - Fill the parameters
##########################################################################

# Model parameters 
nx=3400         # horizontal samples in models 
nz=700          # vertical samples in models
dx=5.0          # horizontal discretization parameter 
dz=5.0          # vertical discretization parameter

# Time parameters
nt=10000        # total samples in time modeling 
dt=0.0005       # temporal discretization parameter

# Cerjan damping parameters 
nabc=50         # samples in Cerjan absorbing boundary condition 
parb=0.0045     # parameter to use in exponential function of damp

# End-On Acquisition Geometry parameters
# ns=357          # number of shots in modeling 
# ds=25           # sources spacing
# nr=677          # number of receivers in modeling
# dr=25           # receivers spacing
# spread=320      # active receivers per shot
# mOffset=100     # minimum offset in acquisition geometry

# OBN Acquisition Geometry parameters
ns=200          # number of shots in modeling 
ds=85           # sources spacing
nr=101          # number of receivers in modeling
dr=170          # receivers spacing
wb=100          # water bottim depth

# Source parameters
fcut=50         # maximum frequency of Ricker source in Hz
nsrc=600        # total samples of source 

# Models filename
vpPath="model/marmousi2_vp_700x3400_dh5.bin"

# Acquisition geometry
xrecPath="parameters/xrec.bin"
xsrcPath="parameters/xsrc.bin"

####################################################################### 
# Processing - running auxiliary codes to build parameters 
#######################################################################
echo "Pre-contitioning parameters:"

xrecPath="parameters/xrec.bin"; xsrcPath="parameters/xsrc.bin" 
python3 auxCodes/buildOBNGeometry.py $nr $dr $ns $ds $dx $xsrcPath $xrecPath  
echo -e "\nAcquisition was built..."

inputModel="model/vpInput.bin"
python3 auxCodes/buildBoundaries.py $nx $nz $nabc $vpPath $inputModel
echo -e "Model was built..."

inputDamp="model/damp.bin"
python3 auxCodes/buildCerjanABC.py $nx $nz $nabc $parb $inputDamp
echo -e "Cerjan condition was built..."

sourceFile="parameters/wavelet.bin"
python3 auxCodes/buildSource.py $dt $nsrc $fcut $sourceFile
echo "Wavelet was built..."

parFileName="parameters/modelingParameters.txt"
echo -e "$nx\n$nz\n$nt\n$dx\n$dz\n$dt\n$nabc\n$nr\n$ns\n$nsrc\n$wb" > $parFileName
echo "Modeling parameters was built..."

pgcc -fast -ta=tesla,cc60 acoustic2D.c -lm -o run.exe
./run.exe $parFileName $inputModel $inputDamp $sourceFile $xsrcPath $xrecPath 

rm run.exe model/damp.bin model/vpInput.bin parameters/*.bin parameters/*.txt




