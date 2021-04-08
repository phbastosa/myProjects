#!/bin/bash

##########################################################################
# User area - Fill the parameters
##########################################################################

# Model parameters 
nx=3400        # horizontal samples in models 
nz=700         # vertical samples in models
dx=5.0         # horizontal discretization parameter 
dz=5.0         # vertical discretization parameter

# Time parameters
nt=10000       # total samples in time modeling 
dt=0.0005      # temporal discretization parameter

# Cerjan damping parameters 
nabc=50        # samples in Cerjan absorbing boundary condition 
par=0.0045     # parameter to use in exponential function of damp

# End-On Acquisition Geometry parameters
# ns=357          # number of shots in modeling 
# ds=25           # sources spacing
# nr=677          # number of receivers in modeling
# dr=25           # receivers spacing
# spread=320      # active receivers per shot
# mOffset=100     # minimum offset in acquisition geometry

# OBN Acquisition Geometry parameters
ns=201          # number of shots in modeling 
ds=85           # sources spacing
nr=101          # number of receivers in modeling
dr=170          # receivers spacing
wb=100          # water bottim depth

# Source parameters
fcut=50        # maximum frequency of source Ricker in Hz
nsrc=600       # total samples of source 

# Models filename
vpPath="model/marmousi2_vp_700x3400_dh5.bin"
vsPath="model/marmousi2_vs_700x3400_dh5.bin"
rhoPath="model/marmousi2_rho_700x3400_dh5.bin"

####################################################################### 
# Processing - running auxiliary codes to build parameters 
#######################################################################
echo "Pre-contitioning parameters:"

xrecPath="parameters/xrec.bin"; xsrcPath="parameters/xsrc.bin" 
python3 auxCodes/buildOBNGeometry.py $nr $dr $ns $ds $dx $xsrcPath $xrecPath  
echo -e "\nAcquisition was built..."

vp="parameters/vpInput.bin"; vs="parameters/vsInput.bin"; rho="parameters/rhoInput.bin"
python3 auxCodes/buildBoundaries.py $nx $nz $nabc $vpPath $vsPath $rhoPath $vp $vs $rho  
echo -e "Model was built..."

dampPath="parameters/damp.bin"
python3 auxCodes/buildCerjanABC.py $nx $nz $nabc $par $dampPath
echo -e "Cerjan was built..."

sourcePath="parameters/source.bin"
python3 auxCodes/buildSource.py $fcut $nsrc $dt $sourcePath
echo "Wavelet was built..."

parFileName="parameters/modelingParameters.txt"
echo -e "$nx\n$nz\n$nt\n$dx\n$dz\n$dt\n$nabc\n$nr\n$ns\n$nsrc\n$wb" > $parFileName
echo "Modeling parameters was built..."

pgcc -acc -fast -ta=tesla,cc60 elasticIsotropic2D.c -lm -o run.exe
./run.exe $parFileName $xrecPath $xsrcPath $sourcePath $vp $vs $rho $dampPath 

rm run.exe parameters/*.bin parameters/*.txt