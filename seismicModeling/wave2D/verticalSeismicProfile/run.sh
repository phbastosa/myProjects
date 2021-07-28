#!/bin/bash

##########################################################################
# User area - Fill the parameters
##########################################################################

# Model parameters 
nx=400        # horizontal samples in models 
nz=600         # vertical samples in models
dx=5.0         # horizontal discretization parameter 
dz=5.0         # vertical discretization parameter

# Time parameters
nt=2001       # total samples in time modeling 
dt=0.001      # temporal discretization parameter

# Cerjan damping parameters 
nabc=100        # samples in Cerjan absorbing boundary condition 
par=0.0025      # parameter to use in exponential function of damp

# VSP Acquisition Geometry parameters
xsrc=200
zsrc=50

# Model parameters
depth=[100,150,200,250,300,350,400,450,500,550,600] 
velocities=[1500,1580,1660,1700,1790,1870,1990,2140,2300,2430,2580]

# Source parameters
fcut=50        # maximum frequency of source Ricker in Hz
nsrc=600       # total samples of source 

# Models filename
vpPath="model/planoParalelo_vp_${nz}x${nx}_dh5.bin"
vsPath="model/planoParalelo_vs_${nz}x${nx}_dh5.bin"
rhoPath="model/planoParalelo_rho_${nz}x${nx}_dh5.bin"

####################################################################### 
# Processing - running auxiliary codes to build parameters 
#######################################################################
echo "Pre-contitioning parameters:"

python3 auxCodes/buildModel.py $depth $velocities $nx $nz $vpPath $vsPath $rhoPath
echo -e "\nProperties model was built..."

vp="parameters/vpInput.bin"; vs="parameters/vsInput.bin"; rho="parameters/rhoInput.bin"
python3 auxCodes/buildBoundaries.py $nx $nz $nabc $vpPath $vsPath $rhoPath $vp $vs $rho  
echo -e "Expanded model was built..."

dampPath="parameters/damp.bin"
python3 auxCodes/buildCerjanABC.py $nx $nz $nabc $par $dampPath
echo -e "Cerjan was built..."

sourcePath="parameters/source.bin"
python3 auxCodes/buildSource.py $fcut $nsrc $dt $sourcePath
echo "Wavelet was built..."

parFileName="parameters/modelingParameters.txt"
echo -e "$nx\n$nz\n$nt\n$dx\n$dz\n$dt\n$nabc\n$xsrc\n$zsrc\n$nsrc" > $parFileName
echo "Modeling parameters was built..."

pgcc -acc -fast -ta=tesla,cc60 elasticIsotropic2D.c -lm -o run.exe
./run.exe $parFileName $sourcePath $vp $vs $rho $dampPath 

# rm run.exe parameters/*.bin parameters/*.txt
