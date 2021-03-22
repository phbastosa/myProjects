#!/bin/bash

##########################################################################
# User area - Fill the parameters
##########################################################################

# Model parameters 
nx=1701        # horizontal samples in models 
nz=350         # vertical samples in models
dx=10          # horizontal discretization parameter 
dz=10          # vertical discretization parameter

# Time parameters
nt=10000       # total samples in time modeling 
dt=0.0005      # temporal discretization parameter

# Cerjan damping parameters 
nabc=50        # samples in Cerjan absorbing boundary condition 
par=0.0045     # parameter to use in exponential function of damp

# Acquisition Geometry parameters
ns=1           # number of shots in modeling 
ds=25          # sources spacing
nr=320         # number of receivers in modeling
dr=25          # receivers spacing
spread=320     # active receivers per shot
offsetMin=100  # minimum offset in acquisition geometry

# Source parameters
fcut=30        # maximum frequency of source Ricker in Hz
nsrc=800       # total samples of source 

# Models filename
vpPath="model/marmousi2_vp_1701x350_dh10.bin"
vsPath="model/marmousi2_vs_1701x350_dh10.bin"
rhoPath="model/marmousi2_rho_1701x350_dh10.bin"

####################################################################### 
# Processing - running auxiliary codes to build parameters 
#######################################################################
echo "Pre-contitioning parameters:"

xrecPath="parameters/xrec.bin"; xsrcPath="parameters/xsrc.bin"
python3 auxCodes/buildEndOnGeometry.py $spread $dr $ns $offsetMin $ds $dx $xsrcPath $xrecPath  
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
echo -e "$nx\n$nz\n$nt\n$dx\n$dz\n$dt\n$nabc\n$nr\n$ns\n$spread\n$nsrc" > $parFileName
echo "Parameters was built..."

pgcc -acc -fast -ta=tesla,cc60 elasticIsotropic2D_offshore.c -lm -o run.exe
./run.exe $parFileName $xrecPath $xsrcPath $sourcePath $vp $vs $rho $dampPath 
rm run.exe

rm parameters/*.bin parameters/*.txt

transp n1=$spread n2=$nt <data/seism.bin | ximage n1=$nt d1=$dt n2=$spread d2=$dr perc=99 &
