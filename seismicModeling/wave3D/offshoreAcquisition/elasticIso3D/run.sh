#!/bin/bash

##########################################################################
# User area - Fill the parameters
##########################################################################

# Modeling parameters 

nx=200     # horizontal samples in velocity model 
ny=200     # horizontal samples in velocity model 
nz=200     # vertical samples in velocity model
dx=10      # horizontal discretization parameter 
dy=10      # horizontal discretization parameter
dz=10      # vertical discretization parameter

nt=2000    # Samles in time domain
dt=0.001   # temporal discretization of modeling
abc=50     # samples in Cerjan absorbing boundary condition 
par=0.0045 # parameter to use in exponential function of damp

referenceModel="model/planoParalelo.txt"

vpModelPath="model/vpModel.bin"
vsModelPath="model/vsModel.bin"
rhoModelPath="model/rhoModel.bin"

# Source parameters 

fcut=30    # frequency cutoff of source Ricker in Hz
nsrc=500   # total samples of source 

# Acquisition geometry - Nodes configuration, shots between nodes

horizon=50 # Water bottom horizon location in samples
nrecx=31   # number of receptors in x direction
nrecy=51   # number of receptors in y direction

####################################################################### 
# Processing - running auxiliary codes to build parameters 
######################################################################
# echo "Pre-contitioning parameters:"

# vpInput="model/vpInput.bin"; vsInput="model/vsInput.bin"; rhoInput="model/rhoInput.bin";
# python3 auxiliaries/buildModel.py $referenceModel $nx $ny $nz $vpModelPath $vsModelPath $rhoModelPath
# python3 auxiliaries/buildBoundaries.py $nx $ny $nz $abc $vpModelPath $vsModelPath $rhoModelPath $vpInput $vsInput $rhoInput
# echo -e "\nModel was built..."

# damp="model/damp3D.bin"
# python3 auxiliaries/buildCerjanABC.py $nx $ny $nz $abc $par $damp
# echo -e "Cerjan absorbing condition was built..."

# sourceInput="parameters/source.bin"
# python3 auxiliaries/buildSource.py $dt $nsrc $fcut $sourceInput
# echo "Wavelet was built..."

# xsrc="parameters/xsrcPositions.bin"
# ysrc="parameters/ysrcPositions.bin"
# xrec="parameters/xrecPositions.bin"
# yrec="parameters/yrecPositions.bin"
# python3 auxiliaries/buildGeometry.py $nx $ny $nrecx $nrecy $abc $xrec $yrec $xsrc $ysrc
# echo "Acquisition geometry was built..."

# parFileName="parameters/modelingParameters.txt"
# echo -e "$nx\n$ny\n$nz\n$nt\n$dx\n$dy\n$dz\n$dt\n$abc\n$nrecx\n$nrecy\n$nsrc\n$horizon" > $parFileName
# echo -e "Parameters file for modeling was built...\n\n"

# pgcc -acc -fast -ta=tesla,cc60 elasticIsotropic3D.c -lm -o run.exe
# ./run.exe $parFileName $vpInput $vsInput $rhoInput $damp $sourceInput $xsrc $ysrc $xrec $yrec
# rm run.exe

xplain=15
yplain=25
python3 results/viewSeismograms.py $nrecx $nrecy $nt $xplain $yplain $dt

