#!/bin/bash

##########################################################################
# User area - Fill the parameters
##########################################################################

# Model parameters 
nx=2772        # horizontal samples in velocity model 
nz=240         # vertical samples in velocity model
dx=12.5        # horizontal discretization parameter 
dz=12.5        # vertical discretization parameter

# Time parameters
nt=4000         # total samples in time modeling 
dt=0.0005       # temporal discretization parameter

# Cerjan damping parameters 
nabc=50        # samples in Cerjan absorbing boundary condition 
par=0.0045     # parameter to use in exponential function of damp

# Acquisition Geometry parameters
ns=1064        # number of shots in modeling 
ds=25          # sources spacing
nr=1383        # number of receivers in modeling
dr=25          # receivers spacing
spread=320     # active receivers per shot
offsetMin=100  # minimum offset in acquisition geometry

# Source parameters
fcut=100       # maximum frequency of source Ricker in Hz
nsrc=400       # total samples of source 

modelPath="model/modelEngland_nx2772_nz240_h12.5.bin"

####################################################################### 
# Processing - running auxiliary codes to build parameters 
#######################################################################
echo "Pre-contitioning parameters:"

xrecPath="parameters/xrec.bin"; xsrcPath="parameters/xsrc.bin"
python3 auxCodes/buildEndOnGeometry.py $spread $dr $ns $offsetMin $ds $dx $xsrcPath $xrecPath  
echo -e "\nAcquisition was built..."

inputModel="model/vp_input.bin"
python3 auxCodes/buildBoundaries.py $nx $nz $nabc $modelPath $inputModel
echo -e "Model was built..."

inputDamp="parameters/damp.bin"
python3 auxCodes/buildCerjanABC.py $nx $nz $nabc $par $inputDamp
echo -e "Cerjan ABC was built..."

source="parameters/wavelet.txt"
python3 auxCodes/buildSource.py $dt $nsrc $fcut $source
echo "Wavelet was built..."

parFileName="parameters/modelingParameters.txt"
echo -e "$nx\n$nz\n$nt\n$dx\n$dz\n$dt\n$nabc\n$spread\n$nr\n$ns\n$nsrc\n" > $parFileName
echo "Parameters file for migration was built..."

pgcc -acc -fast -ta=tesla,cc60 RTM.c -lm -o rtm.exe
./rtm.exe $parFileName $inputModel $inputDamp $source $xsrcPath $xrecPath

rm rtm.exe

transp n1=$nx n2=$nz <results/outputImage.bin >results/rawImage.bin
transp n1=$nx n2=$nz <results/outputImageDsumComp.bin >results/imageDsumComp.bin
transp n1=$nx n2=$nz <results/outputImageRsumComp.bin >results/imageRsumComp.bin
transp n1=$nx n2=$nz <results/outputImageSqrtDsumComp.bin >results/imageSqrtDsumComp.bin
transp n1=$nx n2=$nz <results/outputImageSqrtRsumComp.bin >results/imageSqrtRsumComp.bin
transp n1=$nx n2=$nz <results/outputImageDRsumComp.bin >results/imageDRsumComp.bin

rm results/output*
