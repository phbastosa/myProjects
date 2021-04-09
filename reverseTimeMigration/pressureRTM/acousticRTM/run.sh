#!/bin/bash

##########################################################################
# User area - Fill the parameters
##########################################################################

# Model parameters 
nx=3400        # horizontal samples in velocity model 
nz=700         # vertical samples in velocity model
dx=5.0         # horizontal discretization parameter 
dz=5.0         # vertical discretization parameter

# Time parameters
nt=10000        # total samples in time modeling 
dt=0.0005       # temporal discretization parameter

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

modelPath="model/marmousi2_vp_700x3400_dh5_smooth.bin"
dataPath="../../dataSets/acoustic_marmousi2_dh5_dataset.bin"

####################################################################### 
# Processing - running auxiliary codes to build parameters 
#######################################################################
echo "Pre-contitioning parameters:"

xrecPath="parameters/xrec.bin"; xsrcPath="parameters/xsrc.bin" 
python3 auxCodes/buildOBNGeometry.py $nr $dr $ns $ds $dx $xsrcPath $xrecPath  
echo -e "\nAcquisition was built..."

inputModel="parameters/inputModel.bin"
python3 auxCodes/buildBoundaries.py $nx $nz $nabc $modelPath $inputModel
echo -e "Model was built..."

inputDamp="parameters/damp.bin"
python3 auxCodes/buildCerjanABC.py $nx $nz $nabc $par $inputDamp
echo -e "Cerjan ABC was built..."

source="parameters/wavelet.bin"
python3 auxCodes/buildSource.py $dt $nsrc $fcut $source
echo "Wavelet was built..."

parFileName="parameters/modelingParameters.txt"
echo -e "$nx\n$nz\n$nt\n$dx\n$dz\n$dt\n$nabc\n$nr\n$ns\n$nsrc\n$wb" > $parFileName
echo "RTM parameters was built..."

pgcc -acc -fast -ta=tesla,cc60 RTM.c -lm -o rtm.exe
./rtm.exe $parFileName $inputModel $inputDamp $source $xsrcPath $xrecPath $dataPath

transp n1=$nx n2=$nz <results/outputImage.bin >results/rawImage.bin
transp n1=$nx n2=$nz <results/outputImageDsumComp.bin >results/imageDsumComp.bin
transp n1=$nx n2=$nz <results/outputImageRsumComp.bin >results/imageRsumComp.bin
transp n1=$nx n2=$nz <results/outputImageSqrtDsumComp.bin >results/imageSqrtDsumComp.bin
transp n1=$nx n2=$nz <results/outputImageSqrtRsumComp.bin >results/imageSqrtRsumComp.bin
transp n1=$nx n2=$nz <results/outputImageDRsumComp.bin >results/imageDRsumComp.bin

rm results/output* rtm.exe parameters/*.txt parameters/*.bin 