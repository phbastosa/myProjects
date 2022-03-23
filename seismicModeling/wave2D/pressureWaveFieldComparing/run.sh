#!/bin/bash

##########################################################################
# User area - Fill the parameters
##########################################################################

# Model parameters 
nx=400           # horizontal samples in models 
nz=350           # vertical samples in models
dx=10.0           # horizontal discretization parameter 
dz=10.0           # vertical discretization parameter

# Time parameters
nt=1000           # total samples in time modeling 
dt=0.0008         # temporal discretization parameter

# Simple shot acquisition
xsrc=199         #
zsrc=174         #

zrec=119         #

# Source parameters
fcut=10          # maximum frequency of Ricker source in Hz
nsrc=1000         # total samples of source 

####################################################################### 
# Processing - running auxiliary codes to build parameters 
#######################################################################
echo "Pre-contitioning parameters:"

simpGridRicker="parameters/simpGridRicker.bin"
stagGridRicker="parameters/ricker_corr_10.bin"

# python3 auxCodes/sourceGenerator.py $dt $nsrc $fcut $simpGridRicker $stagGridRicker
# echo -e "\nWavelet was built..."

parFileName="parameters/modelingParameters.txt"
echo -e "$nx\n$nz\n$nt\n$dx\n$dz\n$dt\n$nsrc\n$xsrc\n$zsrc\n$zrec" > $parFileName
echo "Modeling parameters was built..."

# pgcc -acc -fast -ta=tesla,cc60 basicModelingCodes/basicAcoustic2D.c -lm -o run.exe
# ./run.exe $parFileName $simpGridRicker

# pgcc -acc -fast -ta=tesla,cc60 basicModelingCodes/basicAcousticVec2D.c -lm -o run.exe
# ./run.exe $parFileName $stagGridRicker

# pgcc -acc -fast -ta=tesla,cc60 basicModelingCodes/basicElasticIsotropic2D.c -lm -o run.exe
gcc basicModelingCodes/basicElasticIsotropic2D.c -lm -o run.exe
./run.exe $parFileName $stagGridRicker

# timeCut=1200
# python3 comparingCode/compareSeismograms.py $nx $dx $nt $dt $timeCut $nsrc
 
rm run.exe