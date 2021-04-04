#!/bin/bash

##########################################################################
# User area - Fill the parameters
##########################################################################

# Model parameters 
nx=500           # horizontal samples in models 
nz=500           # vertical samples in models
dx=5.0           # horizontal discretization parameter 
dz=5.0           # vertical discretization parameter

# Time parameters
nt=2000          # total samples in time modeling 
dt=0.001         # temporal discretization parameter

# Simple shot acquisition
xsrc=250         #
zsrc=300         #
zrec=100         #

# Source parameters
fcut=30          # maximum frequency of Ricker source in Hz
nsrc=400         # total samples of source 

####################################################################### 
# Processing - running auxiliary codes to build parameters 
#######################################################################
echo "Pre-contitioning parameters:"

simpGridRicker="parameters/simpGridRicker.bin"
stagGridRicker="parameters/stagGridRicker.bin"

python3 auxCodes/sourceGenerator.py $dt $nsrc $fcut $simpGridRicker $stagGridRicker
echo -e "\nWavelet was built..."

parFileName="parameters/modelingParameters.txt"
echo -e "$nx\n$nz\n$nt\n$dx\n$dz\n$dt\n$nsrc\n$xsrc\n$zsrc\n$zrec" > $parFileName
echo "Modeling parameters was built..."

pgcc -acc -fast -ta=tesla,cc60 basicModelingCodes/basicAcoustic2D.c -lm -o run.exe
./run.exe $parFileName $simpGridRicker

pgcc -acc -fast -ta=tesla,cc60 basicModelingCodes/basicAcousticVec2D.c -lm -o run.exe
./run.exe $parFileName $stagGridRicker

pgcc -acc -fast -ta=tesla,cc60 basicModelingCodes/basicElasticIsotropic2D.c -lm -o run.exe
./run.exe $parFileName $stagGridRicker

timeCut=1200
python3 comparingCode/compareSeismograms.py $nx $dx $nt $dt $timeCut $nsrc
 
rm run.exe