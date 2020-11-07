#!/bin/bash

##########################################################################
# User area - Fill the parameters
##########################################################################

# Model parameters - Using marmousi 2 resized

nx=1700    # horizontal samples in velocity model 
nz=351     # vertical samples in velocity model
dx=10      # horizontal discretization parameter 
dz=10      # vertical discretization parameter
abc=100    # samples in Cerjan absorbing boundary condition 
par=0.0012 # parameter to use in exponential function of damp

modelPath="model/marmousi2_vp_351x1700_dh10.bin"

# Source parameters 

fcut=30   # frequency cutoff of source Ricker in Hz
nsrc=600  # total samples of source 
dt=0.001  # temporal discretization of modeling

# Acquisition geometry

fs=5000   # first source position in meters
ds=20     # source spacing in meters 
ns=600    # total source in modeling

frec=4900 # nearest source receptor position in meters
drec=10   # receptor spacing in meters 
nrec=490  # receptor group per shot (spread)

####################################################################### 
# Processing - running auxiliary codes to build parameters 
#######################################################################
echo "Pre-contitioning parameters:"

inputModel="model/vp_input.bin"
python3 auxiliaries/buildBoundaries.py $nx $nz $abc $modelPath $inputModel
echo -e "\nModel was built..."

inputDamp="model/damp.bin"
python3 auxiliaries/buildCerjanABC.py $nx $nz $abc $par $inputDamp
echo -e "Cerjan absorbing condition was built..."

source="parameters/wavelet.txt"
python3 auxiliaries/buildRicker.py $dt $nsrc $fcut $source
echo "Wavelet was built..."

# Streamer geometry generation
xshotsFile="parameters/xsrc.txt"
xrecpsFile="parameters/xrec.txt"
zshotsFile="parameters/zsrc.txt"
zrecpsFile="parameters/zrec.txt"
python3 auxiliaries/buildGeometry.py $dx $ns $fs $ds $nrec $frec $drec $xshotsFile $xrecpsFile $abc
echo "Acquisition geometry was built..."

parFileName="parameters/modelingParameters.txt"
echo -e "$nx\n$nz\n$nt\n$dx\n$dz\n$dt\n$abc\n$nrec\n$ns\n$nsrc\n" > $parFileName
echo "Parameters file for migration was built..."

