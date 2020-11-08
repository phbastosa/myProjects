#!/bin/bash

##########################################################################
# User area - Fill the parameters
##########################################################################

# Model parameters - Using marmousi 2 resized

nx=1700    # horizontal samples in velocity model 
nz=351     # vertical samples in velocity model
nt=5000    # 
dx=10      # horizontal discretization parameter 
dz=10      # vertical discretization parameter
dt=0.001   # temporal discretization of modeling
abc=100    # samples in Cerjan absorbing boundary condition 
par=0.0005 # parameter to use in exponential function of damp

modelPath="model/marmousi2_vp_351x1700_dh10.bin"

# Source parameters 

fcut=30    # frequency cutoff of source Ricker in Hz
nsrc=600   # total samples of source 

# Acquisition geometry

fs=5000    # first source position in meters
ds=20      # source spacing in meters 
ns=600     # total source in modeling

frec=4900  # nearest source receptor position in meters
drec=-10   # receptor spacing in meters 
nrec=491   # receptor group per shot (spread)

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

sourceFile="parameters/wavelet.txt"
python3 auxiliaries/buildSource.py $dt $nsrc $fcut $sourceFile
echo "Wavelet was built..."

# Streamer geometry generation
xshotsFile="parameters/xsrc.txt"
xrecpsFile="parameters/xrec.txt"
zshotsFile="parameters/zsrc.txt"
zrecpsFile="parameters/zrec.txt"
python3 auxiliaries/buildGeometry.py $dx $ns $fs $ds $nrec $frec $drec $abc $xshotsFile $zshotsFile $xrecpsFile $zrecpsFile 
echo "Acquisition geometry was built..."

parFileName="parameters/modelingParameters.txt"
echo -e "$nx\n$nz\n$nt\n$dx\n$dz\n$dt\n$abc\n$nrec\n$ns\n$nsrc" > $parFileName
echo "Parameters file for modeling was built..."

gcc acoustic2D.c -lm -O3 -o run.exe
./run.exe $parFileName $inputModel $inputDamp $xshotsFile $zshotsFile $xrecpsFile $zrecpsFile $sourceFile 
rm run.exe