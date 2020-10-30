#/bin/bash

# Model Parameters

nz=400
dz=10

nabc=50
parb=0.0045

modelPar="parameters/modelNodes.txt"
modelPath="model/velocities.bin"
dampPath="model/damp.bin"

nodes=('   0' ' 100' ' 120' ' 140' ' 160' ' 180' ' 200' ' 220' ' 240' ' 280' ' 320' ' 360')
veloc=('1500' '1600' '1700' '2000' '2300' '2500' '1800' '2400' '2600' '2800' '3000' '3200')

if [ -f $modelPar ]; then 
    rm $modelPar 
fi

for index in $(seq 0 ${#nodes[@]}) 
do  
    echo ${nodes[index]} ${veloc[index]} >> $modelPar 
done

python3 code/buildModel.py $modelPar $nz $nabc $parb $modelPath $dampPath

# Source parameters

nsrc=500
dt=0.0005
fcut=30

sourcePath="parameters/source.bin"
python3 code/buildSource.py $dt $nsrc $fcut $sourcePath

# Modeling parameters

nt=7001

parFileName="parameters/modelingParameters.txt"
echo -e "$nz\n$nt\n$dz\n$dt\n$nabc\n$nsrc\n" > $parFileName

gcc code/acoustic1D.c -lm -o wave1D.exe 
./wave1D.exe $parFileName $modelPath $dampPath $sourcePath 
rm wave1D.exe