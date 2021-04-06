#!/bin/bash

file="../marmousi2_vp_700x3400_dh5.bin"
exit="../marmousi2_vp_700x3400_dh5_smooth.bin"

nx=3400
nz=700

times=20

gcc vagarosidade.c -o vag.exe
./vag.exe $file $nx $nz slowness.bin

smooth2 n1=$nz n2=$nx r1=5 r2=5 <slowness.bin >slsmooth.bin

for ii in $(seq 1 $times); do
    smooth2 n1=$nz n2=$nx r1=10 r2=10 <slsmooth.bin >aux.bin
    smooth2 n1=$nz n2=$nx r1=10 r2=10 <aux.bin >slsmooth.bin
done

./vag.exe slsmooth.bin $nx $nz $exit

rm vag.exe aux.bin slsmooth.bin slowness.bin
