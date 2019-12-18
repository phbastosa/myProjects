#!/bin/bash

gcc ../modeling_library/functions_files/io_functions.c ../modeling_library/functions_files/wavelet_functions.c ../modeling_library/functions_files/cerjan_functions.c ../modeling_library/functions_files/wave_equation_FDM_functions.c ../main_projects/general_modeling/acoustic_2D.c -lm -O3 -o acoustic_2D.exe

./acoustic_2D.exe

rm acoustic_2D.exe
