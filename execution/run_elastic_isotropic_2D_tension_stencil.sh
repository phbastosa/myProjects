#!/bin/bash

gcc ../modeling_library/functions_files/io_functions.c ../modeling_library/functions_files/wavelet_functions.c ../modeling_library/functions_files/cerjan_functions.c ../modeling_library/functions_files/wave_equation_FDM_functions.c ../main_projects/general_modeling/elastic_isotropic_2D_tension_stencil.c -lm -O3 -o elastic_isotropic_tension_stencil.exe

./elastic_isotropic_tension_stencil.exe

rm elastic_isotropic_tension_stencil.exe
