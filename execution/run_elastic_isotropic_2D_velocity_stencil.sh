#!/bin/bash

gcc ../modeling_library/functions_files/io_functions.c ../modeling_library/functions_files/wavelet_functions.c ../modeling_library/functions_files/cerjan_functions.c ../modeling_library/functions_files/wave_equation_FDM_functions.c ../main_projects/general_modeling/elastic_isotropic_2D_velocity_stencil.c -lm -O3 -o elastic_isotropic_velocity_stencil.exe

./elastic_isotropic_velocity_stencil.exe

rm elastic_isotropic_velocity_stencil.exe
