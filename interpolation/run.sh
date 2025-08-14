#!/bin/bash

libs="linear.cpp cubic.cpp utils.cpp"

g++ 1d_example.cpp $libs -lm -o ex1d.exe
./ex1d.exe

python3 1d_analysis.py
rm *.bin *.exe

echo -e "\n--------------\n"

g++ 2d_example.cpp $libs -lm -o ex2d.exe
./ex2d.exe

python3 2d_analysis.py
rm *.bin *.exe

echo -e "\n--------------\n"

g++ 3d_example.cpp $libs -lm -o ex3d.exe
./ex3d.exe

python3 3d_analysis.py
rm *.bin *.exe
