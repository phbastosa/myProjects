# ifndef UTILS_HPP
# define UTILS_HPP

# include <cmath>
# include <chrono>
# include <string>
# include <fstream>
# include <iostream>

void linspace(float * x, float xi, float xf, int n);
void export_array(std::string path, float * array, int n);

# endif