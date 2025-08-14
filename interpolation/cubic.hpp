# ifndef CUBIC_INTERPOLATION_HPP
# define CUBIC_INTERPOLATION_HPP

float cubic1d(float P[4], float dx);
float cubic2d(float P[4][4], float dx, float dy);
float cubic3d(float P[4][4][4], float dx, float dy, float dz);

# endif