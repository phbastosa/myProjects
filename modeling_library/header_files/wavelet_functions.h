#ifndef WAVELET_FUNCTIONS_H_DEFINED
#define WAVELET_FUNCTIONS_H_DEFINED

void mem_set(float * vector, int length);
float * build_ricker(int wavelet_size, float f_cut, float dt);
float * build_ricker_printed(int wavelet_size, float f_cut, float dt, char * name);
float * build_first_gaussian_derivative(int wavelet_size, float f_cut, float dt);
float * build_first_gaussian_derivative_printed(int wavelet_size, float f_cut, float dt, char * name);

#endif  // WAVELET_FUNCTIONS_H_DEFINED