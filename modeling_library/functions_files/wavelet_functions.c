# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "../header_files/wavelet_functions.h"

void mem_set(float * vector, int length)
{
    int index;
    for (index = 0; index < length; index++)
    {
        vector[index] = 0.0f;
    }
}

float * build_ricker(int wavelet_size, float f_cut, float dt) 
{
    int ii;
    float aux1,aux2;
    float pi = 4.0*atan(1.0);
    float fc = f_cut/(3 * sqrt(pi));

    float * wavelet = (float *) malloc(wavelet_size*sizeof(float));

    mem_set(wavelet, wavelet_size);   
    for(ii=-wavelet_size/2; ii<wavelet_size/2; ii++) {       
        aux1 = exp(-pi*powf(ii*dt,2.0)*powf(fc,2.0)*powf(pi,2.0));
        aux2 = 1 - 2*pi*powf(ii*dt,2.0)*powf(fc,2.0)*powf(pi,2.0);
        wavelet[ii + wavelet_size/2] = aux1 * aux2; 
    }

    return (wavelet);
}

float * build_ricker_printed(int wavelet_size, float f_cut, float dt, char * name) 
{
    int ii;
    float aux1,aux2;
    float pi = 4.0*atan(1.0);
    float fc = f_cut/(3 * sqrt(pi));

    float * wavelet = (float *) malloc(wavelet_size*sizeof(float));

    mem_set(wavelet, wavelet_size);   
    FILE *fonte = fopen(name,"wb");
    for(ii=-wavelet_size/2; ii<wavelet_size/2; ii++) {       
        aux1 = exp(-pi*powf(ii*dt,2.0)*powf(fc,2.0)*powf(pi,2.0));
        aux2 = 1 - 2*pi*powf(ii*dt,2.0)*powf(fc,2.0)*powf(pi,2.0);
        wavelet[ii + wavelet_size/2] = aux1 * aux2; 
        fwrite(&wavelet[ii + wavelet_size/2], sizeof(float),1,fonte);	                	
    }
    fclose(fonte);

    return (wavelet);
}

float * build_first_gaussian_derivative(int wavelet_size, float f_cut, float dt) 
{
    int ii;
    float aux1,aux2;
    float pi = 4.0*atan(1.0);
    float fc = f_cut/(3 * sqrt(pi));

    float * wavelet = (float *) malloc(wavelet_size*sizeof(float));

    mem_set(wavelet, wavelet_size);   
    for(ii=-wavelet_size/2; ii<wavelet_size/2; ii++) {       
        aux1 = exp(-pi*powf(ii*dt,2.0)*powf(fc,2.0)*powf(pi,2.0));
        wavelet[ii + wavelet_size/2] = -ii*dt*aux1;
    }

    return (wavelet);
}

float * build_first_gaussian_derivative_printed(int wavelet_size, float f_cut, float dt, char * name) 
{
    int ii;
    float aux1,aux2;
    float pi = 4.0*atan(1.0);
    float fc = f_cut/(3 * sqrt(pi));

    float * wavelet = (float *) malloc(wavelet_size*sizeof(float));

    mem_set(wavelet, wavelet_size);   
    FILE *fonte = fopen(name,"wb");
    for(ii=-wavelet_size/2; ii<wavelet_size/2; ii++) {       
        aux1 = exp(-pi*powf(ii*dt,2.0)*powf(fc,2.0)*powf(pi,2.0));
        wavelet[ii + wavelet_size/2] = -ii*dt*aux1;
        fwrite(&wavelet[ii + wavelet_size/2], sizeof(float),1,fonte);	                	
    }
    fclose(fonte);

    return (wavelet);
}
