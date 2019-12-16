# ifndef CERJAN_FUNCTIONS_H_DEFINED
# define CERJAN_FUNCTIONS_H_DEFINED

void Cerjan_1D_acoustic_attenuation(float * P_pre, float * P_fut, float * factor, int abc_layer, int Nx);

void Cerjan_2D_corners(float *up_left, float *up_right, float *down_left, 
                       float *down_right, float *factor, int n_abc);

void Cerjan_2D_acoustic_attenuation(float *P_pre, float *P_fut, float *factor, 
                                    int nxx, int nzz, int abc_layer, float *down_left, 
                                    float *down_right, float *up_left, float *up_right);

void Cerjan_2D_elastic_attenuation(float *U, float *V, float *Txx, float *Tzz,
                                   float *Txz, int nxx, int nzz, float *fator, 
                                   float *up_left, float *up_right, float *down_left, 
                                   float *down_right, int borda);

float * factor_attenuation(float parameter, int n_bondary);

# endif // CERJAN_FUNCTIONS_H_DEFINED