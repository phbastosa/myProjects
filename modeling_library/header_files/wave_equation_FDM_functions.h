# ifndef FDM_FUNCTIONS_H_DEFINED
# define FDM_FUNCTIONS_H_DEFINED

void update_wavefield(float * P_pas, float * P_pre, float * P_fut, int length);

void acoustic_1D_2E2T(float * P_pas, float * P_pre,float * P_fut, float * Vp, int Nx, float dx, float dt);

void acoustic_1D_4E2T(float * P_pas, float * P_pre,float * P_fut, float * Vp, int Nx, float dx, float dt);

void acoustic_1D_8E2T(float * P_pas, float * P_pre,float * P_fut, float * Vp, int Nx, float dx, float dt);

void acoustic_2D_2E2T(int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int z_src,
                      int x_src, int nxx,int nzz, float dx, float dz, float dt);

void acoustic_2D_4E2T(int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int z_src,
                      int x_src, int nxx,int nzz, float dx, float dz, float dt);

void acoustic_2D_8E2T(int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int z_src,
                      int x_src, int nxx,int nzz, float dx, float dz, float dt);

void acoustic_3D_2E2T(int shot_pointer, int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int *z_src,
                      int *y_src, int *x_src, int n_shots, int nxx, int nyy, int nzz, 
                      float dx, float dy, float dz, float dt);

void acoustic_3D_4E2T(int shot_pointer, int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int *z_src,
                      int *y_src, int *x_src, int n_shots, int nxx, int nyy, int nzz, 
                      float dx, float dy, float dz, float dt);

void acoustic_3D_8E2T(int shot_pointer, int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int *z_src,
                      int *y_src, int *x_src, int n_shots, int nxx, int nyy, int nzz, 
                      float dx, float dy, float dz, float dt);

void elastic_isotropic_2D_wave_8E2T_tension_stencil(float *U, float *V, float *Txx, float *Tzz, float *Txz,
    float *rho, float *M, float *L, int nxx, int nzz, float dt, float ds);

void elastic_isotropic_2D_wave_8E2T_velocity_stencil(float *U, float *V, float *Txx, float *Tzz, float *Txz,
    float *rho, float *M, float *L, int nxx, int nzz, float dt, float ds);

# endif