# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include "../../modeling_library/header_files/io_functions.h"
# include "../../modeling_library/header_files/wavelet_functions.h"
# include "../../modeling_library/header_files/cerjan_functions.h"
# include "../../modeling_library/header_files/wave_equation_FDM_functions.h"

void wavefield_set(float *V, float *U, float *Txx, float *Tzz, float *Txz, int nxx, int nzz) {
    mem_set(V,nxx*nzz);
    mem_set(U,nxx*nzz);
    mem_set(Txx,nxx*nzz);
    mem_set(Tzz,nxx*nzz);
    mem_set(Txz,nxx*nzz);
}

int main(int argc, char **argv) {

    printf("\nElastic wave modelling in isotropic media using velocity stencil\n\n");

    char shot_number[200];
    int borda = 200;
    int nsrc = 600;
    int nt = 2000;
    int nxx = 500; 
    int nzz = 500;
    float pi = 4*atan(1); 
    int index, ii, jj, kk, pp;

    float total_time;
    time_t t_0, t_f;
    t_0 = time(NULL);

    float *seismogram = (float *) malloc(nxx*nt*sizeof(float));
    float *snapshot = (float *) malloc(nxx*nzz*sizeof(float));  
    float *vp  = (float *) malloc(nxx*nzz*sizeof(float));  
    float *vs  = (float *) malloc(nxx*nzz*sizeof(float));
    float *rho = (float *) malloc(nxx*nzz*sizeof(float));
    float *M   = (float *) malloc(nxx*nzz*sizeof(float));
    float *L   = (float *) malloc(nxx*nzz*sizeof(float));
    float *U   = (float *) malloc(nxx*nzz*sizeof(float));
    float *V   = (float *) malloc(nxx*nzz*sizeof(float));
    float *Txx = (float *) malloc(nxx*nzz*sizeof(float));
    float *Tzz = (float *) malloc(nxx*nzz*sizeof(float));
    float *Txz = (float *) malloc(nxx*nzz*sizeof(float));
    
    float *ricker = (float *) malloc(nsrc*sizeof(float));
    float *fator  = (float *) malloc(borda*sizeof(float));

    float *up_left = (float *) malloc(borda*borda*sizeof(float));
    float *up_right = (float *) malloc(borda*borda*sizeof(float));
    float *down_left = (float *) malloc(borda*borda*sizeof(float));
    float *down_right = (float *) malloc(borda*borda*sizeof(float));
    
    // vp = import_float32("../models/vp.bin",nxx*nzz);
    // vs = import_float32("../models/vs.bin",nxx*nzz);
    // rho = import_float32("../models/rho.bin",nxx*nzz);

    for(index = 0; index < nxx*nzz; index++) {
        vp[index] = 2000;
        vs[index] = vp[index]/sqrt(3);
        rho[index] = 310*pow(vp[index],0.25);
        M[index] = rho[index]*pow(vs[index],2.0);
        L[index] = rho[index]*pow(vp[index],2.0) - 2*M[index];
    }

    float ds = 2.0;          /* Spatial discretization parameter */  
    float dt = 0.0003;       /* Temporal discretization parameter */ 
    float f_corte = 100.0;   /* Source cutoff frequency */

    int stop = describe_2D_model_stability(vp,vs,nxx,nzz,f_corte,nt,dt,ds);
    if (stop == 1) return 0;

    ricker = build_first_gaussian_derivative(nsrc,f_corte,dt);
    fator = factor_attenuation(0.0008,borda);
    Cerjan_2D_corners(up_left,up_right,down_left,down_right,fator,borda);

    wavefield_set(V,U,Txx,Tzz,Txz,nxx,nzz);

    for(kk = 0; kk < nt; kk++) {    

        if(kk <= nsrc) {
            Tzz[(250)*nxx + 250] += ricker[kk];
            Txx[(250)*nxx + 250] += ricker[kk];                    
        }        
                                    
        elastic_isotropic_2D_wave_8E2T_velocity_stencil(U,V,Txx,Tzz,Txz,rho,M,L,nxx,nzz,dt,ds);

        Cerjan_2D_elastic_attenuation(U,V,Txx,Tzz,Txz,nxx,nzz,fator,up_left,up_right,down_left,down_right,borda);            

        if(kk % 200 == 0) printf("Propagation time = %0.5f seconds\n", kk*dt);

        exporting_2D_snapshots(kk,snapshot,V,vp,"../results/snapshots_velocity_stencil.bin",nsrc,nxx,nzz);

        for(index = 0; index < nxx; index++) 
        {
            seismogram[kk*nxx + index] = V[200*nxx + index];
        }
    }

    exporting_pointer_seismogram("../results/Vz_velocity_stencil.bin",nxx,nt,seismogram);

    t_f = time(NULL);
    total_time = difftime(t_f, t_0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);

    return 0;
} 
