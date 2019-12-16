# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include "../modeling_library/header_files/io_functions.h"
# include "../modeling_library/header_files/wavelet_functions.h"
# include "../modeling_library/header_files/cerjan_functions.h"
# include "../modeling_library/header_files/wave_equation_FDM_functions.h"

int main(int argc, char **argv) 
{
    /* Variables to count the execution time */
    float total_time;  /* Total time */
    time_t t_0, t_f;   /* Initial and final time */
    t_0 = time(NULL);  /* Starting the time */

    int i,j,k;         /* Counter variables */
    int Nx,Nz,Nt;      /* Length of the space and time in modeling */
    int sx, sz;       /* Shot position in the modeling */   
    int amost_f;       /* Total samples of the source */

    float h,dt;        /* Discretization parameters */
    float fcut;        /* Cutoff frequency of the source */                    
    
    FILE *ptr1, *ptr2; /* Input Output pointers */     

    Nx = 500;          /* Length of the model */
    Nz = 500;          /* Length of the model */
    h = 5.0;           /* Spatial discretifation parameter */
    Nt = 3000;         /* Total samples of the time */
    dt = 0.0002;       /* Time discretization parameter */
    sx = 250;          /* Shot position in the modeling */
    sz = 100;          /* Shot position in the modeling */
    fcut = 100;        /* Cutoff frequency of the synthetic seismic source */
    amost_f = 1000;    /* Total samples of the synthetic seismic source */
    int bondary = 100; /* Samples used in the attenuation layer */

    float * ricker = (float *) malloc(amost_f*sizeof(float)); /* Sinthetic seismic source*/
    float * P1 = (float *) malloc(Nx*Nz*sizeof(float));          /* Future wavefield */
    float * P2 = (float *) malloc(Nx*Nz*sizeof(float));          /* Present wavefield */
    float * P3 = (float *) malloc(Nx*Nz*sizeof(float));          /* Past wavefield */
    float * Vp = (float *) malloc(Nx*Nz*sizeof(float));          /* Compressional wave velocities */ 
    float * factor = (float *) malloc(bondary * sizeof(float)); /* Factor to absorb the amplitude of the wave */
    float (*Seism)[Nx] = malloc(sizeof(float[Nt][Nx]));       /* Seismogram, the complete wavefield */

    float * UL = (float *) malloc(bondary * bondary * sizeof(float));
    float * UR = (float *) malloc(bondary * bondary * sizeof(float));
    float * DL = (float *) malloc(bondary * bondary * sizeof(float));
    float * DR = (float *) malloc(bondary * bondary * sizeof(float));

    factor = factor_attenuation(0.004,bondary); 
    Cerjan_2D_corners(UL,UR,DL,DR,factor,bondary);

    ricker = build_ricker(amost_f,fcut,dt);

    mem_set(P1,Nx*Nz);      /* Zeroing the future wavefield */
    mem_set(P2,Nx*Nz);      /* Zeroing the presente wavefield */
    mem_set(P3,Nx*Nz);      /* Zeroing the past wavefield */
    mem_set(Vp,Nx*Nz);      /* Zeroing the velocities */

    for(i = 0; i < Nx*Nz; i++) {
        Vp[i] = 2000;    /* Contant compressional wave velocity */        
    }

    for(k = 0; k < Nt; k++) {  /* Time loop */

        if(k % 200 == 0) printf("Time = %.2f seconds\n",dt*k);

        acoustic_2D_8E2T(k,Vp,P2,P3,P1,ricker,amost_f,sz,sx,Nx,Nz,h,h,dt);    

        Cerjan_2D_acoustic_attenuation(P2,P1,factor,Nx,Nz,bondary,DL,DR,UL,UR);

        update_wavefield(P3,P2,P1,Nx*Nz);

        for(i = 0; i < Nx; i++) {
            Seism[k][i] = P2[sz*Nx + i]; /* Case 2D, we record the surface disturbances */
        }        
    }

    /* Exporting the seismogram */
    ptr2 = fopen("../results/seism_2D.bin","wb");
    for(i = 0; i < Nx; i++) {
        for(j = 0; j < Nt; j++) {
            fwrite(&Seism[j][i], sizeof(float), 1, ptr2);
        }
    }
    fclose(ptr2);

    t_f = time(NULL);                /* Initializing the final time */
    total_time = difftime(t_f, t_0); /* Calculating the difference between the initial and final times */
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);

    return 0;
}