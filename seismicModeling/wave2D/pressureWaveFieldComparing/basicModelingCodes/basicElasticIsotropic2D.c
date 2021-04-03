# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include <string.h>
# include "../auxCodes/funtions.h"

int main(int argc, char **argv) 
{
    float total_time;
    time_t t_0, t_f;
    t_0 = time(NULL);

    int nsrc = 400;
    int nt = 2000;
    int nxx = 500; 
    int nzz = 500;

    float *seismogram = (float *) malloc(nxx*nt*sizeof(float));

    float *vp  = (float *) malloc(nxx*nzz*sizeof(float));  
    float *vs  = (float *) malloc(nxx*nzz*sizeof(float));
    float *rho = (float *) malloc(nxx*nzz*sizeof(float));
    float *M   = (float *) malloc(nxx*nzz*sizeof(float));
    float *L   = (float *) malloc(nxx*nzz*sizeof(float));
    float *Vx  = (float *) malloc(nxx*nzz*sizeof(float));
    float *Vz  = (float *) malloc(nxx*nzz*sizeof(float));
    float *Txx = (float *) malloc(nxx*nzz*sizeof(float));
    float *Tzz = (float *) malloc(nxx*nzz*sizeof(float));
    float *Txz = (float *) malloc(nxx*nzz*sizeof(float));
    
    float *source = (float *) malloc(nsrc*sizeof(float));
    
    importFloatVector(source,nsrc,"staggeredRicker.bin");
    
    for(int index = 0; index < nxx*nzz; index++) 
    {
        vp[index] = 2000;
        vs[index] = vp[index]/sqrt(3);
        rho[index] = 310*pow(vp[index],0.25);
        M[index] = rho[index]*pow(vs[index],2.0);
        L[index] = rho[index]*pow(vp[index],2.0) - 2*M[index];
    }

    float dx = 5.0;           
    float dz = 5.0;
    float dt = 0.001;        

    for(int kk = 0; kk < nt; kk++) 
    {    
        if(kk < nsrc) 
        {
            Tzz[(250)*nxx + 250] += source[kk] / (dx*dz);
            Txx[(250)*nxx + 250] += source[kk] / (dx*dz);                    
        }        
                                    
        FDM8E2T_stressStencil_elasticIsotropic2D(Vx,Vz,Txx,Tzz,Txz,rho,M,L,nxx,nzz,dt,dx,dz);

        if(kk % 200 == 0) printf("Propagation time = %0.5f seconds\n", kk*dt);

        for(int index = 0; index < nxx; index++) 
        {
            seismogram[kk*nxx + index] = (Txx[200*nxx + index] + Tzz[200*nxx + index])/2.0f;
            // seismogram[kk*nxx + index] = Vz[100*nxx + index];
        }
    }

    exportVector(seismogram,nxx*nt,"seismPressure_ricker.bin");
    // exportVector(seismogram,nxx*nt,"seismPressure_intRicker.bin");
    // exportVector(seismogram,nxx*nt,"seismPressure_hankelIntRicker.bin");

    t_f = time(NULL);
    total_time = difftime(t_f, t_0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);

    return 0;
} 
