# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "../auxCodes/funtions.h"

int main(int argc, char **argv) 
{
    printf("\n\nBasic vectorial acoustic modeling\n");

    float totalTime;
    time_t t0, tf;
    t0 = time(NULL);

    int nx, nz, nt, nsrc;
    int xsrc, zsrc, zrec;
    float dx, dz, dt;

    readParameters(&nx,&nz,&nt,&dx,&dz,&dt,&nsrc,&xsrc,&zsrc,&zrec,argv[1]);

    float *seismogram = (float *) malloc(nx*nt*sizeof(float));

    float *vp  = (float *) malloc(nx*nz*sizeof(float));  
    float *rho = (float *) malloc(nx*nz*sizeof(float)); 
    float *Vx  = (float *) malloc(nx*nz*sizeof(float));
    float *Vz  = (float *) malloc(nx*nz*sizeof(float));
    float *P   = (float *) malloc(nx*nz*sizeof(float));

    float *source = (float *) malloc(nsrc*sizeof(float));

    for(int index = 0; index < nx*nz; index++)
    {
        vp[index]  = 1500.0f;
        rho[index] = 1000.0f;
    }

    importFloatVector(source,nsrc,argv[2]);








    tf = time(NULL);
    totalTime = difftime(tf, t0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", totalTime);

    return 0;
}
