# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "../auxCodes/funtions.h"

int main(int argc, char **argv) 
{
    printf("\n\nBasic scalar acoustic modeling\n");

    float totalTime;
    time_t t0, tf;
    t0 = time(NULL);

    int nx, nz, nt, nsrc;
    int xsrc, zsrc, zrec;
    float dx, dz, dt;

    readParameters(&nx,&nz,&nt,&dx,&dz,&dt,&nsrc,&xsrc,&zsrc,&zrec,argv[1]);

    float *seismogram = (float *) malloc(nx*nt*sizeof(float));

    float *vp  = (float *) malloc(nx*nz*sizeof(float));  
    float *Pb  = (float *) malloc(nx*nz*sizeof(float));
    float *Pc  = (float *) malloc(nx*nz*sizeof(float));
    float *Pf  = (float *) malloc(nx*nz*sizeof(float));

    float *source = (float *) malloc(nsrc*sizeof(float));

    for(int index = 0; index < nx*nz; index++) vp[index] = 1500.0f;

    importFloatVector(source,nsrc,argv[2]);

    for(int timePointer = 0; timePointer < nt; timePointer++) 
    {    
        if(timePointer % 200 == 0) printf("Propagation time = %0.5f seconds\n", (timePointer+200)*dt);
                  
        FDM8E2T_acoustic2D(timePointer,vp,Pc,Pb,Pf,source,nsrc,xsrc,zsrc,nx,nz,dx,dz,dt);
        
        getAcousticPressureSeismogram(seismogram,Pf,nt,nx,nz,timePointer,zrec);
        
        waveFieldUpdate(Pb,Pc,Pf,nx*nz);
    }

    exportVector(seismogram,nx*nt,"results/seismograms/pressureScalarAcoustic2D.bin");

    tf = time(NULL);
    totalTime = difftime(tf, t0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", totalTime);

    return 0;
}
