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

    float *vp = (float *) malloc(nx*nz*sizeof(float));  
    float *b  = (float *) malloc(nx*nz*sizeof(float));
    float *K  = (float *) malloc(nx*nz*sizeof(float)); 
    float *Vx = (float *) malloc(nx*nz*sizeof(float));
    float *Vz = (float *) malloc(nx*nz*sizeof(float));
    float *P  = (float *) malloc(nx*nz*sizeof(float));

    float *source = (float *) malloc(nsrc*sizeof(float));

    for(int index = 0; index < nx*nz; index++)
    {
        vp[index] = 1500.0f;
        b[index] = 1.0f / 1000.0f;
        K[index] = (1.0f / b[index]) * powf(vp[index],2.0f);
    }

    importFloatVector(source,nsrc,argv[2]);
    
    # pragma acc enter data copyin(Vx[0:nx*nz],Vz[0:nx*nz],P[0:nx*nz],K[0:nx*nz],b[0:nx*nz])
    # pragma acc enter data copyin(source[0:nsrc],seismogram[0:nx*nt])
    for(int timePointer = 0; timePointer < nt; timePointer++) 
    {    
        if(timePointer % 200 == 0) printf("Propagation time = %0.5f seconds\n", (timePointer+200)*dt);
                  
        FDM8E2T_acousticVec2D(timePointer,Vx,Vz,P,K,b,source,nsrc,xsrc,zsrc,nx,nz,dx,dz,dt);    

        getAcousticPressureSeismogram(seismogram,P,nt,nx,nz,timePointer,zrec);
    }
    # pragma acc exit data delete(Vx[0:nx*nz],Vz[0:nx*nz],P[0:nx*nz],K[0:nx*nz],b[0:nx*nz],source[0:nsrc])
    # pragma acc exit data copyout(seismogram[0:nx*nt])

    exportVector(seismogram,nx*nt,"results/seismograms/pressureVectorialAcoustic2D.bin");

    tf = time(NULL);
    totalTime = difftime(tf, t0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", totalTime);

    return 0;
}
