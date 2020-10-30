# include <time.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include "acoustic1D.h"

int main(int argc, char **argv) 
{
    printf("Scalar 1D wave propagation\n\n");

    float totalTime;  
    time_t ti, tf;    
    ti = time(NULL);  

    int nz, nt, nsrc, nabc;
    float dz, dt;
    
    readParameters(&nz,&nt,&dz,&dt,&nabc,&nsrc,argv[1]);    
    int nzz = nz + 2*nabc;

    float * pas = (float *) malloc(nzz*sizeof(float));   
    float * pre = (float *) malloc(nzz*sizeof(float));   
    float * fut = (float *) malloc(nzz*sizeof(float));   
    float * vel = (float *) malloc(nzz*sizeof(float));   
    float * damp = (float *) malloc(nzz*sizeof(float));
    float * source = (float *) malloc(nsrc*sizeof(float));

    int nSnap = 1000;
    float * snap = (float *) malloc(nz*nSnap*sizeof(float));
    float * seis = (float *) malloc(nt*sizeof(float));

    importVector(vel,nzz,argv[2]);
    importVector(damp,nzz,argv[3]);
    importVector(source,nsrc,argv[4]);

    setToZero(pas,pre,fut,nzz);
    for(int kk = 0; kk < nt; kk++) 
    {
        if(kk % (nt / 10) == 0) printf("Propagation time = %.2f segundos\n",dt*kk);
        acoustic1D_8E2T(pas,pre,fut,vel,damp,source,kk,nsrc,nabc,nzz,dz,dt);    
        getSnapshots1D(snap,fut,nzz,nabc,kk,nt,nSnap);
        getSeismogram1D(seis,fut,nabc,kk);
        updateWavefield(pas,pre,fut,nzz);
    }

    exportVector(snap,nz*nSnap,"results/snapshots.bin");
    exportVector(seis,nt,"results/seismogram.bin");    

    tf = time(NULL);                
    totalTime = difftime(tf, ti);   
    printf("\nExecution time: \033[31m%.0fs\n\033[m", totalTime);

    return 0;
}
