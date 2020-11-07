# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include "auxiliaries/functions.h"

int main(int argc, char **argv) 
{
    float total_time;
    time_t t_0, t_f;
    t_0 = time(NULL);

    int nx,nz,nt,nsrc;
    int absLayer,spread,nShots; 
    float parB,dx,dz,dt; 

    readParameters(&nx,&nz,&nt,&dx,&dz,&dt,&absLayer,&parB,&spread,&nShots,&nsrc,argv[1]);

    int nxx = nx + 2*absLayer;
    int nzz = nz + 2*absLayer;

    float * vels   = (float *) malloc(2*sizeof(float));
    float * damp   = (float *) malloc(nxx*nzz*sizeof(float));     
    float * vp     = (float *) malloc(nxx*nzz*sizeof(float));  
    float * P_pas  = (float *) malloc(nxx*nzz*sizeof(float));   
    float * P_pre  = (float *) malloc(nxx*nzz*sizeof(float));  
    float * P_fut  = (float *) malloc(nxx*nzz*sizeof(float));  
    float * source = (float *) malloc(nsrc*sizeof(float));    

    char seismFile[100], snapsFile[100];
    float * seismogram = (float *) malloc(nt*spread*sizeof(float));
    float * snapshots  = (float *) malloc(nxx*nzz*sizeof(float));

    int *xsrc = (int *) malloc(nShots*sizeof(int));         
    int *zsrc = (int *) malloc(nShots*sizeof(int));         
    int *xrec = (int *) malloc(nShots*spread*sizeof(int));      
    int *zrec = (int *) malloc(nShots*spread*sizeof(int));   

    importVector2D(vp,nxx,nzz,argv[4]);
    importVector2D(damp,nxx,nzz,argv[6]);
    vels = getVelocities(nxx,nzz,vp);

    readIntegerTextParameter(xsrc,argv[3]);
    readIntegerTextParameter(zsrc,argv[3]);
    readIntegerTextParameter(xrec,argv[2]);
    readIntegerTextParameter(zrec,argv[3]);

    readFloatTextParameter(source,argv[5]);

    FILE * seis = fopen("asdfa","ab");
    for(int shot = 0; shot < 1; ++shot) /* Shots loop */ 
    {        
        sprintf(snapsFile,"data/snapshots/",shot+1);   
        FILE * snap = fopen(snapsFile,"ab");
        
        setWaveField(P_pas,P_pre,P_fut,nxx*nzz);

        for(int time = 0; time < nt; ++time) /* Time loop */
        {
            propagationProgress(time,shot,xsrc,nShots,xrec,spread,dx,dz,nt,vels,dt,nxx,nzz,absLayer);
            FDM_8E2T_acoustic2D(shot,time,vp,P_pre,P_pas,P_fut,damp,source,nsrc,zsrc,xsrc,nxx,nzz,dx,dz,dt);
            getSeismograms(seism,P_pre,xrec,zrec,spread,nxx,shot,time);            
            
            waveFieldUpdate(P_pas,P_pre,P_fut,nxx*nzz);
        }
    }

    t_f = time(NULL);                /* Initializing the final time */
    total_time = difftime(t_f, t_0); /* Calculating the difference between the initial and final times */
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);

    return 0;
}
