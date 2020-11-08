# include <time.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include "auxiliaries/functions.h"

int main(int argc, char **argv) 
{
    float total_time;
    time_t t_0, t_f;
    t_0 = time(NULL);

    int nx,nz,nt,nsrc;
    int absLayer,spread,nShots; 
    float parB,dx,dz,dt; 

    readParameters(&nx,&nz,&nt,&dx,&dz,&dt,&absLayer,&spread,&nShots,&nsrc,argv[1]);

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
    float * snapshots  = (float *) malloc(nx*nz*sizeof(float));

    int *xsrc = (int *) malloc(nShots*sizeof(int));         
    int *zsrc = (int *) malloc(nShots*sizeof(int));         
    int *xrec = (int *) malloc(nShots*spread*sizeof(int));      
    int *zrec = (int *) malloc(nShots*spread*sizeof(int));   

    importFloatVector(vp,nxx*nzz,argv[2]);
    importFloatVector(damp,nxx*nzz,argv[3]);

    importIntegerVector(xsrc,nShots,argv[4]);
    importIntegerVector(zsrc,nShots,argv[5]);
    importIntegerVector(xrec,nShots*spread,argv[6]);
    importIntegerVector(zrec,nShots*spread,argv[7]);

    importFloatVector(source,nsrc,argv[8]);

    // vels = getVelocities(nxx,nzz,vp);

    // for(int shot = 0; shot < 1; ++shot) /* Shots loop */ 
    // {        
    //     setWaveField(P_pas,P_pre,P_fut,nxx*nzz);

    //     sprintf(snapsFile,"data/snapshots/snaps_shot_%i.bin",shot+1);
    //     sprintf(seismFile,"data/seismograms/shot_%i.bin",shot+1);   
                
    //     FILE * snap = fopen(snapsFile,"ab");
    //     for(int time = 0; time < nt; ++time) /* Time loop */
    //     {
    //         propagationProgress(time,shot,xsrc,nShots,xrec,spread,dx,dz,nt,vels,dt,nxx,nzz,absLayer);
    //         FDM_8E2T_acoustic2D(shot,time,vp,P_pre,P_pas,P_fut,damp,source,nsrc,zsrc,xsrc,nxx,nzz,dx,dz,dt);
    //         getSeismograms(seismogram,P_pre,xrec,zrec,spread,nxx,shot,time);            
    //         getSnapshots(snap,snapshots,P_pre,vp,nxx,nzz,absLayer,time,nt,50,0);
    //         waveFieldUpdate(P_pas,P_pre,P_fut,nxx*nzz);
    //     }
    //     exportVector(seismogram,nt*spread,seismFile);
    //     fclose(snap);
    // }

    t_f = time(NULL);                
    total_time = difftime(t_f, t_0); 
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);

    return 0;
}
