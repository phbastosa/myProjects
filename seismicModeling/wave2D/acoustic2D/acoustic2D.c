# include <time.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include "auxCodes/functions.h"

int main(int argc, char **argv) 
{
    float total_time;
    time_t t_0, t_f;
    t_0 = time(NULL);

    int nx,nz,nt,nsrc,wbh;
    int absLayer,nrecs,nshot; 
    float dx,dz,dt; 

    readParameters(&nx,&nz,&nt,&dx,&dz,&dt,&absLayer,&nrecs,&nshot,&nsrc,&wbh,argv[1]);

    int nxx = nx + 2*absLayer;
    int nzz = nz + 2*absLayer;

    float * vels   = (float *) malloc(2*sizeof(float));
    float * damp   = (float *) malloc(nxx*nzz*sizeof(float));     
    float * vp     = (float *) malloc(nxx*nzz*sizeof(float));  
    float * P_pas  = (float *) malloc(nxx*nzz*sizeof(float));   
    float * P_pre  = (float *) malloc(nxx*nzz*sizeof(float));  
    float * P_fut  = (float *) malloc(nxx*nzz*sizeof(float));  
    float * source = (float *) malloc(nsrc*sizeof(float));    

    float * seismogram = (float *) malloc(nt*nrecs*sizeof(float));
    float * data = (float *) malloc(nshot*nt*nrecs*sizeof(float));

    int *xsrc = (int *) malloc(nshot*sizeof(int));         
    int *xrec = (int *) malloc(nrecs*sizeof(int));      
    int *topo = (int *) malloc(nxx*sizeof(int));   

    importFloatVector(vp,nxx*nzz,argv[2]);
    importFloatVector(damp,nxx*nzz,argv[3]);
    importFloatVector(source,nsrc,argv[4]);

    importIntegerVector(xsrc,nshot,argv[5]);
    importIntegerVector(xrec,nrecs,argv[6]);

    vels = getVelocities(nxx,nzz,vp);

    ajustCoordinates(xrec,xsrc,topo,absLayer,nxx,nrecs,nshot);

    for(int shotPointer = 0; shotPointer < nshot; ++shotPointer) 
    {        
        setWaveField(P_pas,P_pre,P_fut,nxx*nzz);
                
        #pragma acc enter data copyin(vp[0:nxx*nzz],P_pre[0:nxx*nzz],P_pas[0:nxx*nzz],P_fut[0:nxx*nzz],damp[0:nxx*nzz])
        #pragma acc enter data copyin(seismogram[0:nrecs*nt],xsrc[0:nshot],xrec[0:nrecs],topo[0:nxx],source[0:nsrc])

        for(int timePointer = 0; timePointer < nt; ++timePointer) 
        {
            modelingStatus(shotPointer,timePointer,xsrc,nshot,xrec,nrecs,dx,dz,nt,vels,dt,nxx,nzz,absLayer);

            FDM_8E2T_acoustic2D(shotPointer,timePointer,vp,P_pre,P_pas,P_fut,source,nsrc,topo,xsrc,nxx,nzz,dx,dz,dt,nshot);

            cerjanAbsorbingBoundaryCondition(P_pas,P_pre,P_fut,damp,nxx*nzz);

            getSeismograms(seismogram,P_fut,xrec,wbh,nrecs,nxx,nzz,nt,nrecs,shotPointer,timePointer);            

            waveFieldUpdate(P_pas,P_pre,P_fut,nxx*nzz);
        }
        
        #pragma acc exit data delete(vp[0:nxx*nzz],P_pre[0:nxx*nzz],P_pas[0:nxx*nzz],P_fut[0:nxx*nzz])
        #pragma acc exit data delete(xsrc[0:nshot],xrec[0:nrecs],topo[0:nxx],source[0:nsrc],damp[0:nxx*nzz])
        #pragma acc exit data copyout(seismogram[0:nrecs*nt])        

        joiningSeismograms(seismogram,data,nrecs,nt,shotPointer);
    }
    exportVector(data,nshot*nt*nrecs,"data/acoustic_marmousi2_dh5_dataset.bin");

    t_f = time(NULL);                
    total_time = difftime(t_f, t_0); 
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);

    return 0;
}
