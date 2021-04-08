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
    int nabc,nrecs,nshot; 
    float dx,dz,dt; 

    readParameters(&nx,&nz,&nt,&dx,&dz,&dt,&nabc,&nrecs,&nshot,&nsrc,&wbh,argv[1]);

    int nxx = nx + 2*nabc;
    int nzz = nz + 2*nabc;

    float *vp   = (float *) malloc(nxx*nzz*sizeof(float));  
    float *rho  = (float *) malloc(nxx*nzz*sizeof(float));
    float *K    = (float *) malloc(nxx*nzz*sizeof(float));
    float *b    = (float *) malloc(nxx*nzz*sizeof(float));
    float *Vx   = (float *) malloc(nxx*nzz*sizeof(float));
    float *Vz   = (float *) malloc(nxx*nzz*sizeof(float));
    float *P    = (float *) malloc(nxx*nzz*sizeof(float));
    float *damp = (float *) malloc(nxx*nzz*sizeof(float));
    
    float *source = (float *) malloc(nsrc*sizeof(float));
    
    int *seaTop = (int *) malloc(nxx*sizeof(int));
    int *seaBot = (int *) malloc(nxx*sizeof(int));
    
    int *xrec = (int *) malloc(nrecs*sizeof(int));
    int *xsrc = (int *) malloc(nshot*sizeof(int));
    
    float *seismogram = (float *) malloc(nrecs*nt*sizeof(float));
    float *data = (float *) malloc(nshot*nrecs*nt*sizeof(float));

    importIntegerVector(xrec,nrecs,argv[2]);
    importIntegerVector(xsrc,nshot,argv[3]);
    
    ajustCoordinates(xrec,xsrc,seaTop,seaBot,wbh,nabc,nx,nrecs,nshot);

    importFloatVector(source,nsrc,argv[4]);

    importFloatVector(vp,nxx*nzz,argv[5]);
    importFloatVector(rho,nxx*nzz,argv[6]);
    importFloatVector(damp,nxx*nzz,argv[7]);

    for(int index = 0; index < nxx*nzz; index++)
    {
        b[index] = 1.0f / rho[index];
        K[index] = rho[index] * powf(vp[index],2.0f);
    }

    float *vels = getVelocities(nxx,nzz,vp);

    for(int shotPointer = 0; shotPointer < nshot; shotPointer++)
    {
        wavefield_set(Vx,Vz,P,nxx,nzz); 

        # pragma acc enter data copyin(rho[0:nxx*nzz],b[0:nxx*nzz],Vx[0:nxx*nzz],Vz[0:nxx*nzz],P[0:nxx*nzz],source[0:nsrc])
        # pragma acc enter data copyin(damp[0:nxx*nzz],seismogram[0:nt*nrecs],seaTop[0:nxx],seaBot[0:nxx],xrec[0:nrecs],xsrc[0:nshot])
        
        for(int timePointer = 0; timePointer < nt; timePointer++) 
        {    
            modelingStatus(shotPointer,timePointer,xsrc,nshot,xrec,nrecs,dx,dz,nt,vels,dt,nxx,nzz,nabc);
            
            FDM8E2T_acousticVec2D(shotPointer,timePointer,Vx,Vz,P,K,b,source,nsrc,xsrc,seaTop,nxx,nzz,dx,dz,dt,nshot);    

            cerjanAbsorbingBoundaryCondition(Vx,Vz,P,damp,nxx*nzz);          

            getSeismogram(seismogram,P,xrec,seaBot,nrecs,nxx,nzz,nt,shotPointer,timePointer);
        }

        # pragma acc exit data delete(rho[0:nxx*nzz],b[0:nxx*nzz],Vx[0:nxx*nzz],Vz[0:nxx*nzz],P[0:nxx*nzz],source[0:nsrc])
        # pragma acc exit data delete(damp[0:nxx*nzz],seaTop[0:nxx],seaBot[0:nxx],xrec[0:nrecs],xsrc[0:nshot])
        # pragma acc exit data copyout(seismogram[0:nt*nrecs])

        joiningSeismograms(seismogram,data,nrecs,nt,shotPointer);
    }

    exportVector(data,nshot*nrecs*nt,"data/acousticVec_marmousi2_dh5_dataset.bin");

    t_f = time(NULL);
    total_time = difftime(t_f, t_0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);

    return 0;
} 
