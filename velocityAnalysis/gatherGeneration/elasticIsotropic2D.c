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
    int nabc,nshot,nrecs; 
    float dx,dz,dt; 

    readParameters(&nx,&nz,&nt,&dx,&dz,&dt,&nabc,&nrecs,&nshot,&nsrc,&wbh,argv[1]);
    
    int nxx = nx + 2*nabc;
    int nzz = nz + 2*nabc;
 
    float *vp   = (float *) malloc(nxx*nzz*sizeof(float));  
    float *vs   = (float *) malloc(nxx*nzz*sizeof(float));
    float *rho  = (float *) malloc(nxx*nzz*sizeof(float));
    float *M    = (float *) malloc(nxx*nzz*sizeof(float));
    float *L    = (float *) malloc(nxx*nzz*sizeof(float));
    float *Vx   = (float *) malloc(nxx*nzz*sizeof(float));
    float *Vz   = (float *) malloc(nxx*nzz*sizeof(float));
    float *Txx  = (float *) malloc(nxx*nzz*sizeof(float));
    float *Tzz  = (float *) malloc(nxx*nzz*sizeof(float));
    float *Txz  = (float *) malloc(nxx*nzz*sizeof(float));
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
    importFloatVector(vs,nxx*nzz,argv[6]);
    importFloatVector(rho,nxx*nzz,argv[7]);
    importFloatVector(damp,nxx*nzz,argv[8]);

    for(int index = 0; index < nxx*nzz; index++) 
    {
        M[index] = rho[index]*pow(vs[index],2.0f);
        L[index] = rho[index]*pow(vp[index],2.0f) - 2.0f*M[index];
    }

    float * vels = getVelocities(nxx,nzz,vp);

    for(int shotPointer = 0; shotPointer < nshot; shotPointer++)
    {
        wavefield_set(Vx,Vz,Txx,Tzz,Txz,nxx,nzz);

        #pragma acc enter data copyin(rho[0:nxx*nzz],M[0:nxx*nzz],L[0:nxx*nzz],Vx[0:nxx*nzz],Vz[0:nxx*nzz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz],Txz[0:nxx*nzz])
        #pragma acc enter data copyin(damp[0:nxx*nzz],seismogram[0:nt*nrecs],seaTop[0:nxx],seaBot[0:nxx],xrec[0:nrecs],xsrc[0:nshot],source[0:nsrc])
        
        for(int timePointer = 0; timePointer < nt; timePointer++) 
        {    
            modelingStatus(shotPointer,timePointer,xsrc,nshot,xrec,nrecs,dx,dz,nt,vels,dt,nxx,nzz,nabc);
            
            FDM8E2T_stressStencil_elasticIsotropic2D(Vx,Vz,Txx,Tzz,Txz,rho,M,L,nxx,nzz,dt,dx,dz,seaTop,xsrc,source,nsrc,nshot,timePointer,shotPointer);

            cerjanAbsorbingBoundaryCondition(Vx,Vz,Txx,Tzz,Txz,damp,nxx,nzz);          

            getSeismogram(seismogram,Txx,Tzz,xrec,seaBot,nrecs,nxx,nzz,nt,shotPointer,timePointer);
        }
        
        #pragma acc exit data delete(rho[0:nxx*nzz],M[0:nxx*nzz],L[0:nxx*nzz],Vx[0:nxx*nzz],Vz[0:nxx*nzz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz],Txz[0:nxx*nzz])
        #pragma acc exit data delete(damp[0:nxx*nzz],seaTop[0:nxx],seaBot[0:nxx],xrec[0:nrecs],xsrc[0:nshot],source[0:nsrc])
        #pragma acc exit data copyout(seismogram[0:nt*nrecs])

        joiningSeismograms(seismogram,data,nrecs,nt,shotPointer);
    }

    exportVector(data,nshot*nrecs*nt,"data/elasticIsotropic_marmousi2_dh5_dataset.bin");

    t_f = time(NULL);
    total_time = difftime(t_f, t_0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);

    return 0;
} 
