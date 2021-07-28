# include <time.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include "auxCodes/functions.h"

int main(int argc, char **argv) 
{
    system("clear");
    printf("Vertical seismic profile simulation\n\n");

    float total_time;
    time_t t_0, t_f;
    t_0 = time(NULL);

    int nx,nz,nt,nsrc;
    int nabc,xsrc,zsrc; 
    float dx,dz,dt; 

    readParameters(&nx,&nz,&nt,&dx,&dz,&dt,&nabc,&xsrc,&zsrc,&nsrc,argv[1]);

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
            
    float *seismogram = (float *) malloc(nz*nt*sizeof(float));
    
    xsrc += nabc;
    zsrc += nabc;

    importFloatVector(source,nsrc,argv[2]);

    importFloatVector(vp,nxx*nzz,argv[3]);
    importFloatVector(vs,nxx*nzz,argv[4]);
    importFloatVector(rho,nxx*nzz,argv[5]);
    importFloatVector(damp,nxx*nzz,argv[6]);

    for(int index = 0; index < nxx*nzz; index++) 
    {
        M[index] = rho[index]*pow(vs[index],2.0f);
        L[index] = rho[index]*pow(vp[index],2.0f) - 2.0f*M[index];
    }

    wavefield_set(Vx,Vz,Txx,Tzz,Txz,nxx,nzz);

    #pragma acc enter data copyin(Vx[0:nxx*nzz],Vz[0:nxx*nzz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz],Txz[0:nxx*nzz],L[0:nxx*nzz])
    #pragma acc enter data copyin(damp[0:nxx*nzz],seismogram[0:nt*nz],source[0:nsrc],rho[0:nxx*nzz],M[0:nxx*nzz])
        
    for(int timePointer = 0; timePointer < nt; timePointer++) 
    {                
        if (timePointer % (nt / 10) == 0) printf("Propagação no tempo %i\n",timePointer);

        FDM8E2T_stressStencil_elasticIsotropic2D(Vx,Vz,Txx,Tzz,Txz,rho,M,L,nxx,nzz,dt,dx,dz,zsrc,xsrc,source,nsrc,timePointer);

        cerjanAbsorbingBoundaryCondition(Vx,Vz,Txx,Tzz,Txz,damp,nxx,nzz);          

        getSeismogram(seismogram,Txx,Tzz,xsrc,nz,nabc,nxx,nzz,nt,timePointer);
    }

    #pragma acc exit data delete(Vx[0:nxx*nzz],Vz[0:nxx*nzz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz],Txz[0:nxx*nzz])
    #pragma acc exit data delete(damp[0:nxx*nzz],source[0:nsrc],rho[0:nxx*nzz],M[0:nxx*nzz],L[0:nxx*nzz])        
    #pragma acc exit data copyout(seismogram[0:nt*nz])

    exportVector(seismogram,nz*nt,"data/vsp.bin");

    t_f = time(NULL);
    total_time = difftime(t_f, t_0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);

    return 0;
} 
