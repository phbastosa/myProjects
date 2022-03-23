# include <omp.h>
# include <time.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include "auxiliaries/functions.h"

int main(int argc,char**argv) 
{
    float total_time;
    time_t t_0, t_f;
    t_0 = time(NULL);
    
    int nx,ny,nz,nt,nabc,wbh;
    int nrecx,nrecy,nsrc;
    float dx,dy,dz,dt;

    readParameters(&nx,&ny,&nz,&nt,&dx,&dy,&dz,&dt,&nabc,&nrecx,&nrecy,&nsrc,&wbh,argv[1]);

    int nxx = nx + 2*nabc; 
    int nyy = ny + 2*nabc; 
    int nzz = nz + 2*nabc;
    int nPoints = nxx*nyy*nzz;

    float *vp   = (float *) malloc(nPoints*sizeof(float));
    float *vs   = (float *) malloc(nPoints*sizeof(float));
    float *rho  = (float *) malloc(nPoints*sizeof(float));
    float *M    = (float *) malloc(nPoints*sizeof(float));
    float *L    = (float *) malloc(nPoints*sizeof(float));
    float *damp = (float *) malloc(nPoints*sizeof(float));

    float *Vx  = (float *) malloc(nPoints*sizeof(float));
    float *Vy  = (float *) malloc(nPoints*sizeof(float));
    float *Vz  = (float *) malloc(nPoints*sizeof(float));
    float *Txx = (float *) malloc(nPoints*sizeof(float));
    float *Tyy = (float *) malloc(nPoints*sizeof(float));
    float *Tzz = (float *) malloc(nPoints*sizeof(float));
    float *Txz = (float *) malloc(nPoints*sizeof(float));
    float *Tyz = (float *) malloc(nPoints*sizeof(float));
    float *Txy = (float *) malloc(nPoints*sizeof(float));

    float *P   = (float *) malloc(nPoints*sizeof(float));
    float *Sv  = (float *) malloc(nPoints*sizeof(float));
    float *Shx = (float *) malloc(nPoints*sizeof(float));
    float *Shy = (float *) malloc(nPoints*sizeof(float)); 
    float *Ps  = (float *) malloc(nPoints*sizeof(float));

    importFloatVector(vp,nPoints,argv[2]);
    importFloatVector(vs,nPoints,argv[3]);
    importFloatVector(rho,nPoints,argv[4]);
    importFloatVector(damp,nPoints,argv[5]);

    for (int index = 0; index < nPoints; index++)
    {
        M[index] = rho[index]*pow(vs[index],2.0f);
        L[index] = rho[index]*pow(vp[index],2.0f) - 2.0f*M[index];
    }

    int nrecs = nrecx*nrecy;
    int nshot = nrecx*nrecy;
   
    int *xsrc = (int *) malloc(nshot*sizeof(int)); 
    int *ysrc = (int *) malloc(nshot*sizeof(int)); 
    int *zsrc = (int *) malloc(nshot*sizeof(int));
    int *xrec = (int *) malloc(nrecs*sizeof(int)); 
    int *yrec = (int *) malloc(nrecs*sizeof(int)); 
    int *zrec = (int *) malloc(nrecs*sizeof(int)); 

    for (int i = 0; i < nshot; i++) zsrc[i] = nabc+10;
    for (int i = 0; i < nrecs; i++) zrec[i] = nabc+10;

    float *source = (float *) malloc(nsrc*sizeof(float));

    importFloatVector(source,nsrc,argv[6]);
    importIntegerVector(xsrc,nshot,argv[7]);
    importIntegerVector(ysrc,nshot,argv[8]);
    importIntegerVector(xrec,nrecs,argv[9]);
    importIntegerVector(yrec,nrecs,argv[10]);
    
    float * seismPs = (float *) malloc(nt*nrecs*sizeof(float));
    float * seismVx = (float *) malloc(nt*nrecs*sizeof(float));
    float * seismVy = (float *) malloc(nt*nrecs*sizeof(float));
    float * seismVz = (float *) malloc(nt*nrecs*sizeof(float));
    
    float * seismP   = (float *) malloc(nt*nrecs*sizeof(float));
    float * seismSv  = (float *) malloc(nt*nrecs*sizeof(float));
    float * seismShx = (float *) malloc(nt*nrecs*sizeof(float));
    float * seismShy = (float *) malloc(nt*nrecs*sizeof(float));

    float * zeroOffsetTrace = (float *) malloc(nt*sizeof(float));

    for (int shotPointer = floor(nshot/2); shotPointer < floor(nshot/2 + 1); shotPointer++)
    {
        setWaveField(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,nxx*nyy*nzz);

        #pragma acc enter data copyin(Vx[0:nPoints],Vy[0:nPoints],Vz[0:nPoints],Txx[0:nPoints],Tyy[0:nPoints],Tzz[0:nPoints])
        #pragma acc enter data copyin(Txy[0:nPoints],Txz[0:nPoints],Tyz[0:nPoints],rho[0:nPoints],M[0:nPoints],L[0:nPoints])
        #pragma acc enter data copyin(damp[0:nPoints],source[0:nsrc],xsrc[0:nshot],ysrc[0:nshot],zsrc[0:nshot],xrec[0:nrecs],Ps[0:nPoints])
        #pragma acc enter data copyin(yrec[0:nrecs],zrec[0:nrecs],seismVx[0:nt*nrecs],seismVy[0:nt*nrecs],seismVz[0:nt*nrecs],seismPs[0:nt*nrecs])
        {
            for (int timePointer = 0; timePointer < nt; timePointer++) 
            {
                if (timePointer % (nt/10) == 0) printf("Propagation time = %.3f\n",timePointer*dt);

                FDM8E2T_elasticIsotropic3D(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,rho,M,L,source,xsrc,ysrc,zsrc,nsrc,timePointer,shotPointer,nshot,nxx,nyy,nzz,dx,dy,dz,dt);        
                
                cerjanElasticAbsorbingCondition3D(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,damp,nxx*nyy*nzz);
                    
                getPressureWaveField(Txx,Tyy,Tzz,Ps,nxx*nyy*nzz);
                
                getSeismogram(seismVx,Vx,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
                getSeismogram(seismVy,Vy,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
                getSeismogram(seismVz,Vz,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
                getSeismogram(seismPs,Ps,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
            }
        }
        #pragma acc exit data delete(Vx[0:nPoints],Vy[0:nPoints],Vz[0:nPoints],Txx[0:nPoints],Tyy[0:nPoints],Tzz[0:nPoints],zrec[0:nrecs])
        #pragma acc exit data delete(Txy[0:nPoints],Txz[0:nPoints],Tyz[0:nPoints],rho[0:nPoints],M[0:nPoints],L[0:nPoints],yrec[0:nrecs])
        #pragma acc exit data delete(damp[0:nPoints],source[0:nsrc],xsrc[0:nshot],ysrc[0:nshot],zsrc[0:nshot],xrec[0:nrecs],Ps[0:nPoints])
        #pragma acc exit data copyout(seismVx[0:nt*nrecs],seismVy[0:nt*nrecs],seismVz[0:nt*nrecs],seismPs[0:nt*nrecs])

        exportVector(seismVx,nt*nrecs,"results/seismVx.bin");
        exportVector(seismVy,nt*nrecs,"results/seismVy.bin");
        exportVector(seismVz,nt*nrecs,"results/seismVz.bin");
        exportVector(seismPs,nt*nrecs,"results/seismPs.bin");
    }

    t_f = time(NULL);
    total_time = difftime(t_f, t_0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);    
    return 0;
}
