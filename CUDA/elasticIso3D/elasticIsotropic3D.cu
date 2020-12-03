# include <cuda.h>
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

    int nrecs = nrecx*nrecy;
    int nshot = (nrecx-1)*(nrecy-1);

    /* Host arrays */

    float *vp   = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *vs   = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *rho  = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *M    = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *L    = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *damp = (float *) malloc(nxx*nyy*nzz*sizeof(float));

    float *Vx  = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *Vy  = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *Vz  = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *Txx = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *Tyy = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *Tzz = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *Txz = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *Tyz = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *Txy = (float *) malloc(nxx*nyy*nzz*sizeof(float));

    float *P   = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *Sv  = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *Shx = (float *) malloc(nxx*nyy*nzz*sizeof(float));
    float *Shy = (float *) malloc(nxx*nyy*nzz*sizeof(float)); 
    float *Ps  = (float *) malloc(nxx*nyy*nzz*sizeof(float));

    float *source = (float *) malloc(nsrc*sizeof(float));
    
    int *xsrc = (int *) malloc(nshot*sizeof(int)); 
    int *ysrc = (int *) malloc(nshot*sizeof(int)); 
    int *zsrc = (int *) malloc(nshot*sizeof(int));
    int *xrec = (int *) malloc(nrecs*sizeof(int)); 
    int *yrec = (int *) malloc(nrecs*sizeof(int)); 
    int *zrec = (int *) malloc(nrecs*sizeof(int)); 

    float *seismPs = (float *) malloc(nt*nrecs*sizeof(float));
    float *seismVx = (float *) malloc(nt*nrecs*sizeof(float));
    float *seismVy = (float *) malloc(nt*nrecs*sizeof(float));
    float *seismVz = (float *) malloc(nt*nrecs*sizeof(float));
    
    float *seismP   = (float *) malloc(nt*nrecs*sizeof(float));
    float *seismSv  = (float *) malloc(nt*nrecs*sizeof(float));
    float *seismShx = (float *) malloc(nt*nrecs*sizeof(float));
    float *seismShy = (float *) malloc(nt*nrecs*sizeof(float));

    /* Device arrays */

    float *d_Vx;  cudaMalloc(&d_Vx,nxx*nyy*nzz*sizeof(float));
    float *d_Vy;  cudaMalloc(&d_Vy,nxx*nyy*nzz*sizeof(float));
    float *d_Vz;  cudaMalloc(&d_Vz,nxx*nyy*nzz*sizeof(float));
    float *d_Txx; cudaMalloc(&d_Txx,nxx*nyy*nzz*sizeof(float));
    float *d_Tyy; cudaMalloc(&d_Tyy,nxx*nyy*nzz*sizeof(float));
    float *d_Tzz; cudaMalloc(&d_Tzz,nxx*nyy*nzz*sizeof(float));
    float *d_Txz; cudaMalloc(&d_Txz,nxx*nyy*nzz*sizeof(float));
    float *d_Tyz; cudaMalloc(&d_Tyz,nxx*nyy*nzz*sizeof(float));
    float *d_Txy; cudaMalloc(&d_Txy,nxx*nyy*nzz*sizeof(float));

    float *d_P;   cudaMalloc(&d_P,nxx*nyy*nzz*sizeof(float));
    float *d_Sv;  cudaMalloc(&d_Sv,nxx*nyy*nzz*sizeof(float));
    float *d_Shx; cudaMalloc(&d_Shx,nxx*nyy*nzz*sizeof(float));
    float *d_Shy; cudaMalloc(&d_Shy,nxx*nyy*nzz*sizeof(float)); 
    float *d_Ps;  cudaMalloc(&d_Ps,nxx*nyy*nzz*sizeof(float));

    float *d_rho;  cudaMalloc(&d_rho,nxx*nyy*nzz*sizeof(float));
    float *d_M;    cudaMalloc(&d_M,nxx*nyy*nzz*sizeof(float)); 
    float *d_L;    cudaMalloc(&d_L,nxx*nyy*nzz*sizeof(float)); 
    float *d_damp; cudaMalloc(&d_damp,nxx*nyy*nzz*sizeof(float)); 

    float *d_source; cudaMalloc(&d_source,nsrc*sizeof(float));
    
    int *d_xsrc; cudaMalloc(&d_xsrc,nshot*sizeof(int)); 
    int *d_ysrc; cudaMalloc(&d_ysrc,nshot*sizeof(int)); 
    int *d_zsrc; cudaMalloc(&d_zsrc,nshot*sizeof(int));
    int *d_xrec; cudaMalloc(&d_xrec,nrecs*sizeof(int)); 
    int *d_yrec; cudaMalloc(&d_yrec,nrecs*sizeof(int)); 
    int *d_zrec; cudaMalloc(&d_zrec,nrecs*sizeof(int)); 

    float *d_seismPs; cudaMalloc(&d_seismPs,nt*nrecs*sizeof(float));
    float *d_seismVx; cudaMalloc(&d_seismVx,nt*nrecs*sizeof(float));
    float *d_seismVy; cudaMalloc(&d_seismVy,nt*nrecs*sizeof(float));
    float *d_seismVz; cudaMalloc(&d_seismVz,nt*nrecs*sizeof(float));
    
    float *d_seismP;   cudaMalloc(&d_seismP,nt*nrecs*sizeof(float));
    float *d_seismSv;  cudaMalloc(&d_seismSv,nt*nrecs*sizeof(float));
    float *d_seismShx; cudaMalloc(&d_seismShx,nt*nrecs*sizeof(float));
    float *d_seismShy; cudaMalloc(&d_seismShy,nt*nrecs*sizeof(float));








    // importFloatVector(vp,nxx*nyy*nzz,argv[2]);
    // importFloatVector(vs,nxx*nyy*nzz,argv[3]);
    // importFloatVector(rho,nxx*nyy*nzz,argv[4]);
    // importFloatVector(damp,nxx*nyy*nzz,argv[5]);

    // importFloatVector(source,nsrc,argv[6]);
    // importIntegerVector(xsrc,nshot,argv[7]);
    // importIntegerVector(ysrc,nshot,argv[8]);
    // importIntegerVector(xrec,nrecs,argv[9]);
    // importIntegerVector(yrec,nrecs,argv[10]);

    // for (int index = 0; index < nxx*nyy*nzz; index++)
    // {
    //     M[index] = rho[index]*pow(vs[index],2.0f);
    //     L[index] = rho[index]*pow(vp[index],2.0f) - 2.0f*M[index];
    // }

    // for (int i = 0; i < nshot; i++) zsrc[i] = nabc + 5;
    // for (int i = 0; i < nrecs; i++) zrec[i] = nabc + wbh;
    
    /* Managing memory */
    







    // for (int shotPointer = floor(nshot/2 - nrecx/2); shotPointer < floor(nshot/2 - nrecx/2 + 1); shotPointer++)
    // {
    //     setWaveField(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,nxx*nyy*nzz);

    //     for (int timePointer = 0; timePointer < nt; timePointer++) 
    //     {
    //         if (timePointer % (nt/10) == 0) printf("Propagation time = %.3f\n",timePointer*dt);

    //         FDM8E2T_elasticIsotropic3D(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,rho,M,L,source,xsrc,ysrc,zsrc,nsrc,timePointer,shotPointer,nxx,nyy,nzz,dx,dy,dz,dt);        
    //         cerjanElasticAbsorbingCondition3D(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,damp,nxx*nyy*nzz);
                
    //         getPressureWaveField(Txx,Tyy,Tzz,Ps,nxx*nyy*nzz);
    //         getPWaveField(Vx,Vy,Vz,P,nxx,nyy,nzz,dx,dy,dz);
    //         getSWaveField(Vx,Vy,Vz,Shx,Shy,Sv,nxx,nyy,nzz,dx,dy,dz);

    //         getSeismogram(seismVx,Vx,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
    //         getSeismogram(seismVy,Vy,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
    //         getSeismogram(seismVz,Vz,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
    //         getSeismogram(seismPs,Ps,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);

    //         getSeismogram(seismP,P,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
    //         getSeismogram(seismSv,Sv,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
    //         getSeismogram(seismShx,Shx,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
    //         getSeismogram(seismShy,Shy,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
    //     }
    // }
    
    // exportVector(seismVx,nt*nrecs,(char *)"results/seismVx.bin");
    // exportVector(seismVy,nt*nrecs,(char *)"results/seismVy.bin");
    // exportVector(seismVz,nt*nrecs,(char *)"results/seismVz.bin");
    // exportVector(seismPs,nt*nrecs,(char *)"results/seismPs.bin");

    // exportVector(seismP,nt*nrecs,(char *)"results/seismP.bin");
    // exportVector(seismSv,nt*nrecs,(char *)"results/seismSv.bin");
    // exportVector(seismShx,nt*nrecs,(char *)"results/seismShx.bin");
    // exportVector(seismShy,nt*nrecs,(char *)"results/seismShy.bin");

    t_f = time(NULL);
    total_time = difftime(t_f, t_0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);    
    return 0;
} 