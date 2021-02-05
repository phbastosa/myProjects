# include <cuda.h>
# include <time.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <cuda_runtime.h>
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
    int threads = 1000;
    
    int nrecs = nrecx*nrecy;
    int nshot = (nrecx-1)*(nrecy-1);

    /* Host arrays */

    float *vp   = (float *) malloc(nPoints*sizeof(float));
    float *vs   = (float *) malloc(nPoints*sizeof(float));
    float *rho  = (float *) malloc(nPoints*sizeof(float));
    float *M    = (float *) malloc(nPoints*sizeof(float));
    float *L    = (float *) malloc(nPoints*sizeof(float));
    float *damp = (float *) malloc(nPoints*sizeof(float));

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

    float *d_Vx;  cudaMalloc(&d_Vx,nPoints*sizeof(float));
    float *d_Vy;  cudaMalloc(&d_Vy,nPoints*sizeof(float));
    float *d_Vz;  cudaMalloc(&d_Vz,nPoints*sizeof(float));
    float *d_Txx; cudaMalloc(&d_Txx,nPoints*sizeof(float));
    float *d_Tyy; cudaMalloc(&d_Tyy,nPoints*sizeof(float));
    float *d_Tzz; cudaMalloc(&d_Tzz,nPoints*sizeof(float));
    float *d_Txz; cudaMalloc(&d_Txz,nPoints*sizeof(float));
    float *d_Tyz; cudaMalloc(&d_Tyz,nPoints*sizeof(float));
    float *d_Txy; cudaMalloc(&d_Txy,nPoints*sizeof(float));

    float *d_P;   cudaMalloc(&d_P,nPoints*sizeof(float));
    float *d_Sv;  cudaMalloc(&d_Sv,nPoints*sizeof(float));
    float *d_Shx; cudaMalloc(&d_Shx,nPoints*sizeof(float));
    float *d_Shy; cudaMalloc(&d_Shy,nPoints*sizeof(float)); 
    float *d_Ps;  cudaMalloc(&d_Ps,nPoints*sizeof(float));

    float *d_rho;  cudaMalloc(&d_rho,nPoints*sizeof(float));
    float *d_M;    cudaMalloc(&d_M,nPoints*sizeof(float)); 
    float *d_L;    cudaMalloc(&d_L,nPoints*sizeof(float)); 
    float *d_damp; cudaMalloc(&d_damp,nPoints*sizeof(float)); 

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

    /* Host inputs and calculations */

    importFloatVector(vp,nPoints,argv[2]);
    importFloatVector(vs,nPoints,argv[3]);
    importFloatVector(rho,nPoints,argv[4]);
    importFloatVector(damp,nPoints,argv[5]);

    importFloatVector(source,nsrc,argv[6]);
    importIntegerVector(xsrc,nshot,argv[7]);
    importIntegerVector(ysrc,nshot,argv[8]);
    importIntegerVector(xrec,nrecs,argv[9]);
    importIntegerVector(yrec,nrecs,argv[10]);

    for (int index = 0; index < nPoints; index++)
    {
        M[index] = rho[index]*pow(vs[index],2.0f);
        L[index] = rho[index]*pow(vp[index],2.0f) - 2.0f*M[index];
    }

    for (int i = 0; i < nshot; i++) zsrc[i] = nabc + 5;
    for (int i = 0; i < nrecs; i++) zrec[i] = nabc + wbh;
    
    /* Managing memory */
    
    cudaMemcpy(d_rho,rho,nPoints*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_M,M,nPoints*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_L,L,nPoints*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_damp,damp,nPoints*sizeof(float),cudaMemcpyHostToDevice);

    cudaMemcpy(d_source,source,nsrc*sizeof(float),cudaMemcpyHostToDevice);

    cudaMemcpy(d_xsrc,xsrc,nshot*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_ysrc,ysrc,nshot*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_zsrc,zsrc,nshot*sizeof(int),cudaMemcpyHostToDevice);
    
    cudaMemcpy(d_xrec,xrec,nrecs*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_yrec,yrec,nrecs*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_zrec,zrec,nrecs*sizeof(int),cudaMemcpyHostToDevice);

    printf("Wave propagation in elastic isotropic 3D media\n\n");    
    for (int shotPointer = floor(nshot/2 - nrecx/2); shotPointer < floor(nshot/2 - nrecx/2 + 1); shotPointer++)
    {
        setWaveField<<<nPoints/threads,threads>>>(d_Vx,d_Vy,d_Vz,d_Txx,d_Tyy,d_Tzz,d_Txz,d_Tyz,d_Txy,nPoints);

        for (int timePointer = 0; timePointer < nt; timePointer++) 
        {
            if (timePointer % (nt/10) == 0) printf("Propagation time = %.3f\n",timePointer*dt);

            computeStress<<<nPoints/threads,threads>>>(d_Vx,d_Vy,d_Vz,d_Txx,d_Tyy,d_Tzz,d_Txy,d_Txz,d_Tyz,d_M,d_L,d_source,d_xsrc,d_ysrc,d_zsrc,nsrc,timePointer,shotPointer,nxx,nyy,nzz,dx,dy,dz,dt);        
            cudaDeviceSynchronize();

            computeVelocity<<<nPoints/threads,threads>>>(d_Vx,d_Vy,d_Vz,d_Txx,d_Tyy,d_Tzz,d_Txy,d_Txz,d_Tyz,d_rho,nxx,nyy,nzz,dx,dy,dz,dt);    
            cudaDeviceSynchronize();

            cerjanElasticABC3D<<<nPoints/threads,threads>>>(d_Vx,d_Vy,d_Vz,d_Txx,d_Tyy,d_Tzz,d_Txy,d_Txz,d_Tyz,d_damp,nPoints);
            cudaDeviceSynchronize();    

            getPressureWaveField<<<nPoints/threads,threads>>>(d_Txx,d_Tyy,d_Tzz,d_Ps,nPoints);
            cudaDeviceSynchronize();    
            
            getPWaveField<<<nPoints/threads,threads>>>(d_Vx,d_Vy,d_Vz,d_P,nxx,nyy,nzz,dx,dy,dz);
            cudaDeviceSynchronize();    
          
            getSWaveField<<<nPoints/threads,threads>>>(d_Vx,d_Vy,d_Vz,d_Shx,d_Shy,d_Sv,nxx,nyy,nzz,dx,dy,dz);
            cudaDeviceSynchronize();    

            getSeismogram<<<1,nrecs>>>(d_seismVx,d_Vx,d_xrec,d_yrec,d_zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer);
            cudaDeviceSynchronize();    

    //         getSeismogram(seismVy,Vy,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
    //         getSeismogram(seismVz,Vz,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
    //         getSeismogram(seismPs,Ps,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);

    //         getSeismogram(seismP,P,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
    //         getSeismogram(seismSv,Sv,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
    //         getSeismogram(seismShx,Shx,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
    //         getSeismogram(seismShy,Shy,xrec,yrec,zrec,nrecx,nrecy,nrecs,nt,nxx,nyy,nzz,timePointer,nsrc,dt);
        }
    }
    
    cudaMemcpy(seismVx,d_seismVx,nt*nrecs*sizeof(float),cudaMemcpyDeviceToHost);

    exportVector(seismVx,nt*nrecs,(char *)"results/seismVx.bin");
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