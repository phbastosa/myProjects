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
    
    int nx,ny,nz,nt,nabc,wbh;
    int nrecx,nrecy,nsrc;
    float dx,dy,dz,dt;

    readParameters(&nx,&ny,&nz,&nt,&dx,&dy,&dz,&dt,&nabc,&nrecx,&nrecy,&nsrc,&wbh,argv[1]);

    int nxx = nx + 2*nabc; int nyy = ny + 2*nabc; int nzz = nz + 2*nabc;
    
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
    float *Pss = (float *) malloc(nxx*nyy*nzz*sizeof(float));

    importFloatVector(vp,nxx*nyy*nzz,argv[2]);
    importFloatVector(vs,nxx*nyy*nzz,argv[3]);
    importFloatVector(rho,nxx*nyy*nzz,argv[4]);
    importFloatVector(damp,nxx*nyy*nzz,argv[5]);

    for (int index = 0; index < nxx*nyy*nzz; index++)
    {
        M[index] = rho[index]*pow(vs[index],2.0f);
        L[index] = rho[index]*pow(vp[index],2.0f) - 2.0f*M[index];
    }

    int nrecs = nrecx*nrecy;
    int nshot = (nrecx-1)*(nrecy-1);
   
    int *xsrc = (int *) malloc(nshot*sizeof(int)); 
    int *ysrc = (int *) malloc(nshot*sizeof(int)); 
    int *zsrc = (int *) malloc(nshot*sizeof(int));
    int *xrec = (int *) malloc(nrecs*sizeof(int)); 
    int *yrec = (int *) malloc(nrecs*sizeof(int)); 
    int *zrec = (int *) malloc(nrecs*sizeof(int)); 

    for (int i = 0; i < nshot; i++) zsrc[i] = nabc;
    for (int i = 0; i < nrecs; i++) zrec[i] = wbh;

    float *source = (float *) malloc(nsrc*sizeof(float));
    
    importFloatVector(source,nsrc,argv[6]);
    importIntegerVector(xsrc,nshot,argv[7]);
    importIntegerVector(ysrc,nshot,argv[8]);
    importIntegerVector(xrec,nrecs,argv[9]);
    importIntegerVector(yrec,nrecs,argv[10]);
    
    for (int shotPointer = 0; shotPointer < 1; shotPointer++)
    {
        setWaveField(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,nxx*nyy*nzz);

        for (int timePointer = 0; timePointer < nt; timePointer++) 
        {
            FDM_8E2T_elastic_Isotropic3D(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,rho,M,L,source,xsrc,ysrc,zsrc,nsrc,timePointer,shotPointer,nxx,nyy,nzz,dx,dy,dz,dt);        
            cerjanElasticAbsorbingCondition3D(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,damp,nxx,nyy,nzz);

            // Calcular ondas P
            // Calcular ondas Sv, Shx e Shy
            // Calcular campo escalar de pressão 

            // Calcular diferentes planos de snapshots XZ, YZ e XY
            // calcular sismograma 3D com a geometria de nodes
        }
    }
    
    t_f = time(NULL);
    total_time = difftime(t_f, t_0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);    
    return 0;
} 