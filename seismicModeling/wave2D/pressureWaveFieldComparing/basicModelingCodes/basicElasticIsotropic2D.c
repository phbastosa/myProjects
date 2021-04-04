# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "../auxCodes/funtions.h"

int main(int argc, char **argv) 
{
    printf("\n\nBasic elastic isotropic modeling\n");

    float totalTime;
    time_t t0, tf;
    t0 = time(NULL);

    int nx, nz, nt, nsrc;
    int xsrc, zsrc, zrec;
    float dx, dz, dt;

    readParameters(&nx,&nz,&nt,&dx,&dz,&dt,&nsrc,&xsrc,&zsrc,&zrec,argv[1]);

    float *seismogram = (float *) malloc(nx*nt*sizeof(float));

    float *vp  = (float *) malloc(nx*nz*sizeof(float));  
    float *vs  = (float *) malloc(nx*nz*sizeof(float));
    float *rho = (float *) malloc(nx*nz*sizeof(float));
    float *M   = (float *) malloc(nx*nz*sizeof(float));
    float *L   = (float *) malloc(nx*nz*sizeof(float));
    float *Vx  = (float *) malloc(nx*nz*sizeof(float));
    float *Vz  = (float *) malloc(nx*nz*sizeof(float));
    float *Txx = (float *) malloc(nx*nz*sizeof(float));
    float *Tzz = (float *) malloc(nx*nz*sizeof(float));
    float *Txz = (float *) malloc(nx*nz*sizeof(float));
    
    float *source = (float *) malloc(nsrc*sizeof(float));
    
    for(int index = 0; index < nx*nz; index++) 
    {
        vp[index]  = 1500.0f;
        vs[index]  = 0.0f;
        rho[index] = 1000.0f;
        M[index]   = rho[index]*powf(vs[index],2.0f);
        L[index]   = rho[index]*powf(vp[index],2.0f) - 2.0f*M[index];
    }

    importFloatVector(source,nsrc,argv[2]);

    # pragma acc enter data copyin(Vx[0:nx*nz],Vz[0:nx*nz],Txx[0:nx*nz],Tzz[0:nx*nz],Txz[0:nx*nz])
    # pragma acc enter data copyin(rho[0:nx*nz],M[0:nx*nz],L[0:nx*nz],source[0:nsrc],seismogram[0:nx*nt])
    for(int timePointer = 0; timePointer < nt; timePointer++) 
    {    
        if(timePointer % 200 == 0) printf("Propagation time = %0.5f seconds\n", (timePointer+200)*dt);
                  
        FDM8E2T_stressStencil_elasticIsotropic2D(Vx,Vz,Txx,Tzz,Txz,rho,M,L,nx,nz,dt,dx,dz,timePointer,source,nsrc,zsrc,xsrc);
        // FDM8E2T_velocityStencil_elasticIsotropic2D(Vx,Vz,Txx,Tzz,Txz,rho,M,L,nx,nz,dt,dx,dz,timePointer,source,nsrc,zsrc,xsrc);

        getElasticIsotropicPressureSeismogram(seismogram,Txx,Tzz,nt,nx,nz,timePointer,zrec);
    }
    # pragma acc exit data delete(Vx[0:nx*nz],Vz[0:nx*nz],Txx[0:nx*nz],Tzz[0:nx*nz],Txz[0:nx*nz])
    # pragma acc exit data delete(rho[0:nx*nz],M[0:nx*nz],L[0:nx*nz],source[0:nsrc])
    # pragma acc exit data copyout(seismogram[0:nx*nt])

    exportVector(seismogram,nx*nt,"results/seismograms/pressureElasticIsotropic2D.bin");

    tf = time(NULL);
    totalTime = difftime(tf, t0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", totalTime);

    return 0;
} 
