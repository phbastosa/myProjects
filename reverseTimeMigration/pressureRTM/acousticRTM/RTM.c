# include <time.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "auxCodes/functions.h"

int main(int argc,char **argv) 
{
    float total_time;
    time_t t_0, t_f;
    t_0 = time(NULL);

    int nx,nz,nt,nsrc,nshot;
    int nabc,spread,nrecs; 
    float dx,dz,dt; 

    readParameters(&nx,&nz,&nt,&dx,&dz,&dt,&nabc,&spread,&nrecs,&nshot,&nsrc,argv[1]);

    int nxx = nx + 2*nabc;
    int nzz = nz + 2*nabc;

    int sampleInterval = 50;   /* Used to slice wave field */
    
    float *vels  = (float *) malloc(2*sizeof(float));
    float *vp    = (float *) malloc(nxx*nzz*sizeof(float));   /* Matrix of p wave velocities */
    float *P_pas = (float *) malloc(nxx*nzz*sizeof(float));   /* Matrix of past wave field */ 
    float *P_pre = (float *) malloc(nxx*nzz*sizeof(float));   /* Matrix of present wave field */
    float *P_fut = (float *) malloc(nxx*nzz*sizeof(float));   /* Matrix of future wave field */

    float *direct_field = (float *) malloc(nx*nz*(nt/sampleInterval)*sizeof(float)); /* Incident wave field matrix */ 
    float *direct_field_sum = (float *) malloc(nx*nz*sizeof(float));                 /* Summation of direct wave field */
    float *reverse_field_sum = (float *) malloc(nx*nz*sizeof(float));                /* Summation of reverse wave field */

    float *seism  = (float *) malloc(nt*spread*sizeof(float));  /* Seismogram matrix to use in migration */
    float *source = (float *) malloc(nsrc*sizeof(float));       /* Source Vector */
    
    float *image  = (float *) malloc(nx*nz*sizeof(float));      /* Final image matrix */
    float *aImage = (float *) malloc(nx*nz*sizeof(float));      /* Final image matrix */
    float *lImage = (float *) malloc(nx*nz*sizeof(float));      /* Final image matrix */
    
    int *xsrc = (int *) malloc(nshot*sizeof(int));              /* Horizontal shooting position array */
    int *xrec = (int *) malloc(nrecs*sizeof(int));              /* Horizontal receiver position array */    
    int *topo = (int *) malloc(nxx*sizeof(int));

    float *damp  = (float *) malloc(nxx*nzz*sizeof(float));  /* */

    importFloatVector(vp,nxx*nzz,argv[2]);
    importFloatVector(damp,nxx*nzz,argv[3]);
    importFloatVector(source,nsrc,argv[4]);
    
    importIntegerVector(xsrc,nshot,argv[5]);
    importIntegerVector(xrec,nrecs,argv[6]);

    ajustCoordinates(xrec,xsrc,topo,nabc,nxx,nrecs,nshot);

    getVelocities(nxx,nzz,vp,vels);

    memSet(image,nx*nz);                /* Zeroing the migrated image */
    memSet(direct_field_sum,nx*nz);     /* Zeroing summating direct field matrix */
    memSet(reverse_field_sum,nx*nz);    /* Zeroing summating reverse field matrix */

    FILE *read = fopen(argv[7],"rb");   /* Opening seismic data file */

    for(int shotPointer = 0; shotPointer < nshot; ++shotPointer) /* Shots loop */ 
    {
        memSet(P_pas,nxx*nzz);                          /* Zeroing past wave field */  
        memSet(P_pre,nxx*nzz);                          /* Zeroing present wave field */
        memSet(P_fut,nxx*nzz);                          /* Zeroing future wave field */
        memSet(seism,spread*nt);                        /* Zeroing seismogram per shot */
        memSet(direct_field,nx*nz*(nt/sampleInterval)); /* Zeroing wave fields per shot */ 

        fread(seism,sizeof(float),spread*nt,read);      /* Reading seismic data file */

        #pragma acc enter data copyin(P_pas[0:nxx*nzz],P_pre[0:nxx*nzz],P_fut[0:nxx*nzz],vp[0:nxx*nzz],seism[0:nt*spread],damp[0:nxx*nzz])
        #pragma acc enter data copyin(xsrc[0:nshot],xrec[0:nrecs],topo[0:nxx],source[0:nsrc],direct_field[0:nx*nz*(nt/sampleInterval)])
        #pragma acc enter data copyin(direct_field_sum[0:nx*nz],reverse_field_sum[0:nx*nz],image[0:nx*nz])
        for(int timePointer = 0; timePointer < nt; ++timePointer) /* Time loop */
        {
            propagationProgress(timePointer,shotPointer,xsrc,nshot,xrec,spread,dx,dz,nt,vels,dt,nxx,nzz,nabc);

            FDM_8E2T_acoustic2D(shotPointer,timePointer,vp,P_pre,P_pas,P_fut,source,nsrc,topo,xsrc,nshot,nxx,nzz,dx,dz,dt);
            
            cerjanAbsorbingBoundaryCondition(P_pas,P_pre,P_fut,damp,nxx*nzz);

            getSquareSumField(P_fut,direct_field_sum,nxx,nzz,nabc);

            getWaveField(timePointer,P_fut,direct_field,nxx,nzz,nt,nabc,sampleInterval);    
            
            waveFieldUpdate(P_pas,P_pre,P_fut,nxx*nzz);
        }

        waveFieldSet(P_pas,P_pre,P_fut,nxx*nzz);

        for(int timePointer = nt-1; timePointer >= 0; --timePointer) /* Temporal depropagation loop */
        {                  
            depropagationProgress(nt-timePointer,shotPointer,xsrc,nshot,xrec,spread,dx,dz,nt,vels,dt,nxx,nzz,nabc);

            FDM_8E2T_acoustic_depropagation(shotPointer,timePointer,vp,P_pre,P_pas,P_fut,seism,spread,xrec,nrecs,topo,nxx,nzz,nt,dx,dz,dt);

            getSquareSumField(P_fut,reverse_field_sum,nxx,nzz,nabc);

            crossCorrelation(timePointer,image,direct_field,P_fut,nxx,nzz,nabc,nt,sampleInterval);

            waveFieldUpdate(P_pas,P_pre,P_fut,nxx*nzz);
        }
        #pragma acc exit data delete(P_pas[0:nxx*nzz],P_pre[0:nxx*nzz],P_fut[0:nxx*nzz],vp[0:nxx*nzz],seism[0:nt*spread],damp[0:nxx*nzz])
        #pragma acc exit data delete(xsrc[0:nshot],xrec[0:nrecs],topo[0:nxx],source[0:nsrc],direct_field[0:nx*nz*(nt/sampleInterval)])
        #pragma acc exit data copyout(direct_field_sum[0:nx*nz],reverse_field_sum[0:nx*nz],image[0:nx*nz])
    }
    
    laplaciano(lImage,image,nx,nz,dx,dz);
    exportVector2D(lImage,nx*nz,"results/outputImage.bin"); /* writing the output image file */
    
    DsumCompensation(aImage,image,direct_field_sum,nx,nz);
    laplaciano(lImage,aImage,nx,nz,dx,dz);
    exportVector2D(lImage,nx*nz,"results/outputImageDsumComp.bin");
    
    RsumCompensation(aImage,image,reverse_field_sum,nx,nz);
    laplaciano(lImage,aImage,nx,nz,dx,dz);
    exportVector2D(lImage,nx*nz,"results/outputImageRsumComp.bin");

    sqrtDsumCompensation(aImage,image,direct_field_sum,nx,nz);
    laplaciano(lImage,aImage,nx,nz,dx,dz);    
    exportVector2D(lImage,nx*nz,"results/outputImageSqrtDsumComp.bin");   
    
    sqrtRsumCompensation(aImage,image,reverse_field_sum,nx,nz);
    laplaciano(lImage,aImage,nx,nz,dx,dz);        
    exportVector2D(lImage,nx*nz,"results/outputImageSqrtRsumComp.bin");
    
    DRsumCompensation(aImage,image,direct_field_sum,reverse_field_sum,nx,nz);
    laplaciano(lImage,aImage,nx,nz,dx,dz);            
    exportVector2D(lImage,nx*nz,"results/outputImageDRsumComp.bin");
    
    printf("\nThe migrated image was successfully written!\n");

    t_f = time(NULL);
    total_time = difftime(t_f, t_0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);

    return 0;
}
