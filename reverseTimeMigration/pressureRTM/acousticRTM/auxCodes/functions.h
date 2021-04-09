# ifndef FUNCTIONS_H_DEFINED
# define FUNCTIONS_H_DEFINED

void memSet(float * pointer, int n)
{
    for(int i = 0; i < n; i++) pointer[n] = 0.0f;
}

void importFloatVector(float * vector, int nPoints, char filename[])
{
    FILE * read = fopen((const char *) filename,"rb");
    fread(vector,sizeof(float),nPoints,read);
    fclose(read);
}

void importIntegerVector(int * vector, int nPoints, char filename[])
{
    FILE * read = fopen((const char *) filename,"rb");
    fread(vector,sizeof(int),nPoints,read);
    fclose(read);
}

void exportVector2D(float * vector, int nPoints, char filename[])
{
    FILE * write = fopen((const char *) filename, "wb");
    fwrite(vector, sizeof(float), nPoints, write);
    fclose(write);
}

void readParameters(int *nx,int *nz,int *nt,float *dx,float *dz,float *dt,int *abc,int *nrecs,int *nshot,int *nsrc,int *wbh, char filename[])
{
    FILE * arq = fopen((const char *) filename,"r"); 
    if(arq != NULL) 
    {
        fscanf(arq,"%i",nx); fscanf(arq,"%i",nz); fscanf(arq,"%i",nt); 
        fscanf(arq,"%f",dx); fscanf(arq,"%f",dz); fscanf(arq,"%f",dt); 
        fscanf(arq,"%i",abc); fscanf(arq,"%i",nrecs); fscanf(arq,"%i",nshot); 
        fscanf(arq,"%i",nsrc); fscanf(arq,"%i",wbh);  
    } 
    fclose(arq);
}

void getVelocities(int nxx, int nzz, float * vp,float *vels)
{
    float v_max = vp[0];
    float v_min = vp[0];

    for(int index = 0; index < nxx*nzz; ++index) 
    {
        int ii = floor(index / nxx);  /* line indicator */
        int jj = index % nxx;         /* column indicator */  

        if(v_max < vp[ii*nxx + jj]) v_max = vp[ii*nxx + jj]; 
        if(v_min > vp[ii*nxx + jj]) v_min = vp[ii*nxx + jj]; 
    }

    vels[0] = v_max; 
    vels[1] = v_min;    
}

/* */
void ajustCoordinates(int *xrec, int *xsrc, int *seaTop,int *seaBot, int wbh, int nabc, int nxx, int nrec, int nsrc)
{
    for (int ii = 0; ii < nrec; ii++) xrec[ii] += nabc;
    for (int ii = 0; ii < nsrc; ii++) xsrc[ii] += nabc;    

    for (int ii = 0; ii < nxx; ii++) seaTop[ii] = 4 + nabc;    
    for (int ii = 0; ii < nxx; ii++) seaBot[ii] = nabc + wbh;    
}

void printStatus(float * vels, float dx, float dz, float dt, int nxx, int nzz, int abc, int nt)
{
    printf("Reverse time migration\n\n");
    printf("Parameters:\n");
    printf("   Highest velocity in the model = %.1f m/s\n",vels[0]);
    printf("   Slowest velocity in the model = %.1f m/s\n",vels[1]);
    printf("   Spatial X discratization = %.1f meters\n",dx);
    printf("   Spatial Z discratization = %.1f meters\n",dz);
    printf("   Temporal discratization = %.4f seconds\n",dt);
    printf("   Total modeling time = %.2f seconds\n",dt*nt);
    printf("   Horizontal length of model = %.0f meters\n",dx*(nxx - 2*abc));
    printf("   Vertical length of model = %.0f meters\n\n",dz*(nzz - 2*abc));
}

void shotStatus(int tt, int * xsrc, int n_shot, int * xrec, int spread, float dx, 
                float dz, int nt, float * vels, float dt, int nxx, int nzz, int abc)
{
    system("clear");        
    printStatus(vels,dx,dz,dt,nxx,nzz,abc,nt);    
    printf("Migration status:\n");
    printf("   Shot position: %.1f meters\n",(xsrc[tt]-abc)*dx);
    printf("   Recivers position: %.1f - %.1f meters\n",(float)(xrec[0]-abc)*dx,(float)(xrec[spread-1]-abc)*dx);
    printf("   Total progress: %.2f %%\n",(float) tt/n_shot * 100.0f);
    printf("\nMigrated seismograms: %i of %i\n",tt,n_shot);
}   

void propagationProgress(int timePointer,int tt, int * xsrc, int n_shot, int * xrec, int spread, float dx, 
                         float dz, int nt, float * vels, float dt, int nxx, int nzz, int abc)
{
    if(timePointer % (nt / 100) == 0)
    {
        shotStatus(tt,xsrc,n_shot,xrec,spread,dx,dz,nt,vels,dt,nxx,nzz,abc);
        printf("\n    Propagation progress: %.2f %%\n",(float) timePointer/nt * 100.0f);
        printf("  Depropagation progress: 0.00 %%\n");
    }
}

void depropagationProgress(int timePointer,int tt, int * xsrc, int n_shot, int * xrec, int spread, float dx, 
                           float dz, int nt, float * vels, float dt, int nxx, int nzz, int abc)
{
    if(timePointer % (nt / 100) == 0)
    {    
        shotStatus(tt,xsrc,n_shot,xrec,spread,dx,dz,nt,vels,dt,nxx,nzz,abc);
        printf("\n    Propagation progress: 100.0 %%\n");
        printf("  Depropagation progress: %.2f %%\n",(float) timePointer/nt * 100.0f);
    }
}

void FDM_8E2T_acoustic2D(int shotPointer, int timePointer, float *vp, float *P_pre, float *P_pas, float *P_fut, float *source, int nsrc, int *topo, int *xsrc,int nshot,int nxx,int nzz, float dx, float dz, float dt)
{
    #pragma acc parallel loop present(P_pas[0:nxx*nzz],P_pre[0:nxx*nzz],P_fut[0:nxx*nzz],vp[0:nxx*nzz],topo[0:nxx],xsrc[0:nshot],source[0:nsrc])
    for(int index = 0; index < nxx*nzz; ++index) /* Spatial loop */ 
    {
        int ii = floor(index / nxx);  /* Line indicator */
        int jj = index % nxx;         /* Column indicator */  
        
        if((timePointer < nsrc) && (index == 0))
        { 
            P_pre[topo[xsrc[shotPointer]]*nxx + xsrc[shotPointer]] += source[timePointer]; 
        }

        if((ii > 3) && (ii < nzz-4) && (jj > 3) && (jj < nxx-4)) 
        {
            float d2_Px2 = (- 9.0f*(P_pre[ii*nxx + (jj-4)] + P_pre[ii*nxx + (jj+4)])
                        +   128.0f*(P_pre[ii*nxx + (jj-3)] + P_pre[ii*nxx + (jj+3)])
                        -  1008.0f*(P_pre[ii*nxx + (jj-2)] + P_pre[ii*nxx + (jj+2)])
                        +  8064.0f*(P_pre[ii*nxx + (jj+1)] + P_pre[ii*nxx + (jj-1)])
                        - 14350.0f*(P_pre[ii*nxx + jj]))/(5040.0f*powf(dx,2));

            float d2_Pz2 = (- 9.0f*(P_pre[(ii-4)*nxx + jj] + P_pre[(ii+4)*nxx + jj])
                        +   128.0f*(P_pre[(ii-3)*nxx + jj] + P_pre[(ii+3)*nxx + jj])
                        -  1008.0f*(P_pre[(ii-2)*nxx + jj] + P_pre[(ii+2)*nxx + jj])
                        +  8064.0f*(P_pre[(ii-1)*nxx + jj] + P_pre[(ii+1)*nxx + jj])
                        - 14350.0f*(P_pre[ii*nxx + jj]))/(5040.0f*powf(dz,2));

            P_fut[ii*nxx + jj] = (powf(dt,2.0f)*powf(vp[ii*nxx + jj],2.0f)*(d2_Px2 + d2_Pz2)) + 2.0f*P_pre[ii*nxx + jj] - P_pas[ii*nxx + jj];
        }
    }
}

void cerjanAbsorbingBoundaryCondition(float *P_pas, float *P_pre, float *P_fut, float *damp, int nPoints)
{
    #pragma acc parallel loop present(P_pas[0:nPoints],P_pre[0:nPoints],P_fut[0:nPoints],damp[0:nPoints])
    for (int index = 0; index < nPoints; index++)
    {
        P_pas[index] *= damp[index];
        P_pre[index] *= damp[index];
        P_fut[index] *= damp[index];
    }
}

void FDM_8E2T_acoustic_depropagation(int shotPointer, int timePointer, float *vp, float *P_pre, float *P_pas, float *P_fut, float *seism, int spread, int *xrec, int nrecs, int *topo, int nxx, int nzz,int nt, float dx, float dz, float dt)
{
    #pragma acc parallel loop present(P_pre[0:nxx*nzz],topo[0:nxx],xrec[0:nrecs],seism[0:nt*spread])
    for(int index = 0; index < spread; ++index)
    {
        P_pre[topo[xrec[index]]*nxx + xrec[index]] += seism[timePointer*spread + index]; /* Applying source*/
    }

    #pragma acc parallel loop present(P_pas[0:nxx*nzz],P_pre[0:nxx*nzz],P_fut[0:nxx*nzz],vp[0:nxx*nzz])
    for(int index = 0; index < nxx*nzz; ++index) /* Spatial loop */ 
    {
        int ii = floor(index / nxx);  /* Line indicator */
        int jj = index % nxx;         /* Column indicator */  

        if((ii > 3) && (ii < nzz-4) && (jj > 3) && (jj < nxx-4)) 
        {
            float d2_Px2 = (- 9.0f*(P_pre[ii*nxx + (jj-4)] + P_pre[ii*nxx + (jj+4)])
                        +   128.0f*(P_pre[ii*nxx + (jj-3)] + P_pre[ii*nxx + (jj+3)])
                        -  1008.0f*(P_pre[ii*nxx + (jj-2)] + P_pre[ii*nxx + (jj+2)])
                        +  8064.0f*(P_pre[ii*nxx + (jj+1)] + P_pre[ii*nxx + (jj-1)])
                        - 14350.0f*(P_pre[ii*nxx + jj]))/(5040.0f*powf(dx,2));

            float d2_Pz2 = (- 9.0f*(P_pre[(ii-4)*nxx + jj] + P_pre[(ii+4)*nxx + jj])
                        +   128.0f*(P_pre[(ii-3)*nxx + jj] + P_pre[(ii+3)*nxx + jj])
                        -  1008.0f*(P_pre[(ii-2)*nxx + jj] + P_pre[(ii+2)*nxx + jj])
                        +  8064.0f*(P_pre[(ii-1)*nxx + jj] + P_pre[(ii+1)*nxx + jj])
                        - 14350.0f*(P_pre[ii*nxx + jj]))/(5040.0f*powf(dz,2));

            P_fut[ii*nxx + jj] = (powf(dt,2)*powf(vp[ii*nxx + jj],2)*(d2_Px2 + d2_Pz2)) + 2.0f*P_pre[ii*nxx + jj] - P_pas[ii*nxx + jj];
        }
    }
}

void getSquareSumField(float * waveField, float * waveFieldSum, int nxx, int nzz, int abc)
{
    #pragma acc parallel loop present(waveField[0:nxx*nzz],waveFieldSum[0:(nxx-2*abc)*(nzz-2*abc)])
    for(int index = 0; index < nxx*nzz; ++index)
    {
        int ii = floor(index / nxx);  /* line indicator */
        int jj = index % nxx;         /* column indicator */                              

        /* Collect the sum of wavefield to applies compensations to final image */
        if((ii >= abc) && (ii < nzz - abc) && (jj >= abc) && (jj < nxx - abc)) 
        {                            
            waveFieldSum[(ii-abc)*(nxx - 2*abc)  + (jj-abc)] += powf(waveField[ii*nxx + jj],2.0f);
        }
    }
}

void getWaveField(int kk,float * waveField,float * waveFieldToGet,int nxx,int nzz,int nt,int abc,int sampleInterval)
{
    if(kk % sampleInterval == 0)
    {   
        #pragma acc parallel loop present(waveField[0:nxx*nzz],waveFieldToGet[0:(nxx-2*abc)*(nzz-2*abc)*(nt/sampleInterval)]) 
        for(int index = 0; index < nxx*nzz; ++index)
        {
            int ii = floor(index / nxx);  /* line indicator */
            int jj = index % nxx;         /* column indicator */  

            /* Getting wave field to use in cross correlation */
            if((ii >= abc) && (ii < nzz - abc) && (jj >= abc) && (jj < nxx - abc)) 
            {
                waveFieldToGet[(kk/sampleInterval)*(nzz-2*abc)*(nxx-2*abc) + (ii-abc)*(nxx-2*abc) + (jj-abc)] = waveField[ii*nxx + jj];
            }
        }
    }
}

void waveFieldUpdate(float * pas, float * pre, float * fut, int nPoints)
{
    #pragma acc parallel loop present(pas[0:nPoints],pre[0:nPoints],fut[0:nPoints])
    for(int index = 0; index < nPoints; ++index)
    {
        pas[index] = pre[index];         
        pre[index] = fut[index];         
    }
}

void waveFieldSet(float * pas, float * pre, float * fut, int nPoints)
{
    #pragma acc parallel loop present(pas[0:nPoints],pre[0:nPoints],fut[0:nPoints])
    for(int index = 0; index < nPoints; ++index)
    {
        pas[index] = 0.0f;
        pre[index] = 0.0f;         
        fut[index] = 0.0f;         
    }
}

void crossCorrelation(int kk, float * image, float * directField, float * reverseField, int nxx, int nzz, int abc, int nt, int sampleInterval)
{
    #pragma acc parallel loop present(image[0:(nxx-2*abc)*(nzz-2*abc)],directField[0:(nxx-2*abc)*(nzz-2*abc)*(nt/sampleInterval)],reverseField[0:nxx*nzz])
    if(kk % sampleInterval == 0)
    {
        for(int index = 0; index < nxx*nzz; ++index)
        {
            int ii = floor(index / nxx);  /* line indicator */
            int jj = index % nxx;         /* column indicator */                              

            /* Cross Correlation image condition */
            if((ii >= abc) && (ii < nzz - abc) && (jj >= abc) && (jj < nxx - abc)) 
            {                            
                image[(ii-abc)*(nxx-2*abc) + (jj-abc)] += directField[(kk/sampleInterval)*(nxx-2*abc)*(nzz-2*abc) + (ii-abc)*(nxx-2*abc) + (jj-abc)] * reverseField[ii*nxx + jj];
            }
        }
    }   
}

/* */
void DsumCompensation(float * aImage,float * image, float * directFieldSum, int nx, int nz)
{    
    for(int index = 0; index < nx*nz; ++index)
    {
        int ii = floor(index / nx);  /* line indicator */
        int jj = index % nx;         /* column indicator */                              

        aImage[ii*nx + jj] = image[ii*nx + jj] / directFieldSum[ii*nx + jj];
    }
}

/* */
void RsumCompensation(float * aImage,float * image,float * reverseFieldSum,int nx, int nz)
{    
    for(int index = 0; index < nx*nz; ++index)
    {
        int ii = floor(index / nx);  /* line indicator */
        int jj = index % nx;         /* column indicator */                              

        aImage[ii*nx + jj] = image[ii*nx + jj] / reverseFieldSum[ii*nx + jj];
    }
}

/* */
void sqrtDsumCompensation(float *aImage,float * image, float * directFieldSum, int nx, int nz)
{
    for(int index = 0; index < nx*nz; ++index)
    {
        int ii = floor(index / nx);  /* line indicator */
        int jj = index % nx;         /* column indicator */                              

        aImage[ii*nx + jj] = image[ii*nx + jj] /  sqrt(directFieldSum[ii*nx + jj]);
    }
}

/* */
void sqrtRsumCompensation(float *aImage,float * image, float * reverseFieldSum, int nx, int nz)
{
    for(int index = 0; index < nx*nz; ++index)
    {
        int ii = floor(index / nx);  /* line indicator */
        int jj = index % nx;         /* column indicator */                              

        aImage[ii*nx + jj] = image[ii*nx + jj] / sqrt(reverseFieldSum[ii*nx + jj]);
    }
}

/* */
void DRsumCompensation(float *aImage, float * image, float * directFieldSum,float * reverseFieldSum, int nx, int nz)
{
    for(int index = 0; index < nx*nz; ++index)
    {
        int ii = floor(index / nx);  /* line indicator */
        int jj = index % nx;         /* column indicator */                              

        aImage[ii*nx + jj] = image[ii*nx + jj] / (directFieldSum[ii*nx + jj] * reverseFieldSum[ii*nx + jj]);
    }
}

/* */

void laplaciano(float *aImage, float *image, int nx, int nz, float dx,float dz)
{
    for(int index = 0; index < nx*nz; index++)
    {
        int ii = floor(index / nx);
        int jj = index % nx;

        // float d2_U_dx2 = (image[(ii-1)*nx + jj] - 2.0f*image[ii*nx + jj] + image[(ii+1)*nx + jj])/powf(dx,2.0f);
        float d2_U_dx2 = 0.0f;

        float d2_U_dz2 = (image[ii*nx + (jj-1)] - 2.0f*image[nx*ii + jj] + image[(ii+1)*nx + jj])/powf(dz,2.0f);

        aImage[ii*nx + jj] = d2_U_dx2 + d2_U_dz2; 
    }    
}

# endif
