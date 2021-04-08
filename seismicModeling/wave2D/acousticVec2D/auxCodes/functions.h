# ifndef FUNCTIONS_H_DEFINED
# define FUNCTIONS_H_DEFINED

/* Function that make the memmory block to zero */
void mem_set(float *v, int nPoints) 
{
    for(int index = 0; index < nPoints; ++index) v[index] = 0.0f;
}

/* */
void importFloatVector(float *vector, int nPoints, char filename[])
{
    FILE * read = fopen((const char *) filename,"rb");
    fread(vector,sizeof(float),nPoints,read);
    fclose(read);
}

/* */
void importIntegerVector(int *vector, int nPoints, char filename[])
{
    FILE * read = fopen((const char *) filename,"rb");
    fread(vector,sizeof(int),nPoints,read);
    fclose(read);
}

/* */
void wavefield_set(float *Vx, float *Vz, float *P, int nxx, int nzz) 
{
    mem_set(Vx,nxx*nzz);
    mem_set(Vz,nxx*nzz);
    mem_set(P,nxx*nzz);
}

/* */
void exportVector(float *vector, int nPoints, char filename[])
{
    FILE * write = fopen((const char *) filename, "wb");
    fwrite(vector, sizeof(float), nPoints, write);
    fclose(write);
}

/* */
void readParameters(int *nx, int *nz, int *nt, float *dx, float *dz, float *dt, int *abc, int *nrec, int *nshot, int *nsrc, int *wbh, char filename[])
{
    FILE * arq = fopen((const char *) filename,"r"); 
    if(arq != NULL) 
    {
        fscanf(arq,"%i",nx); fscanf(arq,"%i",nz); fscanf(arq,"%i",nt); 
        fscanf(arq,"%f",dx); fscanf(arq,"%f",dz); fscanf(arq,"%f",dt); 
        fscanf(arq,"%i",abc); fscanf(arq,"%i",nrec); fscanf(arq,"%i",nshot); 
        fscanf(arq,"%i",nsrc); fscanf(arq,"%i",wbh); 
    } 
    fclose(arq);
}

/* */
float * getVelocities(int nxx, int nzz, float *vp)
{
    float v_max = vp[0];
    float v_min = 99999;
    
    float * vels = (float *) malloc(2*sizeof(float));

    for(int index = 0; index < nxx*nzz; ++index) 
    {
        int ii = floor(index / nxx);  
        int jj = index % nxx;           

        if(v_max < vp[ii*nxx + jj]) v_max = vp[ii*nxx + jj];         
        if(v_min > vp[ii*nxx + jj]) v_min = vp[ii*nxx + jj]; 
    }

    vels[0] = v_max; vels[1] = v_min;    
    return vels;
}

/* */
void printStatus(float *vels, float dx, float dz, float dt, int nxx, int nzz, int abc, int nt)
{
    printf("Synthetic seismic acquisition\n\n");
    printf("Parameters:\n");
    printf("   Highest velocity in the model = %.1f m/s\n",vels[0]);
    printf("   Slowest velocity in the model = %.1f m/s\n",vels[1]);
    printf("   Spatial X discratization = %.3f m\n",dx);
    printf("   Spatial Z discratization = %.3f m\n",dz);
    printf("   Temporal discratization = %.4f s\n",dt);
    printf("   Total modeling time = %.2f s\n",dt*nt);
    printf("   Horizontal length of model = %.0f m\n",dx*(nxx - 2*abc));
    printf("   Vertical length of model = %.0f m\n\n",dz*(nzz - 2*abc));
}

/* */
void modelingStatus(int shot, int time, int *xsrc, int n_shot, int *xrec, int spread, float dx, 
                    float dz, int nt, float *vels, float dt, int nxx, int nzz, int abc)
{
    if (time % (nt / 100) == 0)
    {
        system("clear");        
        printStatus(vels,dx,dz,dt,nxx,nzz,abc,nt);    

        printf("------------------------------------------------------\n");
        printf("   Propagation progress: %.2f %%\n",(float) time/nt * 100.0f);
        printf("------------------------------------------------------\n\n");

        printf("Modeling status:\n");
        printf("   Shot position: %.1f meters\n",(xsrc[shot]-abc)*dx);
        printf("   Recivers position: %.1f - %.1f meters\n",(xrec[0]-abc)*dx,(xrec[spread-1]-abc)*dx);
        printf("   Total progress: %.2f %%\n",(float) shot/n_shot * 100.0f);
        printf("\nExported seismograms: %i of %i\n",shot,n_shot);
    }
}   

/* */
void ajustCoordinates(int *xrec, int *xsrc, int *seaTop, int *seaBot, int wbh, int nabc, int nx, int nrec, int nsrc)
{
    for (int ii = 0; ii < nrec; ii++) xrec[ii] += nabc;
    for (int ii = 0; ii < nsrc; ii++) xsrc[ii] += nabc;    

    for (int ii = 0; ii < nx+2*nabc; ii++) seaTop[ii] = nabc + 4;    
    for (int ii = 0; ii < nx+2*nabc; ii++) seaBot[ii] = nabc + wbh;    
}

void FDM8E2T_acousticVec2D(int shot, int time, float *Vx, float *Vz, float *P, float *K, float *b, float *source,
                           int nsrc, int *xsrc, int *topo, int nxx, int nzz, float dx, float dz, float dt, int nshot)
{
    # pragma acc parallel loop present(Vx[0:nxx*nzz],Vz[0:nxx*nzz],P[0:nxx*nzz],K[0:nxx*nzz],source[0:nsrc],topo[0:nxx],xsrc[0:nshot])
    for(int index = 0; index < nxx*nzz; index++) 
    {                  
        int ii = (int) index / nxx;      // indicador de linhas  (direção z)
        int jj = (int) index % nxx;      // indicador de colunas (direção x)

        if((index == 0) && (time < nsrc)) 
        {
            P[topo[xsrc[shot]]*nxx + xsrc[shot]] += source[time] / (dx*dz);
        }

        if((ii >= 3) && (ii < nzz-3) && (jj >= 3) && (jj < nxx-4)) 
        {
            float dVx_dx = (75.0f*(Vx[(jj-3) + ii*nxx] - Vx[(jj+4) + ii*nxx]) + 
                          1029.0f*(Vx[(jj+3) + ii*nxx] - Vx[(jj-2) + ii*nxx]) +
                          8575.0f*(Vx[(jj-1) + ii*nxx] - Vx[(jj+2) + ii*nxx]) + 
                        128625.0f*(Vx[(jj+1) + ii*nxx] - Vx[jj + ii*nxx]))/(dx*107520.0f);

            float dVz_dz = (75.0f*(Vz[jj + (ii-3)*nxx] - Vz[jj + (ii+4)*nxx]) +   
                          1029.0f*(Vz[jj + (ii+3)*nxx] - Vz[jj + (ii-2)*nxx]) +
                          8575.0f*(Vz[jj + (ii-1)*nxx] - Vz[jj + (ii+2)*nxx]) +
                        128625.0f*(Vz[jj + (ii+1)*nxx] - Vz[jj + ii*nxx]))/(dz*107520.0f);     

            P[index] += -dt*K[index]*(dVx_dx + dVz_dz);
        }
    }

    # pragma acc parallel loop present(Vx[0:nxx*nzz],Vz[0:nxx*nzz],P[0:nxx*nzz],b[0:nxx*nzz])
    for(int index = 0; index < nxx*nzz; index++) 
    {                  
        int ii = (int) index / nxx;      // indicador de linhas  (direção z)
        int jj = (int) index % nxx;      // indicador de colunas (direção x)
       
        if((ii >= 0) && (ii < nzz) && (jj > 3) && (jj < nxx-3)) 
        {
            float dP_dx = (75.0f*(P[(jj-4) + ii*nxx] - P[(jj+3) + ii*nxx]) +
                         1029.0f*(P[(jj+2) + ii*nxx] - P[(jj-3) + ii*nxx]) +
                         8575.0f*(P[(jj-2) + ii*nxx] - P[(jj+1) + ii*nxx]) +
                       128625.0f*(P[jj + ii*nxx]     - P[(jj-1) + ii*nxx]))/(dx*107520.0f);

            float bx = 0.5f*(b[jj+1 + ii*nxx] + b[jj + ii*nxx]);

            Vx[index] += -dt*bx*dP_dx;
        }

        if((ii > 3) && (ii < nzz-3) && (jj >= 0) && (jj < nxx)) 
        {
            float dP_dz = (75.0f*(P[jj + (ii-4)*nxx] - P[jj + (ii+3)*nxx]) + 
                         1029.0f*(P[jj + (ii+2)*nxx] - P[jj + (ii-3)*nxx]) +
                         8575.0f*(P[jj + (ii-2)*nxx] - P[jj + (ii+1)*nxx]) +
                       128625.0f*(P[jj + ii*nxx]     - P[jj + (ii-1)*nxx]))/(dz*107520.0f);

            float bz = 0.5f*(b[jj + (ii+1)*nxx] + b[jj + ii*nxx]);

            Vz[index] += -dt*bz*dP_dz;
        }
    }
}

/* */
void cerjanAbsorbingBoundaryCondition(float *Vx, float *Vz, float *P, float *damp, int nPoints)
{
    #pragma acc parallel loop present(Vx[0:nPoints],Vz[0:nPoints],P[0:nPoints],damp[0:nPoints])
    for (int index = 0; index < nPoints; index++)
    {
        Vx[index] *= damp[index];
        Vz[index] *= damp[index];
        P[index] *= damp[index];
    }
}

/* */
void getSeismogram(float *seism, float *P, int *xrec, int *seaBot, int nrec, int nxx, int nzz, int nt, int shotPointer, int timePointer)
{
    #pragma acc parallel loop present(seism[0:nt*nrec],P[0:nxx*nzz],xrec[0:nrec],seaBot[0:nxx])
    for (int ii = 0; ii < nrec; ii++)
    {
        seism[timePointer*nrec + ii] = P[seaBot[xrec[ii]]*nxx + xrec[ii]];
    }
}

/* */
void joiningSeismograms(float *seism, float *seismograms, int spread, int nt, int shotPointer)
{
    for (int s = 0; s < spread; s++)
    {
        for (int t = 0; t < nt; t++)
        {
            seismograms[shotPointer*spread*nt + t*spread + s] = seism[t*spread + s];        
        }
    }
}

# endif
