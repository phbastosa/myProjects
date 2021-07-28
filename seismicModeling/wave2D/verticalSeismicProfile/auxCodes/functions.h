# ifndef FUNCTIONS_H_DEFINED
# define FUNCTIONS_H_DEFINED

/* */
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
void wavefield_set(float *Vx, float *Vz, float *Txx, float *Tzz, float *Txz, int nxx, int nzz) 
{
    mem_set(Vx,nxx*nzz);
    mem_set(Vz,nxx*nzz);
    mem_set(Txx,nxx*nzz);
    mem_set(Tzz,nxx*nzz);
    mem_set(Txz,nxx*nzz);
}

/* */
void exportVector(float *vector, int nPoints, char filename[])
{
    FILE * write = fopen((const char *) filename, "wb");
    fwrite(vector, sizeof(float), nPoints, write);
    fclose(write);
}

/* */
void readParameters(int *nx, int *nz, int *nt, float *dx, float *dz, float *dt, int *abc, int *xsrc, int *zsrc, int *nsrc, char filename[])
{
    FILE * arq = fopen((const char *) filename,"r"); 
    if(arq != NULL) 
    {
        fscanf(arq,"%i",nx); fscanf(arq,"%i",nz); fscanf(arq,"%i",nt); 
        fscanf(arq,"%f",dx); fscanf(arq,"%f",dz); fscanf(arq,"%f",dt); 
        fscanf(arq,"%i",abc); fscanf(arq,"%i",xsrc); fscanf(arq,"%i",zsrc); 
        fscanf(arq,"%i",nsrc);  
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
void ajustCoordinates(int *xrec, int *xsrc, int *seaTop,int *seaBot, int wbh, int nabc, int nx, int nrec, int nsrc)
{
    for (int ii = 0; ii < nrec; ii++) xrec[ii] += nabc;
    for (int ii = 0; ii < nsrc; ii++) xsrc[ii] += nabc;    

    for (int ii = 0; ii < nx+2*nabc; ii++) seaTop[ii] = 4 + nabc;    
    for (int ii = 0; ii < nx+2*nabc; ii++) seaBot[ii] = wbh + nabc;    
}

/* */
void FDM8E2T_stressStencil_elasticIsotropic2D(float *Vx, float *Vz, float *Txx, float *Tzz, float *Txz, float *rho, float *M, float *L, int nxx, int nzz, float dt, float dx, float dz, int zsrc, int xsrc, float *source, int nsrc, int timePointer) 
{   
    #pragma acc parallel loop present(Vx[0:nxx*nzz],Vz[0:nxx*nzz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz],Txz[0:nxx*nzz],L[0:nxx*nzz],M[0:nxx*nzz],source[0:nsrc])
    for(int index = 0; index < nxx*nzz; index++) 
    {                  
        int ii = (int) index / nxx;      // indicador de linhas  (direção z)
        int jj = (int) index % nxx;      // indicador de colunas (direção x)

        if((timePointer < nsrc) && (index == 0))
        {
            Txx[zsrc*nxx + xsrc] += source[timePointer] / (dx*dz); 
            Tzz[zsrc*nxx + xsrc] += source[timePointer] / (dx*dz); 
        }

        if((ii >= 3) && (ii < nzz-4) && (jj >= 3) && (jj < nxx-4)) 
        {           
            float dVx_dx = (75.0f*(Vx[(jj-3) + ii*nxx] - Vx[(jj+4) + ii*nxx]) + 
                          1029.0f*(Vx[(jj+3) + ii*nxx] - Vx[(jj-2) + ii*nxx]) +
                          8575.0f*(Vx[(jj-1) + ii*nxx] - Vx[(jj+2) + ii*nxx]) + 
                        128625.0f*(Vx[(jj+1) + ii*nxx] - Vx[jj + ii*nxx]))/(dx*107520.0f);

            float dVz_dz = (75.0f*(Vz[jj + (ii-3)*nxx] - Vz[jj + (ii+4)*nxx]) +   
                          1029.0f*(Vz[jj + (ii+3)*nxx] - Vz[jj + (ii-2)*nxx]) +
                          8575.0f*(Vz[jj + (ii-1)*nxx] - Vz[jj + (ii+2)*nxx]) +
                        128625.0f*(Vz[jj + (ii+1)*nxx] - Vz[jj + ii*nxx]))/(dz*107520.0f);     

            Txx[index] += dt*((L[index] + 2.0f*M[index])*dVx_dx + L[index]*dVz_dz);   
            Tzz[index] += dt*((L[index] + 2.0f*M[index])*dVz_dz + L[index]*dVx_dx);
        }
    
        if((ii > 3) && (ii < nzz-3) && (jj > 3) && (jj < nxx-3)) 
        {
            float dVz_dx = (75.0f*(Vz[(jj-4) + ii*nxx] - Vz[(jj+3) + ii*nxx]) +
                          1029.0f*(Vz[(jj+2) + ii*nxx] - Vz[(jj-3) + ii*nxx]) +
                          8575.0f*(Vz[(jj-2) + ii*nxx] - Vz[(jj+1) + ii*nxx]) +
                        128625.0f*(Vz[jj + ii*nxx]     - Vz[(jj-1) + ii*nxx]))/(dx*107520.0f);

            float dVx_dz = (75.0f*(Vx[jj + (ii-4)*nxx] - Vx[jj + (ii+3)*nxx]) +
                          1029.0f*(Vx[jj + (ii+2)*nxx] - Vx[jj + (ii-3)*nxx]) +
                          8575.0f*(Vx[jj + (ii-2)*nxx] - Vx[jj + (ii+1)*nxx]) +
                        128625.0f*(Vx[jj + ii*nxx]     - Vx[jj + (ii-1)*nxx]))/(dz*107520.0f);

            float M_int = powf(0.25*(1.0f/M[jj + (ii+1)*nxx] + 1.0f/M[(jj+1) + ii*nxx] + 1.0f/M[(jj+1) + (ii+1)*nxx] + 1.0f/M[jj + ii*nxx]),-1); 

            Txz[index] += dt*M_int*(dVx_dz + dVz_dx);            
        }          
    }

    #pragma acc parallel loop present(Vx[0:nxx*nzz],Vz[0:nxx*nzz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz],Txz[0:nxx*nzz],rho[0:nxx*nzz])
    for(int index = 0; index < nxx*nzz; index++) 
    {              
        int ii = (int) index / nxx;      // indicador de linhas  (direção z)
        int jj = (int) index % nxx;      // indicador de colunas (direção x)

        if((ii >= 3) && (ii < nzz-4) && (jj > 3) && (jj < nxx-3)) 
        {
            float dTxx_dx = (75.0f*(Txx[(jj-4) + ii*nxx] - Txx[(jj+3) + ii*nxx]) +
                           1029.0f*(Txx[(jj+2) + ii*nxx] - Txx[(jj-3) + ii*nxx]) +
                           8575.0f*(Txx[(jj-2) + ii*nxx] - Txx[(jj+1) + ii*nxx]) +
                         128625.0f*(Txx[jj + ii*nxx]     - Txx[(jj-1) + ii*nxx]))/(dx*107520.0f);

            float dTxz_dz = (75.0f*(Txz[jj + (ii-3)*nxx] - Txz[jj + (ii+4)*nxx]) +
                           1029.0f*(Txz[jj + (ii+3)*nxx] - Txz[jj + (ii-2)*nxx]) + 
                           8575.0f*(Txz[jj + (ii-1)*nxx] - Txz[jj + (ii+2)*nxx]) +
                         128625.0f*(Txz[jj + (ii+1)*nxx] - Txz[jj + ii*nxx]))/(dz*107520.0f);

            float rho_int = 0.5f*(rho[(jj+1) + ii*nxx] + rho[jj + ii*nxx]);

            Vx[index] += dt/rho_int*(dTxx_dx + dTxz_dz);  
        }
      
        if((ii > 3) && (ii < nzz-3) && (jj >= 3) && (jj < nxx-4)) 
        {
            float dTxz_dx = (75.0f*(Txz[(jj-3) + ii*nxx] - Txz[(jj+4) + ii*nxx]) +
                           1029.0f*(Txz[(jj+3) + ii*nxx] - Txz[(jj-2) + ii*nxx]) +
                           8575.0f*(Txz[(jj-1) + ii*nxx] - Txz[(jj+2) + ii*nxx]) +
                         128625.0f*(Txz[(jj+1) + ii*nxx] - Txz[jj + ii*nxx]))/(dx*107520.0f);

            float dTzz_dz = (75.0f*(Tzz[jj + (ii-4)*nxx] - Tzz[jj + (ii+3)*nxx]) + 
                           1029.0f*(Tzz[jj + (ii+2)*nxx] - Tzz[jj + (ii-3)*nxx]) +
                           8575.0f*(Tzz[jj + (ii-2)*nxx] - Tzz[jj + (ii+1)*nxx]) +
                         128625.0f*(Tzz[jj + ii*nxx]     - Tzz[jj + (ii-1)*nxx]))/(dz*107520.0f);

            float rho_int = 0.5f*(rho[jj + (ii+1)*nxx] + rho[jj + ii*nxx]);

            Vz[index] += dt/rho_int*(dTxz_dx + dTzz_dz); 
        }
    }
}

/* */
void cerjanAbsorbingBoundaryCondition(float *Vx, float *Vz, float *Txx, float *Tzz, float *Txz, float *damp, int nxx, int nzz)
{
    #pragma acc parallel loop present(Vx[0:nxx*nzz],Vz[0:nxx*nzz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz],Txz[0:nxx*nzz],damp[0:nxx*nzz])
    for (int index = 0; index < nxx*nzz; index++)
    {
        Vx[index] *= damp[index];
        Vz[index] *= damp[index];
        Txx[index] *= damp[index];
        Tzz[index] *= damp[index];
        Txz[index] *= damp[index];
    }
}

/* */
void getSeismogram(float *seism, float *Txx, float *Tzz, int xsrc, int nz, int nabc, int nxx, int nzz, int nt, int timePointer)
{
    #pragma acc parallel loop present(seism[0:nt*nz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz])
    for (int ii = nabc; ii < nz + nabc; ii++)
    {
        seism[timePointer*nz + ii-nabc] = (Txx[ii*nxx + xsrc] + Tzz[ii*nxx + xsrc]) / 2.0f;
    }
}

# endif
