# ifndef FUNCTIONS_H_DEFINED
# define FUNCTIONS_H_DEFINED

void importFloatVector(float * vector, int nPoints, char filename[])
{
    FILE * read = fopen((const char *) filename,"rb");
    fread(vector,sizeof(float),nPoints,read);
    fclose(read);
}

void exportVector(float * vector, int nPoints, char filename[])
{
    FILE * write = fopen((const char *) filename, "wb");
    fwrite(vector, sizeof(float), nPoints, write);
    fclose(write);
}

void readParameters(int *nx, int *nz, int *nt, float *dx, float *dz, float *dt, int *nsrc, int *xsrc, int *zsrc, int *zrec, char filename[])
{
    FILE * arq = fopen((const char *) filename,"r"); 
    if(arq != NULL) 
    {
        fscanf(arq,"%i",nx); fscanf(arq,"%i",nz); fscanf(arq,"%i",nt); 
        fscanf(arq,"%f",dx); fscanf(arq,"%f",dz); fscanf(arq,"%f",dt); 
        fscanf(arq,"%i",nsrc); fscanf(arq,"%i",xsrc); fscanf(arq,"%i",zsrc); 
        fscanf(arq,"%i",zrec);  
    } 
    fclose(arq);
}

void FDM8E2T_acoustic2D(int timePointer, float *vp, float *P_pre, float *P_pas, float *P_fut, float *source, 
                        int nsrc, int xsrc, int zsrc, int nxx, int nzz, float dx, float dz, float dt)
{
    # pragma acc parallel loop present(vp[0:nxx*nzz],P_pas[0:nxx*nzz],P_pre[0:nxx*nzz],P_fut[0:nxx*nzz],source[0:nsrc])
    for(int index = 0; index < nxx*nzz; ++index) 
    {
        int ii = floor(index / nxx);  /* Line indicator */
        int jj = index % nxx;         /* Column indicator */  
        
        if((timePointer < nsrc) && (index == 0))
        { 
            P_pre[zsrc*nxx + xsrc] += source[timePointer] / (dx*dz); 
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

            P_fut[ii*nxx + jj] = powf(dt,2.0f)*powf(vp[ii*nxx + jj],2.0f) * (d2_Px2 + d2_Pz2) 
                               + 2.0f*P_pre[ii*nxx + jj] - P_pas[ii*nxx + jj];
        }
    }
}

void waveFieldUpdate(float * pas, float * pre, float * fut, int nPoints)
{
    # pragma acc parallel loop present(pas[0:nPoints],pre[0:nPoints],fut[0:nPoints])
    for(int index = 0; index < nPoints; ++index)
    {
        pas[index] = pre[index];         
        pre[index] = fut[index];         
    }
}

void FDM8E2T_stressStencil_elasticIsotropic2D(float *Vx, float *Vz, float *Txx, float *Tzz, float *Txz, float *rho, 
                                              float *M, float *L, int nxx, int nzz, float dt, float dx, float dz, 
                                              int timePointer, float *source, int nsrc, int zsrc, int xsrc) 
{   
    # pragma acc parallel loop present(Vx[0:nxx*nzz],Vz[0:nxx*nzz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz],Txz[0:nxx*nzz],L[0:nxx*nzz],M[0:nxx*nzz],source[0:nsrc])
    for(int index = 0; index < nxx*nzz; index++) 
    {                  
        int ii = (int) index / nxx;      // indicador de linhas  (direção z)
        int jj = (int) index % nxx;      // indicador de colunas (direção x)

        if((index == 0) && (timePointer < nsrc)) 
        {
            Txx[zsrc*nxx + xsrc] += source[timePointer] / (dx*dz);
            Tzz[zsrc*nxx + xsrc] += source[timePointer] / (dx*dz);   
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

            float M_int = powf(0.25f*(1.0f/M[jj + (ii+1)*nxx] + 1.0f/M[(jj+1) + ii*nxx] + 1.0f/M[(jj+1) + (ii+1)*nxx] + 1.0f/M[jj + ii*nxx]),-1.0f); 

            Txz[index] += dt*M_int*(dVx_dz + dVz_dx);            
        }          
    }

    # pragma acc parallel loop present(Vx[0:nxx*nzz],Vz[0:nxx*nzz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz],Txz[0:nxx*nzz],rho[0:nxx*nzz])
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

void FDM8E2T_velocityStencil_elasticIsotropic2D(float *Vx, float *Vz, float *Txx, float *Tzz, float *Txz, float *rho, 
                                                float *M, float *L, int nxx, int nzz, float dt, float dx, float dz,
                                                int timePointer, float *source, int nsrc, int zsrc, int xsrc) 
{
    # pragma acc parallel loop present(Vx[0:nxx*nzz],Vz[0:nxx*nzz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz],Txz[0:nxx*nzz],rho[0:nxx*nzz])
    for(int index = 0; index < nxx*nzz; index++) 
    {              
        int ii = (int) index / nxx;      // indicador de linhas  (direção z)
        int jj = (int) index % nxx;      // indicador de colunas (direção x)

        if((ii >= 3) && (ii < nzz-4) && (jj >= 3) && (jj < nxx-4)) 
        {
            float dTxx_dx = (75.0f*(Txx[(jj-3) + ii*nxx] - Txx[(jj+4) + ii*nxx]) +
                           1029.0f*(Txx[(jj+3) + ii*nxx] - Txx[(jj-2) + ii*nxx]) +
                           8575.0f*(Txx[(jj-1) + ii*nxx] - Txx[(jj+2) + ii*nxx]) +
                         128625.0f*(Txx[(jj+1) + ii*nxx] - Txx[jj + ii*nxx]))/(dx*107520.0f);

            float dTxz_dz = (75.0f*(Txz[jj + (ii-3)*nxx] - Txz[jj + (ii+4)*nxx]) +
                           1029.0f*(Txz[jj + (ii+3)*nxx] - Txz[jj + (ii-2)*nxx]) + 
                           8575.0f*(Txz[jj + (ii-1)*nxx] - Txz[jj + (ii+2)*nxx]) +
                         128625.0f*(Txz[jj + (ii+1)*nxx] - Txz[jj + ii*nxx]))/(dz*107520.0f);

            Vx[index] += dt/rho[index]*(dTxx_dx + dTxz_dz);  
        }
      
        if((ii > 3) && (ii < nzz-3) && (jj > 3) && (jj < nxx-3)) 
        {
            float dTxz_dx = (75.0f*(Txz[(jj-4) + ii*nxx] - Txz[(jj+3) + ii*nxx]) +
                           1029.0f*(Txz[(jj+2) + ii*nxx] - Txz[(jj-3) + ii*nxx]) +
                           8575.0f*(Txz[(jj-2) + ii*nxx] - Txz[(jj+1) + ii*nxx]) +
                         128625.0f*(Txz[jj + ii*nxx]     - Txz[(jj-1) + ii*nxx]))/(dx*107520.0f);

            float dTzz_dz = (75.0f*(Tzz[jj + (ii-4)*nxx] - Tzz[jj + (ii+3)*nxx]) + 
                           1029.0f*(Tzz[jj + (ii+2)*nxx] - Tzz[jj + (ii-3)*nxx]) +
                           8575.0f*(Tzz[jj + (ii-2)*nxx] - Tzz[jj + (ii+1)*nxx]) +
                         128625.0f*(Tzz[jj + ii*nxx]     - Tzz[jj + (ii-1)*nxx]))/(dz*107520.0f);

            float rho_int = powf(0.25f*(1.0f/rho[jj + ii*nxx] + 1.0f/rho[jj + (ii+1)*nxx] + 1.0f/rho[(jj+1) + (ii+1)*nxx] + 1.0f/rho[(jj+1) + ii*nxx]),-1.0f);

            Vz[index] += dt/rho_int*(dTxz_dx + dTzz_dz); 
        }
    }

    # pragma acc parallel loop present(Vx[0:nxx*nzz],Vz[0:nxx*nzz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz],Txz[0:nxx*nzz],L[0:nxx*nzz],M[0:nxx*nzz],source[0:nsrc])
    for(int index = 0; index < nxx*nzz; index++) 
    {                  
        int ii = (int) index / nxx;      // indicador de linhas  (direção z)
        int jj = (int) index % nxx;      // indicador de colunas (direção x)

        if((index == 0) && (timePointer < nsrc)) 
        {
            Txx[zsrc*nxx + xsrc] += source[timePointer] / (dx*dz);
            Tzz[zsrc*nxx + xsrc] += source[timePointer] / (dx*dz);   
        }

        if((ii >= 3) && (ii < nzz-4) && (jj > 3) && (jj < nxx-3)) 
        {
            float dVx_dx = (75.0f*(Vx[(jj-4) + ii*nxx] - Vx[(jj+3) + ii*nxx]) + 
                          1029.0f*(Vx[(jj+2) + ii*nxx] - Vx[(jj-3) + ii*nxx]) +
                          8575.0f*(Vx[(jj-2) + ii*nxx] - Vx[(jj+1) + ii*nxx]) + 
                        128625.0f*(Vx[jj + ii*nxx]     - Vx[(jj-1) + ii*nxx]))/(dx*107520.0f);

            float dVz_dz = (75.0f*(Vz[jj + (ii-3)*nxx] - Vz[jj + (ii+4)*nxx]) +   
                          1029.0f*(Vz[jj + (ii+3)*nxx] - Vz[jj + (ii-2)*nxx]) +
                          8575.0f*(Vz[jj + (ii-1)*nxx] - Vz[jj + (ii+2)*nxx]) +
                        128625.0f*(Vz[jj + (ii+1)*nxx] - Vz[jj + ii*nxx]))/(dz*107520.0f);     

            float L_int = 0.5f*(L[(jj+1) + ii*nxx] + L[jj + ii*nxx]);
            float M_int = 0.5f*(M[(jj+1) + ii*nxx] + M[jj + ii*nxx]);

            Txx[index] += dt*((L_int + 2.0f*M_int)*dVx_dx + L_int*dVz_dz);   
            Tzz[index] += dt*((L_int + 2.0f*M_int)*dVz_dz + L_int*dVx_dx);
        }

        if((ii >= 3) && (ii < nzz-4) && (jj >= 3) && (jj < nxx-4)) 
        {
            float dVz_dx = (75.0f*(Vz[(jj-3) + ii*nxx] - Vz[(jj+4) + ii*nxx]) +
                          1029.0f*(Vz[(jj+3) + ii*nxx] - Vz[(jj-2) + ii*nxx]) +
                          8575.0f*(Vz[(jj-1) + ii*nxx] - Vz[(jj+2) + ii*nxx]) +
                        128625.0f*(Vz[(jj+1) + ii*nxx] - Vz[jj + ii*nxx]))/(dx*107520.0f);

            float dVx_dz = (75.0f*(Vx[jj + (ii-4)*nxx] - Vx[jj + (ii+3)*nxx]) +
                          1029.0f*(Vx[jj + (ii+2)*nxx] - Vx[jj + (ii-3)*nxx]) +
                          8575.0f*(Vx[jj + (ii-2)*nxx] - Vx[jj + (ii+1)*nxx]) +
                        128625.0f*(Vx[jj + ii*nxx]     - Vx[jj + (ii-1)*nxx]))/(dz*107520.0f);

            float M_int = 0.5f*(M[jj + (ii+1)*nxx] + M[jj + ii*nxx]); 

            Txz[index] += dt*M_int*(dVx_dz + dVz_dx);            
        }      
    }
}

void FDM8E2T_shearStencil_elasticIsotropic2D(float *Vx, float *Vz, float *Txx, float *Tzz, float *Txz, float *rho, 
                                              float *M, float *L, int nxx, int nzz, float dt, float dx, float dz,
                                              int timePointer, float *source, int nsrc, int zsrc, int xsrc) 
{
    # pragma acc parallel loop present(Vx[0:nxx*nzz],Vz[0:nxx*nzz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz],Txz[0:nxx*nzz],L[0:nxx*nzz],M[0:nxx*nzz],source[0:nsrc])
    for(int index = 0; index < nxx*nzz; index++) 
    {                  
        int ii = (int) index / nxx;      /* Line indicator */
        int jj = (int) index % nxx;      /* Column indicator */

        if((timePointer < nsrc) && (index == 0))
        {
            Txx[zsrc*nxx + xsrc] += source[timePointer] / (dx*dz); 
            Tzz[zsrc*nxx + xsrc] += source[timePointer] / (dx*dz); 
        }

        if((ii > 3) && (ii < nzz-3) && (jj > 3) && (jj < nxx-3)) 
        {
            float dVx_dx = (75.0f*(Vx[(jj-4) + ii*nxx] - Vx[(jj+3) + ii*nxx]) + 
                          1029.0f*(Vx[(jj+2) + ii*nxx] - Vx[(jj-3) + ii*nxx]) +
                          8575.0f*(Vx[(jj-2) + ii*nxx] - Vx[(jj+1) + ii*nxx]) + 
                        128625.0f*(Vx[jj + ii*nxx]     - Vx[(jj-1) + ii*nxx]))/(107520.0f*dx);

            float dVz_dz = (75.0f*(Vz[jj + (ii-4)*nxx] - Vz[jj + (ii+3)*nxx]) +   
                          1029.0f*(Vz[jj + (ii+2)*nxx] - Vz[jj + (ii-3)*nxx]) +
                          8575.0f*(Vz[jj + (ii-2)*nxx] - Vz[jj + (ii+1)*nxx]) +
                        128625.0f*(Vz[jj + ii*nxx]     - Vz[jj + (ii-1)*nxx]))/(107520.0f*dz);     

            float L_int = powf(0.25f*(1.0f/L[jj + (ii+1)*nxx] + 1.0f/L[jj + ii*nxx] + 1.0f/L[(jj+1) + ii*nxx] + 1.0f/L[(jj+1) + (ii+1)*nxx]),-1.0f);
            float M_int = powf(0.25f*(1.0f/M[jj + (ii+1)*nxx] + 1.0f/M[jj + ii*nxx] + 1.0f/M[(jj+1) + ii*nxx] + 1.0f/M[(jj+1) + (ii+1)*nxx]),-1.0f);

            Txx[index] += dt*((L_int + 2.0f*M_int)*dVx_dx + L_int*dVz_dz);   
            Tzz[index] += dt*((L_int + 2.0f*M_int)*dVz_dz + L_int*dVx_dx);
        }

        if((ii >= 3) && (ii < nzz-4) && (jj >= 3) && (jj < nxx-4)) 
        { 
            float dVx_dz = (75*(Vx[jj + (ii-3)*nxx] - Vx[jj + (ii+4)*nxx]) +
                          1029*(Vx[jj + (ii+3)*nxx] - Vx[jj + (ii-2)*nxx]) +
                          8575*(Vx[jj + (ii-1)*nxx] - Vx[jj + (ii+2)*nxx]) +
                        128625*(Vx[jj + (ii+1)*nxx] - Vx[jj + ii*nxx]))/(107520.0f*dz);

            float dVz_dx = (75*(Vz[(jj-3) + ii*nxx] - Vz[(jj+4) + ii*nxx]) +
                          1029*(Vz[(jj+3) + ii*nxx] - Vz[(jj-2) + ii*nxx]) +
                          8575*(Vz[(jj-1) + ii*nxx] - Vz[(jj+2) + ii*nxx]) +
                        128625*(Vz[(jj+1) + ii*nxx] - Vz[jj + ii*nxx]))/(107520.0f*dx);

            Txz[index] += dt*M[index]*(dVx_dz + dVz_dx);            
        }      
    }

    # pragma acc parallel loop present(Vx[0:nxx*nzz],Vz[0:nxx*nzz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz],Txz[0:nxx*nzz],rho[0:nxx*nzz])
    for(int index = 0; index < nxx*nzz; index++) 
    {              
        int ii = (int) index / nxx;      /* Line indicator */
        int jj = (int) index % nxx;      /* Column indicator */ 

        if((ii > 3) && (ii < nzz-3) && (jj >= 3) && (jj < nxx-4)) 
        {
            float dTxx_dx = (75.0f*(Txx[(jj-3) + ii*nxx] - Txx[(jj+4) + ii*nxx]) +
                           1029.0f*(Txx[(jj+3) + ii*nxx] - Txx[(jj-2) + ii*nxx]) +
                           8575.0f*(Txx[(jj-1) + ii*nxx] - Txx[(jj+2) + ii*nxx]) +
                         128625.0f*(Txx[(jj+1) + ii*nxx] - Txx[jj + ii*nxx]))/(107520.0f*dx);

            float dTxz_dz = (75.0f*(Txz[jj + (ii-4)*nxx] - Txz[jj + (ii+3)*nxx]) +
                           1029.0f*(Txz[jj + (ii+2)*nxx] - Txz[jj + (ii-3)*nxx]) + 
                           8575.0f*(Txz[jj + (ii-2)*nxx] - Txz[jj + (ii+1)*nxx]) +
                         128625.0f*(Txz[jj + ii*nxx]     - Txz[jj + (ii-1)*nxx]))/(107520.0f*dz);

            float rho_int = 0.5f*(rho[jj + ii*nxx] + rho[(jj+1) + ii*nxx]);

            Vx[index] += dt/rho_int*(dTxx_dx + dTxz_dz);  
        }

        if((ii >= 3) && (ii < nzz-4) && (jj > 3) && (jj < nxx-3)) 
        {
            float dTxz_dx = (75.0f*(Txz[(jj-4) + ii*nxx] - Txz[(jj+3) + ii*nxx]) +
                           1029.0f*(Txz[(jj+2) + ii*nxx] - Txz[(jj-3) + ii*nxx]) +
                           8575.0f*(Txz[(jj-2) + ii*nxx] - Txz[(jj+1) + ii*nxx]) +
                         128625.0f*(Txz[jj + ii*nxx]     - Txz[(jj-1) + ii*nxx]))/(107520.0f*dx);

            float dTzz_dz = (75.0f*(Tzz[jj + (ii-3)*nxx] - Tzz[jj + (ii+4)*nxx]) + 
                           1029.0f*(Tzz[jj + (ii+3)*nxx] - Tzz[jj + (ii-2)*nxx]) +
                           8575.0f*(Tzz[jj + (ii-1)*nxx] - Tzz[jj + (ii+2)*nxx]) +
                         128625.0f*(Tzz[jj + (ii+1)*nxx] - Tzz[jj + ii*nxx]))/(107520.0f*dz);

            float rho_int = 0.5f*(rho[jj + ii*nxx] + rho[jj + (ii+1)*nxx]);

            Vz[index] += dt/rho_int*(dTxz_dx + dTzz_dz); 
        }
    }
}

void getAcousticPressureSeismogram(float *seism, float *P, int nt, int nxx, int nzz, int timePointer, int zrec)
{
    # pragma acc parallel loop present(seism[0:nxx*nt],P[0:nxx*nzz])
    for(int index = 0; index < nxx; index++) 
    {
        seism[timePointer*nxx + index] = P[zrec*nxx + index];
    }
}

void getElasticIsotropicPressureSeismogram(float *seism, float *Txx, float *Tzz, int nt, int nxx, int nzz, int timePointer, int zrec)
{
    # pragma acc parallel loop present(seism[0:nxx*nt],Txx[0:nxx*nzz],Tzz[0:nxx*nzz])
    for(int index = 0; index < nxx; index++) 
    {
        seism[timePointer*nxx + index] = (Txx[zrec*nxx + index] + Tzz[zrec*nxx + index])/2.0f;
    }
}

# endif
