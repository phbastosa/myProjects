# ifndef FUNCTIONS_H_DEFINED
# define FUNCTIONS_H_DEFINED

void memSet(float*pointer,int n)
{
    for(int i = 0; i < n; i++) pointer[i] = 0.0f;
}

void setWaveField(float*Vx,float*Vy,float*Vz,float*Txx,float*Tyy,float*Tzz,float*Txz,float*Tyz,float*Txy,int n)
{
    memSet(Vx,n);  memSet(Vy,n);  memSet(Vz,n);
    memSet(Txx,n); memSet(Tyy,n); memSet(Tzz,n);
    memSet(Txz,n); memSet(Tyz,n); memSet(Txy,n);
}

void importFloatVector(float*vector,int nPoints,char filename[])
{
    FILE * read = fopen((const char *) filename,"rb");
    fread(vector,sizeof(float),nPoints,read);
    fclose(read);
}

void importIntegerVector(int*vector,int nPoints,char filename[])
{
    FILE * read = fopen((const char *) filename,"rb");
    fread(vector,sizeof(int),nPoints,read);
    fclose(read);
}

void exportVector(float*vector,int nPoints,char filename[])
{
    FILE * write = fopen((const char *) filename, "wb");
    fwrite(vector, sizeof(float), nPoints, write);
    fclose(write);
}

void readParameters(int*nx,int*ny,int*nz,int*nt,float*dx,float*dy,float*dz,float*dt,int*abc,int*nrecx,int*nrecy,int*nsrc,int*wbh,char filename[])
{
    FILE * arq = fopen((const char *) filename,"r"); 
    if(arq != NULL) 
    {
        fscanf(arq,"%i",nx); fscanf(arq,"%i",ny); fscanf(arq,"%i",nz); fscanf(arq,"%i",nt); 
        fscanf(arq,"%f",dx); fscanf(arq,"%f",dy); fscanf(arq,"%f",dz); fscanf(arq,"%f",dt); 
        fscanf(arq,"%i",abc); fscanf(arq,"%i",nrecx); fscanf(arq,"%i",nrecy); 
        fscanf(arq,"%i",nsrc); fscanf(arq,"%i",wbh); 
    } 
    fclose(arq);
} 

void checkIndexes(int nxx,int nyy,int nzz)
{
    int contkk = 0; int contjj = 0; int contii = 0;
    for (int index = 0; index < nxx*nyy*nzz; index++)
    {
        int kk = floor(index/(nxx*nyy));           // indicador de matrizes (direção z)
        int jj = index % nxx;                      // indicador de colunas  (direção x)        
        int ii = floor((index % (nxx*nyy)) / nxx); // indicador de linhas   (direção y)  

        if ((jj == 0) && (ii == 0)) contkk += 1;
        if ((jj == 0) && (kk == 0)) contii += 1;
        if ((kk == 0) && (ii == 0)) contjj += 1;
    }
    printf("matrizes = %i; colunas = %i; linhas = %i\n",contkk,contjj,contii);
}

/* Equação discreta para propagação de ondas em meios elásticos isotrópicos 3D (Graves,1996) com matrizes indexadas. Autor: Paulo Bastos - GISIS/UFF, 22/07/2019 */
void FDM8E2T_elasticIsotropic3D(float*Vx,float*Vy,float*Vz,float*Txx,float*Tyy,float*Tzz,float*Txy,float*Txz,float*Tyz,float*rho,float*M,float*L,float*source,int*xsrc,int*ysrc,int*zsrc,int nsrc,int timePointer,int shotPointer,int nx,int ny,int nz,float dx,float dy,float dz,float dt) 
{
    #pragma omp parallel for
    for(int index = 0; index < nx*ny*nz; index++) 
    {    
        int kk = floor(index/(nx*ny));          // indicador de matrizes (direção z)
        int jj = index % nx;                    // indicador de colunas  (direção x)
        int ii = floor((index % (nx*ny)) / nx); // indicador de linhas   (direção y)  

        if((index == 0) && (timePointer < nsrc))
        {
            Txx[zsrc[shotPointer]*ny*nx + ysrc[shotPointer]*nx + xsrc[shotPointer]] += source[timePointer] / (dx*dy*dz);        
            Tyy[zsrc[shotPointer]*ny*nx + ysrc[shotPointer]*nx + xsrc[shotPointer]] += source[timePointer] / (dx*dy*dz);        
            Tzz[zsrc[shotPointer]*ny*nx + ysrc[shotPointer]*nx + xsrc[shotPointer]] += source[timePointer] / (dx*dy*dz);        
        }

        if((kk >= 3) && (kk < nz-4) && (ii >= 3) && (ii < ny-4) && (jj >= 3) && (jj < nx-4)) 
        {    
            float d_Vx_dx = (75.0f*(Vx[kk*ny*nx + ii*nx + (jj-3)] - Vx[kk*ny*nx + ii*nx + (jj+4)]) +
                           1029.0f*(Vx[kk*ny*nx + ii*nx + (jj+3)] - Vx[kk*ny*nx + ii*nx + (jj-2)]) +
                           8575.0f*(Vx[kk*ny*nx + ii*nx + (jj-1)] - Vx[kk*ny*nx + ii*nx + (jj+2)]) +
                         128625.0f*(Vx[kk*ny*nx + ii*nx + (jj+1)] - Vx[kk*ny*nx + ii*nx + jj]))/(107520.0f*dx);

            float d_Vy_dy = (75.0f*(Vy[kk*ny*nx + (ii-3)*nx + jj] - Vy[kk*ny*nx + (ii+4)*nx + jj]) +
                           1029.0f*(Vy[kk*ny*nx + (ii+3)*nx + jj] - Vy[kk*ny*nx + (ii-2)*nx + jj]) +
                           8575.0f*(Vy[kk*ny*nx + (ii-1)*nx + jj] - Vy[kk*ny*nx + (ii+2)*nx + jj]) +
                         128625.0f*(Vy[kk*ny*nx + (ii+1)*nx + jj] - Vy[kk*ny*nx + ii*nx + jj]))/(107520.0f*dy);

            float d_Vz_dz = (75.0f*(Vz[(kk-3)*ny*nx + ii*nx + jj] - Vz[(kk+4)*ny*nx + ii*nx + jj]) +
                           1029.0f*(Vz[(kk+3)*ny*nx + ii*nx + jj] - Vz[(kk-2)*ny*nx + ii*nx + jj]) +
                           8575.0f*(Vz[(kk-1)*ny*nx + ii*nx + jj] - Vz[(kk+2)*ny*nx + ii*nx + jj]) +
                         128625.0f*(Vz[(kk+1)*ny*nx + ii*nx + jj] - Vz[kk*ny*nx + ii*nx + jj]))/(107520.0f*dz);

            Txx[index] += dt*((L[index] + 2*M[index])*d_Vx_dx + L[index]*(d_Vy_dy + d_Vz_dz));
            Tyy[index] += dt*((L[index] + 2*M[index])*d_Vy_dy + L[index]*(d_Vx_dx + d_Vz_dz));
            Tzz[index] += dt*((L[index] + 2*M[index])*d_Vz_dz + L[index]*(d_Vx_dx + d_Vy_dy));                    
        }

        if((kk >= 0) && (kk < nz) && (ii > 3) && (ii < ny-3) && (jj > 3) && (jj < nx-3)) 
        {
            float d_Vx_dy = (75.0f*(Vx[kk*ny*nx + (ii-4)*nx + jj] - Vx[kk*ny*nx + (ii+3)*nx + jj]) +
                           1029.0f*(Vx[kk*ny*nx + (ii+2)*nx + jj] - Vx[kk*ny*nx + (ii-3)*nx + jj]) +
                           8575.0f*(Vx[kk*ny*nx + (ii-2)*nx + jj] - Vx[kk*ny*nx + (ii+1)*nx + jj]) +
                         128625.0f*(Vx[kk*ny*nx + ii*nx + jj]     - Vx[kk*ny*nx + (ii-1)*nx + jj]))/(107520.0f*dy);

            float d_Vy_dx = (75.0f*(Vy[kk*ny*nx + ii*nx + (jj-4)] - Vy[kk*ny*nx + ii*nx + (jj+3)]) +
                           1029.0f*(Vy[kk*ny*nx + ii*nx + (jj+2)] - Vy[kk*ny*nx + ii*nx + (jj-3)]) +
                           8575.0f*(Vy[kk*ny*nx + ii*nx + (jj-2)] - Vy[kk*ny*nx + ii*nx + (jj+1)]) +
                         128625.0f*(Vy[kk*ny*nx + ii*nx + jj]     - Vy[kk*ny*nx + ii*nx + (jj-1)]))/(107520.0f*dx);

            float M_xy = powf(0.25*(1/M[kk*ny*nx + (ii+1)*nx + (jj+1)] + 1/M[kk*ny*nx + ii*nx + (jj+1)] + 
                                    1/M[kk*ny*nx + (ii+1)*nx + jj]     + 1/M[kk*ny*nx + ii*nx + jj]),-1.0f);

            Txy[index] += dt*M_xy*(d_Vx_dy + d_Vy_dx);
        }

        if((kk > 3) && (kk < nz-3) && (ii >= 0) && (ii < ny) && (jj > 3) && (jj < nx-3)) 
        {
            float d_Vx_dz = (75.0f*(Vx[(kk-4)*ny*nx + ii*nx + jj] - Vx[(kk+3)*ny*nx + ii*nx + jj]) +
                           1029.0f*(Vx[(kk+2)*ny*nx + ii*nx + jj] - Vx[(kk-3)*ny*nx + ii*nx + jj]) +
                           8575.0f*(Vx[(kk-2)*ny*nx + ii*nx + jj] - Vx[(kk+1)*ny*nx + ii*nx + jj]) +
                         128625.0f*(Vx[kk*ny*nx + ii*nx + jj]     - Vx[(kk-1)*ny*nx + ii*nx + jj]))/(107520.0f*dz);

            float d_Vz_dx = (75.0f*(Vz[kk*ny*nx + ii*nx + (jj-4)] - Vz[kk*ny*nx + ii*nx + (jj+3)]) +
                           1029.0f*(Vz[kk*ny*nx + ii*nx + (jj+2)] - Vz[kk*ny*nx + ii*nx + (jj-3)]) +
                           8575.0f*(Vz[kk*ny*nx + ii*nx + (jj-2)] - Vz[kk*ny*nx + ii*nx + (jj+1)]) +
                         128625.0f*(Vz[kk*ny*nx + ii*nx + jj]     - Vz[kk*ny*nx + ii*nx + (jj-1)]))/(107520.0f*dx);

            float M_xz = powf(0.25*(1/M[(kk+1)*ny*nx + ii*nx + (jj+1)] + 1/M[kk*ny*nx + ii*nx + (jj+1)] + 
                                    1/M[(kk+1)*ny*nx + ii*nx + jj]     + 1/M[kk*ny*nx + ii*nx + jj]),-1.0f);

            Txz[index] += dt*M_xz*(d_Vx_dz + d_Vz_dx);
        }
    
        if((kk > 3) && (kk < nz-3) && (ii > 3) && (ii < ny-3) && (jj >= 0) && (jj < nx)) 
        {
            float d_Vy_dz = (75.0f*(Vy[(kk-4)*ny*nx + ii*nx + jj] - Vy[(kk+3)*ny*nx + ii*nx + jj]) +
                           1029.0f*(Vy[(kk+2)*ny*nx + ii*nx + jj] - Vy[(kk-3)*ny*nx + ii*nx + jj]) +
                           8575.0f*(Vy[(kk-2)*ny*nx + ii*nx + jj] - Vy[(kk+1)*ny*nx + ii*nx + jj]) +
                         128625.0f*(Vy[kk*ny*nx + ii*nx + jj]     - Vy[(kk-1)*ny*nx + ii*nx + jj]))/(107520.0f*dz);

            float d_Vz_dy = (75.0f*(Vz[kk*ny*nx + (ii-4)*nx + jj] - Vz[kk*ny*nx + (ii+3)*nx + jj]) +
                           1029.0f*(Vz[kk*ny*nx + (ii+2)*nx + jj] - Vz[kk*ny*nx + (ii-3)*nx + jj]) +
                           8575.0f*(Vz[kk*ny*nx + (ii-2)*nx + jj] - Vz[kk*ny*nx + (ii+1)*nx + jj]) +
                         128625.0f*(Vz[kk*ny*nx + ii*nx + jj]     - Vz[kk*ny*nx + (ii-1)*nx + jj]))/(107520.0f*dy);

            float M_yz = powf(0.25*(1/M[(kk+1)*ny*nx + (ii+1)*nx + jj] + 1/M[(kk+1)*ny*nx + ii*nx + jj] + 
                                    1/M[kk*ny*nx + (ii+1)*nx + jj] +     1/M[kk*ny*nx + ii*nx + jj]),-1.0f);

            Tyz[index] += dt*M_yz*(d_Vy_dz + d_Vz_dy);
        }
    }

    #pragma omp parallel for
    for(int index = 0; index < nx*ny*nz; index++) 
    {    
        int kk = floor(index/(nx*ny));          // indicador de matrizes (direção z)
        int jj = index % nx;                    // indicador de colunas  (direção x)
        int ii = floor((index % (nx*ny)) / nx); // indicador de linhas   (direção y)  

        if((kk >= 3) && (kk < nz-4) && (ii >= 3) && (ii < ny-4) && (jj > 3) && (jj < nx-3)) 
        {
            float d_Txx_dx = (75.0f*(Txx[kk*ny*nx + ii*nx + (jj-4)] - Txx[kk*ny*nx + ii*nx + (jj+3)]) +
                            1029.0f*(Txx[kk*ny*nx + ii*nx + (jj+2)] - Txx[kk*ny*nx + ii*nx + (jj-3)]) +
                            8575.0f*(Txx[kk*ny*nx + ii*nx + (jj-2)] - Txx[kk*ny*nx + ii*nx + (jj+1)]) +
                          128625.0f*(Txx[kk*ny*nx + ii*nx + jj]     - Txx[kk*ny*nx + ii*nx + (jj-1)]))/(107520.0f*dx);

            float d_Txy_dy = (75.0f*(Txy[kk*ny*nx + (ii-3)*nx + jj] - Txy[kk*ny*nx + (ii+4)*nx + jj]) +
                            1029.0f*(Txy[kk*ny*nx + (ii+3)*nx + jj] - Txy[kk*ny*nx + (ii-2)*nx + jj]) +
                            8575.0f*(Txy[kk*ny*nx + (ii-1)*nx + jj] - Txy[kk*ny*nx + (ii+2)*nx + jj]) +
                          128625.0f*(Txy[kk*ny*nx + (ii+1)*nx + jj] - Txy[kk*ny*nx + ii*nx + jj]))/(107520.0f*dy);

            float d_Txz_dz = (75.0f*(Txz[(kk-3)*ny*nx + ii*nx + jj] - Txz[(kk+4)*ny*nx + ii*nx + jj]) +
                            1029.0f*(Txz[(kk+3)*ny*nx + ii*nx + jj] - Txz[(kk-2)*ny*nx + ii*nx + jj]) +
                            8575.0f*(Txz[(kk-1)*ny*nx + ii*nx + jj] - Txz[(kk+2)*ny*nx + ii*nx + jj]) +
                          128625.0f*(Txz[(kk+1)*ny*nx + ii*nx + jj] - Txz[kk*ny*nx + ii*nx + jj]))/(107520.0f*dz);

            float rhox = 0.5f*(rho[kk*ny*nx + ii*nx + (jj+1)] + rho[kk*ny*nx + ii*nx + jj]);

            Vx[index] += dt/rhox*(d_Txx_dx + d_Txy_dy + d_Txz_dz); 
        }
    
        if((kk >= 3) && (kk < nz-3) && (ii > 3) && (ii < ny-3) && (jj >= 3) && (jj < nx-4)) 
        {
            float d_Txy_dx = (75.0f*(Txy[kk*ny*nx + ii*nx + (jj-3)] - Txy[kk*ny*nx + ii*nx + (jj+4)]) +
                            1029.0f*(Txy[kk*ny*nx + ii*nx + (jj+3)] - Txy[kk*ny*nx + ii*nx + (jj-2)]) +
                            8575.0f*(Txy[kk*ny*nx + ii*nx + (jj-1)] - Txy[kk*ny*nx + ii*nx + (jj+2)]) +
                          128625.0f*(Txy[kk*ny*nx + ii*nx + (jj+1)] - Txy[kk*ny*nx + ii*nx + jj]))/(107520.0f*dx);

            float  d_Tyy_dy = (75.0f*(Tyy[kk*ny*nx + (ii-4)*nx + jj] - Tyy[kk*ny*nx + (ii+3)*nx + jj]) +
                             1029.0f*(Tyy[kk*ny*nx + (ii+2)*nx + jj] - Tyy[kk*ny*nx + (ii-3)*nx + jj]) +
                             8575.0f*(Tyy[kk*ny*nx + (ii-2)*nx + jj] - Tyy[kk*ny*nx + (ii+1)*nx + jj]) +
                           128625.0f*(Tyy[kk*ny*nx + ii*nx + jj]     - Tyy[kk*ny*nx + (ii-1)*nx + jj]))/(107520.0f*dy);

            float d_Tyz_dz = (75.0f*(Tyz[(kk-3)*ny*nx + ii*nx + jj] - Tyz[(kk+4)*nz*nx + ii*nx + jj]) +
                            1029.0f*(Tyz[(kk+3)*ny*nx + ii*nx + jj] - Tyz[(kk-2)*ny*nx + ii*nx + jj]) +
                            8575.0f*(Tyz[(kk-1)*ny*nx + ii*nx + jj] - Tyz[(kk+2)*ny*nx + ii*nx + jj]) +
                          128625.0f*(Tyz[(kk+1)*ny*nx + ii*nx + jj] - Tyz[kk*ny*nx + ii*nx + jj]))/(107520.0f*dz);

            float rhoy = 0.5f*(rho[kk*ny*nx + (ii+1)*nx + jj] + rho[kk*ny*nx + ii*nx + jj]);

            Vy[index] += dt/rhoy*(d_Txy_dx + d_Tyy_dy + d_Tyz_dz); 
        }    

        if((kk > 3) && (kk < nz-3) && (ii >= 3) && (ii < ny-4) && (jj >= 3) && (jj < nx-4)) 
        {
            float d_Txz_dx = (75.0f*(Txz[kk*ny*nx + ii*nx + (jj-3)] - Txz[kk*ny*nx + ii*nx + (jj+4)]) +
                            1029.0f*(Txz[kk*ny*nx + ii*nx + (jj+3)] - Txz[kk*ny*nx + ii*nx + (jj-2)]) +
                            8575.0f*(Txz[kk*ny*nx + ii*nx + (jj-1)] - Txz[kk*ny*nx + ii*nx + (jj+2)]) +
                          128625.0f*(Txz[kk*ny*nx + ii*nx + (jj+1)] - Txz[kk*ny*nx + ii*nx + jj]))/(107520.0f*dx);

            float d_Tyz_dy = (75.0f*(Tyz[kk*ny*nx + (ii-3)*nx + jj] - Tyz[kk*ny*nx + (ii+4)*nx + jj]) +
                            1029.0f*(Tyz[kk*ny*nx + (ii+3)*nx + jj] - Tyz[kk*ny*nx + (ii-2)*nx + jj]) +
                            8575.0f*(Tyz[kk*ny*nx + (ii-1)*nx + jj] - Tyz[kk*ny*nx + (ii+2)*nx + jj]) +
                          128625.0f*(Tyz[kk*ny*nx + (ii+1)*nx + jj] - Tyz[kk*ny*nx + ii*nx + jj]))/(107520.0f*dy);

            float d_Tzz_dz = (75.0f*(Tzz[(kk-4)*ny*nx + ii*nx + jj] - Tzz[(kk+3)*ny*nx + ii*nx + jj]) +
                            1029.0f*(Tzz[(kk+2)*ny*nx + ii*nx + jj] - Tzz[(kk-3)*ny*nx + ii*nx + jj]) +
                            8575.0f*(Tzz[(kk-2)*ny*nx + ii*nx + jj] - Tzz[(kk+1)*ny*nx + ii*nx + jj]) +
                          128625.0f*(Tzz[kk*ny*nx + ii*nx + jj]     - Tzz[(kk-1)*ny*nx + ii*nx + jj]))/(107520.0f*dz);

            float rhoz = 0.5f*(rho[(kk+1)*ny*nx + ii*nx + jj] + rho[kk*ny*nx + ii*nx + jj]);

            Vz[index] += dt/rhoz*(d_Txz_dx + d_Tyz_dy + d_Tzz_dz); 
        }
    }    
}

void cerjanElasticAbsorbingCondition3D(float*Vx,float*Vy,float*Vz,float*Txx,float*Tyy,float*Tzz,float*Txy,float*Txz,float*Tyz,float*cubXYZ,int n)
{
    #pragma omp parallel for
    for(int index = 0; index < n; index++)    
    {
        Vx[index] *= cubXYZ[index];
        Vy[index] *= cubXYZ[index];
        Vz[index] *= cubXYZ[index];
        Txx[index] *= cubXYZ[index];
        Txy[index] *= cubXYZ[index];
        Txz[index] *= cubXYZ[index];
        Tyy[index] *= cubXYZ[index];
        Tyz[index] *= cubXYZ[index];
        Tzz[index] *= cubXYZ[index];
    }
}

void getPressureWaveField(float*Txx,float*Tyy,float*Tzz,float*Pss,int nPoints)
{
    #pragma omp parallel for
    for (int index = 0; index < nPoints; index++)
    {
        Pss[index] = (Txx[index] + Tyy[index] + Tzz[index])/3.0f;
    }
}

void getPWaveField(float*Ux,float*Uy,float*Uz,float*P,int nx,int ny,int nz,float dx,float dy,float dz)
{
    #pragma omp parallel for
    for(int index = 0; index < nx*ny*nz; index++) 
    {    
        int kk = floor(index/(nx*ny));          // indicador de matrizes (direção z)
        int jj = index % nx;                    // indicador de colunas  (direção x)
        int ii = floor((index % (nx*ny)) / nx); // indicador de linhas   (direção y)  

        if((kk >= 3) && (kk < nz-4) && (ii >= 3) && (ii < ny-4) && (jj >= 3) && (jj < nx-4)) 
        {    
            float d_Ux_dx = (75.0f*(Ux[kk*ny*nx + ii*nx + (jj-3)] - Ux[kk*ny*nx + ii*nx + (jj+4)]) +
                           1029.0f*(Ux[kk*ny*nx + ii*nx + (jj+3)] - Ux[kk*ny*nx + ii*nx + (jj-2)]) +
                           8575.0f*(Ux[kk*ny*nx + ii*nx + (jj-1)] - Ux[kk*ny*nx + ii*nx + (jj+2)]) +
                         128625.0f*(Ux[kk*ny*nx + ii*nx + (jj+1)] - Ux[kk*ny*nx + ii*nx + jj]))/(107520.0f*dx);

            float d_Uy_dy = (75.0f*(Uy[kk*ny*nx + (ii-3)*nx + jj] - Uy[kk*ny*nx + (ii+4)*nx + jj]) +
                           1029.0f*(Uy[kk*ny*nx + (ii+3)*nx + jj] - Uy[kk*ny*nx + (ii-2)*nx + jj]) +
                           8575.0f*(Uy[kk*ny*nx + (ii-1)*nx + jj] - Uy[kk*ny*nx + (ii+2)*nx + jj]) +
                         128625.0f*(Uy[kk*ny*nx + (ii+1)*nx + jj] - Uy[kk*ny*nx + ii*nx + jj]))/(107520.0f*dy);

            float d_Uz_dz = (75.0f*(Uz[(kk-3)*ny*nx + ii*nx + jj] - Uz[(kk+4)*ny*nx + ii*nx + jj]) +
                           1029.0f*(Uz[(kk+3)*ny*nx + ii*nx + jj] - Uz[(kk-2)*ny*nx + ii*nx + jj]) +
                           8575.0f*(Uz[(kk-1)*ny*nx + ii*nx + jj] - Uz[(kk+2)*ny*nx + ii*nx + jj]) +
                         128625.0f*(Uz[(kk+1)*ny*nx + ii*nx + jj] - Uz[kk*ny*nx + ii*nx + jj]))/(107520.0f*dz);
     
            P[index] = d_Ux_dx + d_Uy_dy + d_Uz_dz;    
        }
    }
}

void getSWaveField(float*Ux,float*Uy,float*Uz,float*Shx,float*Shy,float*Sv,int nx,int ny,int nz,float dx,float dy,float dz)
{
    #pragma omp parallel for
    for(int index = 0; index < nx*ny*nz; index++) 
    {    
        int kk = floor(index/(nx*ny));          // indicador de matrizes (direção z)
        int jj = index % nx;                    // indicador de colunas  (direção x)
        int ii = floor((index % (nx*ny)) / nx); // indicador de linhas   (direção y)  

        if((kk > 3) && (kk < nz-3) && (ii > 3) && (ii < ny-3) && (jj >= 0) && (jj < nx)) 
        {
            float d_Uz_dy = (75.0f*(Uz[kk*ny*nx + (ii-4)*nx + jj] - Uz[kk*ny*nx + (ii+3)*nx + jj]) +
                           1029.0f*(Uz[kk*ny*nx + (ii+2)*nx + jj] - Uz[kk*ny*nx + (ii-3)*nx + jj]) +
                           8575.0f*(Uz[kk*ny*nx + (ii-2)*nx + jj] - Uz[kk*ny*nx + (ii+1)*nx + jj]) +
                         128625.0f*(Uz[kk*ny*nx + ii*nx + jj]     - Uz[kk*ny*nx + (ii-1)*nx + jj]))/(107520.0f*dy);
 
            float d_Uy_dz = (75.0f*(Uy[(kk-4)*ny*nx + ii*nx + jj] - Uy[(kk+3)*ny*nx + ii*nx + jj]) +
                           1029.0f*(Uy[(kk+2)*ny*nx + ii*nx + jj] - Uy[(kk-3)*ny*nx + ii*nx + jj]) +
                           8575.0f*(Uy[(kk-2)*ny*nx + ii*nx + jj] - Uy[(kk+1)*ny*nx + ii*nx + jj]) +
                         128625.0f*(Uy[kk*ny*nx + ii*nx + jj]     - Uy[(kk-1)*ny*nx + ii*nx + jj]))/(107520.0f*dz);

            Shx[index] = d_Uz_dy - d_Uy_dz;
        }

        if((kk > 3) && (kk < nz-3) && (ii >= 0) && (ii < ny) && (jj > 3) && (jj < nx-3)) 
        {
            float d_Ux_dz = (75.0f*(Ux[(kk-4)*ny*nx + ii*nx + jj] - Ux[(kk+3)*ny*nx + ii*nx + jj]) +
                           1029.0f*(Ux[(kk+2)*ny*nx + ii*nx + jj] - Ux[(kk-3)*ny*nx + ii*nx + jj]) +
                           8575.0f*(Ux[(kk-2)*ny*nx + ii*nx + jj] - Ux[(kk+1)*ny*nx + ii*nx + jj]) +
                         128625.0f*(Ux[kk*ny*nx + ii*nx + jj]     - Ux[(kk-1)*ny*nx + ii*nx + jj]))/(107520.0f*dz);

            float d_Uz_dx = (75.0f*(Uz[kk*ny*nx + ii*nx + (jj-4)] - Uz[kk*ny*nx + ii*nx + (jj+3)]) +
                           1029.0f*(Uz[kk*ny*nx + ii*nx + (jj+2)] - Uz[kk*ny*nx + ii*nx + (jj-3)]) +
                           8575.0f*(Uz[kk*ny*nx + ii*nx + (jj-2)] - Uz[kk*ny*nx + ii*nx + (jj+1)]) +
                         128625.0f*(Uz[kk*ny*nx + ii*nx + jj]     - Uz[kk*ny*nx + ii*nx + (jj-1)]))/(107520.0f*dx);

            Shy[index] = d_Ux_dz - d_Uz_dx;
        }
    
        if((kk >= 0) && (kk < nz) && (ii > 3) && (ii < ny-3) && (jj > 3) && (jj < nx-3)) 
        {
            float d_Uy_dx = (75.0f*(Uy[kk*ny*nx + ii*nx + (jj-4)] - Uy[kk*ny*nx + ii*nx + (jj+3)]) +
                           1029.0f*(Uy[kk*ny*nx + ii*nx + (jj+2)] - Uy[kk*ny*nx + ii*nx + (jj-3)]) +
                           8575.0f*(Uy[kk*ny*nx + ii*nx + (jj-2)] - Uy[kk*ny*nx + ii*nx + (jj+1)]) +
                         128625.0f*(Uy[kk*ny*nx + ii*nx + jj]     - Uy[kk*ny*nx + ii*nx + (jj-1)]))/(107520.0f*dx);

            float d_Ux_dy = (75.0f*(Ux[kk*ny*nx + (ii-4)*nx + jj] - Ux[kk*ny*nx + (ii+3)*nx + jj]) +
                           1029.0f*(Ux[kk*ny*nx + (ii+2)*nx + jj] - Ux[kk*ny*nx + (ii-3)*nx + jj]) +
                           8575.0f*(Ux[kk*ny*nx + (ii-2)*nx + jj] - Ux[kk*ny*nx + (ii+1)*nx + jj]) +
                         128625.0f*(Ux[kk*ny*nx + ii*nx + jj]     - Ux[kk*ny*nx + (ii-1)*nx + jj]))/(107520.0f*dy);

            Sv[index] = d_Uy_dx - d_Ux_dy;
        }
    }
}

void getSeismogram(float*seism,float*field,int*xrec,int*yrec,int*zrec,int nrecx,int nrecy,int nrecs,int nt,int nxx,int nyy,int nzz,int timePointer,int nsrc,float dt)
{
    if(timePointer > nsrc/2)
    {
        #pragma omp parallel for
        for (int index = 0; index < nrecs; index++)
        {
            int recx = index % nrecx;                          
            int recy = floor((index % (nrecx*nrecy)) / nrecx);   

            seism[timePointer*nrecx*nrecy + recy*nrecx + recx] = (timePointer - nsrc/2)*dt*field[zrec[index]*nxx*nyy + yrec[index]*nxx + xrec[index]];
        }
    } 
}

# endif
