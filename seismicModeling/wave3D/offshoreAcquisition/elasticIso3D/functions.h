# ifndef FUNCTIONS_H_DEFINED
# define FUNCTIONS_H_DEFINED

//
float getMax(float * array, int n)
{
    float max = 0.0f;

    for (int i = 0; i < n; i++)
    {
        if (array[i] > max)
        {
            max = array[i];
        }
    }

    return max;
}

//
void memSet(float * pointer, int n)
{
    for(int i = 0; i < n; i++) pointer[i] = 0.0f;
}

//
float * transpose(float * matrix, int rows, int columns)
{
    int ii, jj;
    float * trasp = (float *) malloc(rows*columns*sizeof(float));

    for(ii = 0; ii < rows; ii++)
    {
        for(jj = 0; jj < columns; jj++)
        {
            trasp[ii*columns + jj] = matrix[jj*rows + ii];
        }
    }
    return trasp;
}

float* import_float32(char* file, int n_points) 
{
    FILE * fp;
	size_t result;

	fp = fopen(file, "rb");
	if (fp == NULL) {fputs ("File error\n", stderr); exit (1);}

	float* buffer = (float*) malloc(n_points*sizeof(float));

	if (buffer == NULL) {fputs ("Memory error\n",stderr); exit (2);}

	result = fread(buffer,sizeof(float),n_points,fp);
	if (result != n_points) {fputs ("Reading error\n",stderr); exit (3);}

	fclose (fp);
	return (buffer);
}

//
void importVector(float * vector, int nPoints, char filename[])
{
    FILE * read = fopen((const char *) filename,"rb");
    fread(vector,sizeof(float),nPoints,read);
    // fclose(read);
}

//
void exportVector(float * vector, int nPoints, char filename[])
{
    FILE * write = fopen((const char *) filename, "wb");
    fwrite(vector, sizeof(float),nPoints, write);
    // fclose(write);
}

//
void readIntegerTextParameter(int * vector, char filename[])
{
    FILE * rf = fopen((char *) filename, "r");
    if(rf != NULL) {
        int ii = 0;
        while(fscanf(rf,"%i",&vector[ii]) != EOF) {
            ++ii;                    
        }
    } 
    fclose(rf);
}    

//
void readFloatTextParameter(float * vector, char filename[])
{
    FILE * rf = fopen((char *) filename, "r");
    if(rf != NULL) 
    {
        int ii = 0;
        while(fscanf(rf,"%f",&vector[ii]) != EOF) 
        {
            ++ii;                    
        }
    } 
    fclose(rf);
}   

//
void readParameters(int * nxx,int * nzz,int * nt,float * dx,float * dz, float * dt,int * abc_layer,float * fac_abc,float * freq, int * spread,int * n_shot,int * nsrc,int * zsrc,char filename[])
{
    FILE * arq = fopen(filename,"r");
    fscanf(arq,"%i %i %i %f %f %f %i %f %f %i %i %i %i\n",nxx,nzz,nt,dx,dz,dt,abc_layer,fac_abc,freq,spread,n_shot,nsrc,zsrc);
    fclose(arq);
}

//
float * generateRicker(float f_cut, float dt, int nsrc)
{
    float pi = 4.0f*atan(1.0f);
    float fc = f_cut / (3 * sqrt(pi));
    float aux1, aux2, aux3;

    float * ricker = (float *) malloc(nsrc*sizeof(float));

    for(int ii = -nsrc/2; ii < nsrc/2; ii++) 
    {     
    	aux1 = (1.0f - 2.0f*pi*pow(ii*dt,2.0f)*pow(fc,2.0f)*pow(pi,2.0f)); 
        aux2 = exp(-pi*pow(ii*dt,2.0)*pow(fc,2.0)*pow(pi,2.0));
        ricker[ii + nsrc/2] = aux1 * aux2;
    }                                                              

    return ricker;
}

//
float * generateStaggeredSource(float f_cut, float dt, int nsrc)
{
    // float pi = 4.0f * atan(1.0f);
    // float fc = f_cut / (3 * sqrt(pi));
    // float aux1, aux2;

    
    float * ricker = (float *) malloc(nsrc*sizeof(float));
    float * source = (float *) malloc(nsrc*sizeof(float));
    
    // for (int ii = -nsrc/2; ii < nsrc/2; ii++)
    // {
    //     aux1 = ii*dt*pi*fc;
    //     aux2 = exp(-pi*powf(ii*dt,2.0f)*powf(fc,2.0f)*powf(pi,2.0f)); 

    //     source[ii + nsrc/2] = aux1 * aux2;
    // }

    ricker = generateRicker(f_cut,dt,nsrc);

    for (int ii = 0; ii < nsrc; ii++)
    {
        for (int jj = 0; jj < ii+1; jj++)
        {    
            source[ii] += ricker[jj];
        }
    }

    return source;
}

//
void printStatus(float * vels, float dx, float dz, float dt, float freq, int nxx, int nzz, int abc, int nt)
{
    printf("Reverse time migration\n\n");
    printf("Parameters:\n");
    printf("   Highest velocity in the model = %.1f m/s\n",vels[0]);
    printf("   Slowest velocity in the model = %.1f m/s\n",vels[1]);
    printf("   Spatial X discratization = %.1f meters\n",dx);
    printf("   Spatial Z discratization = %.1f meters\n",dz);
    printf("   Temporal discratization = %.4f seconds\n",dt);
    printf("   Source cutoff frequency = %.1f Hz\n",freq);
    printf("   Total modeling time = %.2f seconds\n",dt*nt);
    printf("   Horizontal length of model = %.0f meters\n",dx*(nxx - 2*abc));
    printf("   Vertical length of model = %.0f meters\n\n",dz*(nzz - 2*abc));
}

//
float * getVelocities2D(int nxx, int nzz, float * vp)
{
    int ii, jj, index;
    float v_max = vp[0];
    float v_min = vp[0];
    float * vels = (float *) malloc(2*sizeof(float));

    for(index = 0; index < nxx*nzz; ++index) 
    {
        ii = floor(index / nxx);  /* line indicator */
        jj = index % nxx;         /* column indicator */  

        /* Finding the highest velocity in the model */
        if(v_max < vp[ii*nxx + jj]) v_max = vp[ii*nxx + jj]; 
        
        /* Finding the slowest velocity in the model */
        if(v_min > vp[ii*nxx + jj]) v_min = vp[ii*nxx + jj]; 
    }

    vels[0] = v_max;
    vels[1] = v_min;    
 
    return vels;
}

//
void shotStatus2D(int tt, int * xsrc, int n_shot, int * xrec, int spread, float dx, float dz, int nt, float * vels, float dt, float freq, int nxx, int nzz, int abc)
{
        system("clear");        
        printStatus(vels,dx,dz,dt,freq,nxx,nzz,abc,nt);    
        printf("Migration status:\n");
        printf("   Shot position: %.1f meters\n",xsrc[tt]*dx);
        printf("   Recivers position: %.1f - %.1f meters\n",xrec[tt*spread]*dx,xrec[spread-1 + tt*spread]*dx);
        // printf("   Shot progress %i: %.2f %%",);
        printf("   Total progress: %.2f %%\n",(float) tt/n_shot * 100.0f);
        printf("\nImages cross-correlated: %i of %i\n",tt,n_shot);
}   

//
float * factorAttenuation(int borda, float parameter) 
{
    float * factor = (float *) malloc(borda*sizeof(float));

    for(int index = 0; index < borda; index++) {
        factor[index] = exp(-powf(parameter*(borda-index),2.0f)); 
    }

    return factor;
}

//
void cerjanAcousticAbsorbingCondition1D(float * P_pre, float * P_fut, int nx, float * factor, int borda)
{
    for (int index = 0; index < borda; index++)
    {
        P_pre[index] = P_pre[index] * factor[index];
        P_fut[index] = P_fut[index] * factor[index];

        P_pre[nx-index] = P_pre[nx-index] * factor[index];
        P_fut[nx-index] = P_fut[nx-index] * factor[index];
    }    
}

//
float * cerjanCorners2D(float * factor, int n_abc) 
{
    int ii,kk;
    float * corner = (float *) malloc(n_abc*n_abc*sizeof(float));

    /* Upper left absorbing corner */
    for(ii = n_abc-1; ii >= 0; ii--) 
    {
        corner[ii*n_abc + ii] = factor[ii];    
        for(kk = 1; kk <= n_abc; kk++) 
        {
            corner[(ii+kk)*n_abc + ii] = factor[ii];
            corner[ii*n_abc + ii + kk] = factor[ii];
        }    
    }
    return corner;
}

//
void cerjanAcousticAbsorbingCondition2D(float *P_pre, float *P_fut, float *factor, int nxx, int nzz, int abc_layer, float * corner)
{
    int ii,jj,index;
    #pragma acc parallel loop present(P_pre[0:nxx*nzz],P_fut[0:nxx*nzz],factor[0:abc_layer],corner[0:abc_layer*abc_layer]) 
    for(index = 0; index < nxx*nzz; ++index) 
    {     
        ii = floor(index / nxx);  /* Line indicator */
        jj = index % nxx;         /* Column indicator */

        /* Left and right absorbing layers */ 
        if((ii >= abc_layer) && (ii < nzz-abc_layer) && (jj >= 0) && (jj < abc_layer)) 
        {
            P_pre[ii*nxx+jj] = P_pre[ii*nxx+jj]*factor[jj];            /* Left  */
            P_pre[ii*nxx+(nxx-jj)] = P_pre[ii*nxx+(nxx-jj)]*factor[jj];/* Right */

            P_fut[ii*nxx+jj] = P_fut[ii*nxx+jj]*factor[jj];            /* Left  */     
            P_fut[ii*nxx+(nxx-jj)] = P_fut[ii*nxx+(nxx-jj)]*factor[jj];/* Right */
        }

        /* Up and bottom absorbing layers */
        if((jj >= abc_layer) && (jj < nxx - abc_layer) && (ii >= 0) && (ii < abc_layer)) 
        {
            P_pre[ii*nxx+jj] = P_pre[ii*nxx+jj]*factor[ii];             /*   Up   */
            P_pre[(nzz-ii)*nxx-jj] = P_pre[(nzz-ii)*nxx-jj]*factor[ii]; /* Bottom */ 

            P_fut[ii*nxx+jj] = P_fut[ii*nxx+jj]*factor[ii];             /*   Up   */
            P_fut[(nzz-ii)*nxx-jj] = P_fut[(nzz-ii)*nxx-jj]*factor[ii]; /* Bottom */
        }
        
        /* Upper left absorbing corner */ 
        if((ii >= 0) && (ii < abc_layer) && (jj >= 0) && (jj < abc_layer)) 
        {
            P_pre[ii*nxx + jj] = P_pre[ii*nxx + jj]*corner[ii*abc_layer + jj];
            P_fut[ii*nxx + jj] = P_fut[ii*nxx + jj]*corner[ii*abc_layer + jj];
        }

        /* Upper right absorbing corner */
        if((ii >= 0) && (ii < abc_layer) && (jj >= nxx-abc_layer) && (jj < nxx)) 
        {
            P_pre[ii*nxx + jj] = P_pre[ii*nxx + jj]*corner[ii*abc_layer + (nxx-jj-1)];
            P_fut[ii*nxx + jj] = P_fut[ii*nxx + jj]*corner[ii*abc_layer + (nxx-jj-1)];
        }

        /* Bottom left absorbing corner */
        if((ii >= nzz-abc_layer) && (ii < nzz) && (jj >= 0) && (jj < abc_layer)) 
        {
            P_pre[ii*nxx + jj] = P_pre[ii*nxx + jj]*corner[(nzz-ii-1)*abc_layer + jj];
            P_fut[ii*nxx + jj] = P_fut[ii*nxx + jj]*corner[(nzz-ii-1)*abc_layer + jj];
        }

        /* Botton right absorbing corner */ 
        if((ii >= nzz-abc_layer) && (ii < nzz) && (jj >= nxx-abc_layer) && (jj < nxx)) 
        {
            P_pre[ii*nxx + jj] = P_pre[ii*nxx + jj]*corner[(nzz-ii-1)*abc_layer + (nxx-jj-1)];
            P_fut[ii*nxx + jj] = P_fut[ii*nxx + jj]*corner[(nzz-ii-1)*abc_layer + (nxx-jj-1)];
        }
    }
}

//
void shiftingCoordinates2D(int * xsrc, int * xrec, int nshot, int nrcpt, int absLayer)
{
    for (int ii = 0; ii < nrcpt; ii++)
        xrec[ii] += absLayer;

    for (int ii = 0; ii < nshot; ii++)
        xsrc[ii] += absLayer;
}

//
void FDM_8E2T_acoustic1D(int time_pointer,float * P_pas,float * P_pre,float * P_fut,float * vel,int nx,float dx,float dt,float * source,int nsrc,int sx)
{
    float aux1, aux2, aux3;

    for (int ii = 4; ii < nx-4; ii++)
    {
        if((time_pointer < nsrc) && (ii == 4))
        {
            P_pre[sx] = source[time_pointer];
        }

        aux1 = (- 9.0f*(P_pre[ii-4] + P_pre[ii+4])
            +   128.0f*(P_pre[ii-3] + P_pre[ii+3])
            -  1008.0f*(P_pre[ii-2] + P_pre[ii+2])
            +  8064.0f*(P_pre[ii-1] + P_pre[ii+1])
            - 14350.0f*(P_pre[ii]))/(5040.0f*powf(dx,2));

        aux2 = pow(dt,2.0)*pow(vel[ii],2.0);
        aux3 = 2.0*P_pre[ii] - P_pas[ii];

        P_fut[ii] = aux1 * aux2 + aux3;
    }
}

//
void FDM_8E2T_acoustic2D(int shot_pointer,int time_pointer,float *vp,float *P_pas,float *P_pre,float *P_fut,float *source,int nsrc,int * z_src,int *x_src,int n_shots,int nxx,int nzz,float dx,float dz,float dt)
{
    int ii,jj,index;
    float d2_Px2, d2_Pz2;

    #pragma acc parallel loop present(P_pre[0:nxx*nzz],P_pas[0:nxx*nzz],P_fut[0:nxx*nzz],vp[0:nxx*nzz],x_src[0:n_shots],source[0:nsrc])
    for(index = 0; index < nxx*nzz; ++index) /* Spatial loop */ 
    {
        ii = floor(index / nxx);  /* Line indicator */
        jj = index % nxx;         /* Column indicator */  

        if((time_pointer < nsrc) && (index == 0))
        { 
            /* Applying source */
            P_pre[z_src[shot_pointer]*nxx + x_src[shot_pointer]] += source[time_pointer] / (dx*dz); 
        }

        if((ii > 3) && (ii < nzz-4) && (jj > 3) && (jj < nxx-4)) 
        {
            /* Second derivative of the pressure with respect to x, 
            using eighth-order finite difference discretization */
            d2_Px2 = (- 9.0f*(P_pre[ii*nxx + (jj-4)] + P_pre[ii*nxx + (jj+4)])
                  +   128.0f*(P_pre[ii*nxx + (jj-3)] + P_pre[ii*nxx + (jj+3)])
                  -  1008.0f*(P_pre[ii*nxx + (jj-2)] + P_pre[ii*nxx + (jj+2)])
                  +  8064.0f*(P_pre[ii*nxx + (jj-1)] + P_pre[ii*nxx + (jj+1)])
                  - 14350.0f*(P_pre[ii*nxx + jj]))/(5040.0f*powf(dx,2));

            /* Second derivative of the pressure with respect to z,marmousi_migrado_141x383.bin 
            using eighth-order finite difference discretization */    
            d2_Pz2 = (- 9.0f*(P_pre[(ii-4)*nxx + jj] + P_pre[(ii+4)*nxx + jj])
                  +   128.0f*(P_pre[(ii-3)*nxx + jj] + P_pre[(ii+3)*nxx + jj])
                  -  1008.0f*(P_pre[(ii-2)*nxx + jj] + P_pre[(ii+2)*nxx + jj])
                  +  8064.0f*(P_pre[(ii-1)*nxx + jj] + P_pre[(ii+1)*nxx + jj])
                  - 14350.0f*(P_pre[ii*nxx + jj]))/(5040.0f*powf(dz,2));

            /* Calculating the future wave field through the previous fields */
            P_fut[ii*nxx + jj] = (powf(dt,2)*powf(vp[ii*nxx + jj],2)
                            * (d2_Px2 + d2_Pz2)) + 2.0f*P_pre[ii*nxx + jj] 
                            -  P_pas[ii*nxx + jj];
        }
    }
}

//
void FDM_8E2T_acoustic2D_depropagation(int shot_pointer, int time_pointer, float *vp, float *P_pre, float *P_pas, float *P_fut, float *seism, int z_src, int nt, int spread, int *x_src, int *xrec, int nxx, int nzz, float dx, float dz, float dt)
{
    int ii,jj,index;
    float d2_Px2, d2_Pz2;
    
    #pragma acc parallel loop present(seism[0:spread*nt],P_pre[0:nxx*nzz],xrec[0:spread*nt])
    for(index = 0; index <spread; index++)
    {
        P_pre[z_src*nxx + xrec[shot_pointer*spread + index]] = seism[time_pointer*spread + index]; /* Applying source*/
    }

    #pragma acc parallel loop present(P_pre[0:nxx*nzz],P_pas[0:nxx*nzz],P_fut[0:nxx*nzz],vp[0:nxx*nzz])
    for(index = 0; index < nxx*nzz; ++index) /* Spatial loop */ 
    {
        ii = floor(index / nxx);  /* Line indicator */
        jj = index % nxx;         /* Column indicator */  

        if((ii > 3) && (ii < nzz-4) && (jj > 3) && (jj < nxx-4)) 
        {
            /* Second derivative of the pressure with respect to x, 
            using eighth-order finite difference discretization */
            d2_Px2 = (- 9.0f*(P_pre[ii*nxx + (jj-4)] + P_pre[ii*nxx + (jj+4)])
                  +   128.0f*(P_pre[ii*nxx + (jj-3)] + P_pre[ii*nxx + (jj+3)])
                  -  1008.0f*(P_pre[ii*nxx + (jj-2)] + P_pre[ii*nxx + (jj+2)])
                  +  8064.0f*(P_pre[ii*nxx + (jj+1)] + P_pre[ii*nxx + (jj-1)])
                  - 14350.0f*(P_pre[ii*nxx + jj]))/(5040.0f*powf(dx,2));

            /* Second derivative of the pressure with respect to z, 
            using eighth-order finite difference discretization */    
            d2_Pz2 = (- 9.0f*(P_pre[(ii-4)*nxx + jj] + P_pre[(ii+4)*nxx + jj])
                  +   128.0f*(P_pre[(ii-3)*nxx + jj] + P_pre[(ii+3)*nxx + jj])
                  -  1008.0f*(P_pre[(ii-2)*nxx + jj] + P_pre[(ii+2)*nxx + jj])
                  +  8064.0f*(P_pre[(ii-1)*nxx + jj] + P_pre[(ii+1)*nxx + jj])
                  - 14350.0f*(P_pre[ii*nxx + jj]))/(5040.0f*powf(dz,2));

            /* Calculating the future wave field through the previous fields */
            P_fut[ii*nxx + jj] = (powf(dt,2)*powf(vp[ii*nxx + jj],2)
                            * (d2_Px2 + d2_Pz2)) + 2.0f*P_pre[ii*nxx + jj] 
                            -  P_pas[ii*nxx + jj];
        }
    }
}

//
void FDM_8E2T_t_elastic_Isotropic2D(float *Vx, float *Vz, float *Txx, float *Tzz, float *Txz, float *rho, float *M, float *L, int nxx, int nzz, float dt, float dx, float dz) 
{
    int index, ii, jj;
    float rho_int, L_int, M_int;
    float d_Txx_dx, d_Txz_dz, d_Txz_dx, d_Tzz_dz;
    float d_Vx_dx, d_Vz_dz, d_Vx_dz, d_Vz_dx;
    
    for(index = 0; index < nxx*nzz; index++) 
    {                  
        ii = (int) index / nxx;      // indicador de linhas  (direção z)
        jj = (int) index % nxx;      // indicador de colunas (direção x)

        if((ii >= 3) && (ii < nzz-3) && (jj >= 3) && (jj < nxx-4)) 
        {
            d_Vx_dx = (75.0f*(Vx[(jj-3) + ii*nxx] - Vx[(jj+4) + ii*nxx]) + 
                     1029.0f*(Vx[(jj+3) + ii*nxx] - Vx[(jj-2) + ii*nxx]) +
                     8575.0f*(Vx[(jj-1) + ii*nxx] - Vx[(jj+2) + ii*nxx]) + 
                   128625.0f*(Vx[(jj+1) + ii*nxx] - Vx[jj + ii*nxx]))/(dx*107520.0f);

            d_Vz_dz = (75.0f*(Vz[jj + (ii-3)*nxx] - Vz[jj + (ii+4)*nxx]) +   
                     1029.0f*(Vz[jj + (ii+3)*nxx] - Vz[jj + (ii-2)*nxx]) +
                     8575.0f*(Vz[jj + (ii-1)*nxx] - Vz[jj + (ii+2)*nxx]) +
                   128625.0f*(Vz[jj + (ii+1)*nxx] - Vz[jj + ii*nxx]))/(dz*107520.0f);     

            Txx[index] += dt*((L[index] + 2*M[index])*d_Vx_dx + L[index]*d_Vz_dz);   
        
            Tzz[index] += dt*((L[index] + 2*M[index])*d_Vz_dz + L[index]*d_Vx_dx);
        }
    }

    for(index = 0; index < nxx*nzz; index++) 
    {                  
        ii = (int) index / nxx;      // indicador de linhas  (direção z)
        jj = (int) index % nxx;      // indicador de colunas (direção x)

        if((ii > 3) && (ii < nzz-3) && (jj > 3) && (jj < nxx-3)) 
        {
            d_Vz_dx = (75.0f*(Vz[(jj-4) + ii*nxx] - Vz[(jj+3) + ii*nxx]) +
                     1029.0f*(Vz[(jj+2) + ii*nxx] - Vz[(jj-3) + ii*nxx]) +
                     8575.0f*(Vz[(jj-2) + ii*nxx] - Vz[(jj+1) + ii*nxx]) +
                   128625.0f*(Vz[jj + ii*nxx]     - Vz[(jj-1) + ii*nxx]))/(dx*107520.0f);

            d_Vx_dz = (75.0f*(Vx[jj + (ii-4)*nxx] - Vx[jj + (ii+3)*nxx]) +
                     1029.0f*(Vx[jj + (ii+2)*nxx] - Vx[jj + (ii-3)*nxx]) +
                     8575.0f*(Vx[jj + (ii-2)*nxx] - Vx[jj + (ii+1)*nxx]) +
                   128625.0f*(Vx[jj + ii*nxx]     - Vx[jj + (ii-1)*nxx]))/(dz*107520.0f);

            M_int = powf(0.25*(1/M[jj + (ii+1)*nxx] + 1/M[(jj+1) + ii*nxx] + 1/M[(jj+1) + (ii+1)*nxx] + 1/M[jj + ii*nxx]),-1); 

            Txz[index] += dt * M_int * (d_Vx_dz + d_Vz_dx);            
        }      
    }

    for(index = 0; index < nxx*nzz; index++) 
    {              
        ii = (int) index / nxx;      // indicador de linhas  (direção z)
        jj = (int) index % nxx;      // indicador de colunas (direção x)

        if((ii >= 3) && (ii < nzz-4) && (jj > 3) && (jj < nxx-3)) 
        {
            d_Txx_dx = (75.0f*(Txx[(jj-4) + ii*nxx] - Txx[(jj+3) + ii*nxx]) +
                      1029.0f*(Txx[(jj+2) + ii*nxx] - Txx[(jj-3) + ii*nxx]) +
                      8575.0f*(Txx[(jj-2) + ii*nxx] - Txx[(jj+1) + ii*nxx]) +
                    128625.0f*(Txx[jj + ii*nxx]     - Txx[(jj-1) + ii*nxx]))/(dx*107520.0f);

            d_Txz_dz = (75.0f*(Txz[jj + (ii-3)*nxx] - Txz[jj + (ii+4)*nxx]) +
                      1029.0f*(Txz[jj + (ii+3)*nxx] - Txz[jj + (ii-2)*nxx]) + 
                      8575.0f*(Txz[jj + (ii-1)*nxx] - Txz[jj + (ii+2)*nxx]) +
                    128625.0f*(Txz[jj + (ii+1)*nxx] - Txz[jj + ii*nxx]))/(dz*107520.0f);

            rho_int = 0.5f * (rho[(jj+1) + ii*nxx] + rho[jj + ii*nxx]);

            Vx[index] += dt/rho_int*(d_Txx_dx + d_Txz_dz);  
        }
    }

    for(index = 0; index < nxx*nzz; index++) 
    {              
        ii = (int) index / nxx;      // indicador de linhas  (direção z)
        jj = (int) index % nxx;      // indicador de colunas (direção x)
      
        if((ii > 3) && (ii < nzz-3) && (jj >= 3) && (jj < nxx-4)) 
        {
            d_Txz_dx = (75.0f*(Txz[(jj-3) + ii*nxx] - Txz[(jj+4) + ii*nxx]) +
                      1029.0f*(Txz[(jj+3) + ii*nxx] - Txz[(jj-2) + ii*nxx]) +
                      8575.0f*(Txz[(jj-1) + ii*nxx] - Txz[(jj+2) + ii*nxx]) +
                    128625.0f*(Txz[(jj+1) + ii*nxx] - Txz[jj + ii*nxx]))/(dx*107520.0f);

            d_Tzz_dz = (75.0f*(Tzz[jj + (ii-4)*nxx] - Tzz[jj + (ii+3)*nxx]) + 
                      1029.0f*(Tzz[jj + (ii+2)*nxx] - Tzz[jj + (ii-3)*nxx]) +
                      8575.0f*(Tzz[jj + (ii-2)*nxx] - Tzz[jj + (ii+1)*nxx]) +
                    128625.0f*(Tzz[jj + ii*nxx]     - Tzz[jj + (ii-1)*nxx]))/(dz*107520.0f);

            rho_int = 0.5f * (rho[jj + (ii+1)*nxx] + rho[jj + ii*nxx]);

            Vz[index] += dt/rho_int*(d_Txz_dx + d_Tzz_dz); 
        }
    }
}

//
void FDM_8E2T_v_elastic_Isotropic2D(float *Vx, float *Vz, float *Txx, float *Tzz, float *Txz, float *rho, float *M, float *L, int nxx, int nzz, float dt, float dx, float dz) 
{
    int index, ii, jj;
    float rho_int, L_int, M_int;
    float d_Txx_dx, d_Txz_dz, d_Txz_dx, d_Tzz_dz;
    float d_Vx_dx, d_Vz_dz, d_Vx_dz, d_Vz_dx;

    for(index = 0; index < nxx*nzz; index++) 
    {              
        ii = (int) index / nxx;      // indicador de linhas  (direção z)
        jj = (int) index % nxx;      // indicador de colunas (direção x)

        if((ii >= 3) && (ii < nzz-4) && (jj >= 3) && (jj < nxx-4)) 
        {
            d_Txx_dx = (75.0f*(Txx[(jj-3) + ii*nxx] - Txx[(jj+4) + ii*nxx]) +
                      1029.0f*(Txx[(jj+3) + ii*nxx] - Txx[(jj-2) + ii*nxx]) +
                      8575.0f*(Txx[(jj-1) + ii*nxx] - Txx[(jj+2) + ii*nxx]) +
                    128625.0f*(Txx[(jj+1) + ii*nxx] - Txx[jj + ii*nxx]))/(dx*107520.0f);

            d_Txz_dz = (75.0f*(Txz[jj + (ii-3)*nxx] - Txz[jj + (ii+4)*nxx]) +
                      1029.0f*(Txz[jj + (ii+3)*nxx] - Txz[jj + (ii-2)*nxx]) + 
                      8575.0f*(Txz[jj + (ii-1)*nxx] - Txz[jj + (ii+2)*nxx]) +
                    128625.0f*(Txz[jj + (ii+1)*nxx] - Txz[jj + ii*nxx]))/(dz*107520.0f);

            Vx[index] += dt/rho[index]*(d_Txx_dx + d_Txz_dz);  
        }
    }

    for(index = 0; index < nxx*nzz; index++) 
    {              
        ii = (int) index / nxx;      // indicador de linhas  (direção z)
        jj = (int) index % nxx;      // indicador de colunas (direção x)
      
        if((ii > 3) && (ii < nzz-3) && (jj > 3) && (jj < nxx-3)) 
        {
            d_Txz_dx = (75.0f*(Txz[(jj-4) + ii*nxx] - Txz[(jj+3) + ii*nxx]) +
                      1029.0f*(Txz[(jj+2) + ii*nxx] - Txz[(jj-3) + ii*nxx]) +
                      8575.0f*(Txz[(jj-2) + ii*nxx] - Txz[(jj+1) + ii*nxx]) +
                    128625.0f*(Txz[jj + ii*nxx]     - Txz[(jj-1) + ii*nxx]))/(dx*107520.0f);

            d_Tzz_dz = (75.0f*(Tzz[jj + (ii-4)*nxx] - Tzz[jj + (ii+3)*nxx]) + 
                      1029.0f*(Tzz[jj + (ii+2)*nxx] - Tzz[jj + (ii-3)*nxx]) +
                      8575.0f*(Tzz[jj + (ii-2)*nxx] - Tzz[jj + (ii+1)*nxx]) +
                    128625.0f*(Tzz[jj + ii*nxx]     - Tzz[jj + (ii-1)*nxx]))/(dz*107520.0f);

            rho_int = powf(0.25*(1/rho[jj + ii*nxx] + 1/rho[jj + (ii+1)*nxx] + 1/rho[(jj+1) + (ii+1)*nxx] + 1/rho[(jj+1) + ii*nxx]),-1);

            Vz[index] += dt/rho_int*(d_Txz_dx + d_Tzz_dz); 
        }
    }
    
    for(index = 0; index < nxx*nzz; index++) 
    {                  
        ii = (int) index / nxx;      // indicador de linhas  (direção z)
        jj = (int) index % nxx;      // indicador de colunas (direção x)

        if((ii >= 3) && (ii < nzz-4) && (jj > 3) && (jj < nxx-3)) 
        {
            d_Vx_dx = (75.0f*(Vx[(jj-4) + ii*nxx] - Vx[(jj+3) + ii*nxx]) + 
                     1029.0f*(Vx[(jj+2) + ii*nxx] - Vx[(jj-3) + ii*nxx]) +
                     8575.0f*(Vx[(jj-2) + ii*nxx] - Vx[(jj+1) + ii*nxx]) + 
                   128625.0f*(Vx[jj + ii*nxx]     - Vx[(jj-1) + ii*nxx]))/(dx*107520.0f);

            d_Vz_dz = (75.0f*(Vz[jj + (ii-3)*nxx] - Vz[jj + (ii+4)*nxx]) +   
                     1029.0f*(Vz[jj + (ii+3)*nxx] - Vz[jj + (ii-2)*nxx]) +
                     8575.0f*(Vz[jj + (ii-1)*nxx] - Vz[jj + (ii+2)*nxx]) +
                   128625.0f*(Vz[jj + (ii+1)*nxx] - Vz[jj + ii*nxx]))/(dz*107520.0f);     

            L_int = 0.5*(L[(jj+1) + ii*nxx] + L[jj + ii*nxx]);
            M_int = 0.5*(M[(jj+1) + ii*nxx] + M[jj + ii*nxx]);

            Txx[index] += dt*((L_int + 2*M_int)*d_Vx_dx + L_int*d_Vz_dz);   
        
            Tzz[jj+ ii*nxx] += dt*((L_int + 2*M_int)*d_Vz_dz + L_int*d_Vx_dx);
        }
    }

    for(index = 0; index < nxx*nzz; index++) 
    {                  
        ii = (int) index / nxx;      // indicador de linhas  (direção z)
        jj = (int) index % nxx;      // indicador de colunas (direção x)

        if((ii >= 3) && (ii < nzz-4) && (jj >= 3) && (jj < nxx-4)) 
        {
            d_Vz_dx = (75.0f*(Vz[(jj-3) + ii*nxx] - Vz[(jj+4) + ii*nxx]) +
                     1029.0f*(Vz[(jj+3) + ii*nxx] - Vz[(jj-2) + ii*nxx]) +
                     8575.0f*(Vz[(jj-1) + ii*nxx] - Vz[(jj+2) + ii*nxx]) +
                   128625.0f*(Vz[(jj+1) + ii*nxx] - Vz[jj + ii*nxx]))/dx;

            d_Vx_dz = (75.0f*(Vx[jj + (ii-4)*nxx] - Vx[jj + (ii+3)*nxx]) +
                     1029.0f*(Vx[jj + (ii+2)*nxx] - Vx[jj + (ii-3)*nxx]) +
                     8575.0f*(Vx[jj + (ii-2)*nxx] - Vx[jj + (ii+1)*nxx]) +
                   128625.0f*(Vx[jj + ii*nxx]     - Vx[jj + (ii-1)*nxx]))/dz;

            M_int = 0.5*(M[jj + (ii+1)*nxx] + M[jj + ii*nxx]); 

            Txz[index] += M_int*dt/107520.0f*(d_Vx_dz + d_Vz_dx);            
        }      
    }
}

//
void cerjanElasticAbsorbingCondition2D(float *Vx,float *Vz,float *Txx,float *Tzz,float *Txz, float *factor, int nxx, int nzz, int abc_layer, float * corner)
{
    int ii,jj,index;
    #pragma acc parallel loop present(P_pre[0:nxx*nzz],P_fut[0:nxx*nzz],factor[0:abc_layer],corner[0:abc_layer*abc_layer]) 
    for(index = 0; index < nxx*nzz; ++index) 
    {     
        ii = floor(index / nxx);  /* Line indicator */
        jj = index % nxx;         /* Column indicator */

        /* Left and right absorbing layers */ 
        if((ii >= abc_layer) && (ii < nzz-abc_layer) && (jj >= 0) && (jj < abc_layer)) 
        {
            Vx[ii*nxx+jj] *= factor[jj];         /* Left  */
            Vz[ii*nxx+jj] *= factor[jj];         /* Left  */
            Txx[ii*nxx+jj] *= factor[jj];        /* Left  */
            Tzz[ii*nxx+jj] *= factor[jj];        /* Left  */
            Txz[ii*nxx+jj] *= factor[jj];        /* Left  */
                    
            Vx[ii*nxx+(nxx-jj)] *= factor[jj];   /* Right */
            Vz[ii*nxx+(nxx-jj)] *= factor[jj];   /* Right */
            Txx[ii*nxx+(nxx-jj)] *= factor[jj];  /* Right */
            Tzz[ii*nxx+(nxx-jj)] *= factor[jj];  /* Right */
            Txz[ii*nxx+(nxx-jj)] *= factor[jj];  /* Right */
        }

        /* Up and bottom absorbing layers */
        if((jj >= abc_layer) && (jj < nxx - abc_layer) && (ii >= 0) && (ii < abc_layer)) 
        {
            Vx[ii*nxx+jj] *= factor[ii];         /*   Up   */
            Vz[ii*nxx+jj] *= factor[ii];         /*   Up   */
            Txx[ii*nxx+jj] *= factor[ii];        /*   Up   */
            Tzz[ii*nxx+jj] *= factor[ii];        /*   Up   */
            Txz[ii*nxx+jj] *= factor[ii];        /*   Up   */

            Vx[(nzz-ii)*nxx-jj] *= factor[ii];   /* Bottom */ 
            Vz[(nzz-ii)*nxx-jj] *= factor[ii];   /* Bottom */ 
            Txx[(nzz-ii)*nxx-jj] *= factor[ii];  /* Bottom */ 
            Tzz[(nzz-ii)*nxx-jj] *= factor[ii];  /* Bottom */ 
            Txz[(nzz-ii)*nxx-jj] *= factor[ii];  /* Bottom */ 
        }
        
        /* Upper left absorbing corner */ 
        if((ii >= 0) && (ii < abc_layer) && (jj >= 0) && (jj < abc_layer)) 
        {
            Vx[ii*nxx + jj] *= corner[ii*abc_layer + jj];
            Vz[ii*nxx + jj] *= corner[ii*abc_layer + jj];
            Txx[ii*nxx + jj] *= corner[ii*abc_layer + jj];
            Tzz[ii*nxx + jj] *= corner[ii*abc_layer + jj];
            Txz[ii*nxx + jj] *= corner[ii*abc_layer + jj];
        }

        /* Upper right absorbing corner */
        if((ii >= 0) && (ii < abc_layer) && (jj >= nxx-abc_layer) && (jj < nxx)) 
        {
            Vx[ii*nxx + jj] *= corner[ii*abc_layer + (nxx-jj-1)];
            Vz[ii*nxx + jj] *= corner[ii*abc_layer + (nxx-jj-1)];
            Txx[ii*nxx + jj] *= corner[ii*abc_layer + (nxx-jj-1)];
            Tzz[ii*nxx + jj] *= corner[ii*abc_layer + (nxx-jj-1)];
            Txz[ii*nxx + jj] *= corner[ii*abc_layer + (nxx-jj-1)];
        }

        /* Bottom left absorbing corner */
        if((ii >= nzz-abc_layer) && (ii < nzz) && (jj >= 0) && (jj < abc_layer)) 
        {
            Vx[ii*nxx + jj] *= corner[(nzz-ii-1)*abc_layer + jj];
            Vz[ii*nxx + jj] *= corner[(nzz-ii-1)*abc_layer + jj];
            Txx[ii*nxx + jj] *= corner[(nzz-ii-1)*abc_layer + jj];
            Tzz[ii*nxx + jj] *= corner[(nzz-ii-1)*abc_layer + jj];
            Txz[ii*nxx + jj] *= corner[(nzz-ii-1)*abc_layer + jj];
        }

        /* Botton right absorbing corner */ 
        if((ii >= nzz-abc_layer) && (ii < nzz) && (jj >= nxx-abc_layer) && (jj < nxx)) 
        {
            Vx[ii*nxx + jj] *= corner[(nzz-ii-1)*abc_layer + (nxx-jj-1)];
            Vz[ii*nxx + jj] *= corner[(nzz-ii-1)*abc_layer + (nxx-jj-1)];
            Txx[ii*nxx + jj] *= corner[(nzz-ii-1)*abc_layer + (nxx-jj-1)];
            Tzz[ii*nxx + jj] *= corner[(nzz-ii-1)*abc_layer + (nxx-jj-1)];
            Txz[ii*nxx + jj] *= corner[(nzz-ii-1)*abc_layer + (nxx-jj-1)];
        }
    }
}

//
void FDM_8E2T_acoustic3D(float *vp,float *P_pas,float *P_pre,float *P_fut,int nxx,int nyy,int nzz,float dx,float dy,float dz,float dt)
{
    int ii,jj,kk,index;
    float d2_Px2,d2_Py2,d2_Pz2;

    for(index = 0; index < nxx*nyy*nzz; ++index)   /* Spatial loop */ 
    {
        kk = floor(index/(nxx*nzz));               // indicador de matrizes (direção y)
        jj = index % nxx;                          // indicador de colunas  (direção x)
        ii = floor((index % (nxx*nzz)) / nxx);     // indicador de linhas   (direção z)  

        if((kk > 3) && (kk < nyy-4) && (ii > 3) && (ii < nzz-4) && (jj > 3) && (jj < nxx-4)) 
        {
            /* Second derivative of the pressure with respect to x, 
            using eighth-order finite difference discretization */
            d2_Px2 = (- 9.0f*(P_pre[kk*nxx*nzz + ii*nxx + (jj-4)] + P_pre[kk*nxx*nzz + ii*nxx + (jj+4)])
                  +   128.0f*(P_pre[kk*nxx*nzz + ii*nxx + (jj-3)] + P_pre[kk*nxx*nzz + ii*nxx + (jj+3)])
                  -  1008.0f*(P_pre[kk*nxx*nzz + ii*nxx + (jj-2)] + P_pre[kk*nxx*nzz + ii*nxx + (jj+2)])
                  +  8064.0f*(P_pre[kk*nxx*nzz + ii*nxx + (jj-1)] + P_pre[kk*nxx*nzz + ii*nxx + (jj+1)])
                  - 14350.0f*(P_pre[kk*nxx*nzz + ii*nxx + jj]))/(5040.0f*powf(dx,2));

            /* Second derivative of the pressure with respect to y, 
            using eighth-order finite difference discretization */
            d2_Py2 = (- 9.0f*(P_pre[(kk-4)*nxx*nzz + ii*nxx + jj] + P_pre[(kk+4)*nxx*nzz + ii*nxx + jj])
                  +   128.0f*(P_pre[(kk-3)*nxx*nzz + ii*nxx + jj] + P_pre[(kk+3)*nxx*nzz + ii*nxx + jj])
                  -  1008.0f*(P_pre[(kk-2)*nxx*nzz + ii*nxx + jj] + P_pre[(kk+2)*nxx*nzz + ii*nxx + jj])
                  +  8064.0f*(P_pre[(kk-1)*nxx*nzz + ii*nxx + jj] + P_pre[(kk+1)*nxx*nzz + ii*nxx + jj])
                  - 14350.0f*(P_pre[kk*nxx*nzz + ii*nxx + jj]))/(5040.0f*powf(dx,2));

            /* Second derivative of the pressure with respect to z, 
            using eighth-order finite difference discretization */    
            d2_Pz2 = (- 9.0f*(P_pre[kk*nxx*nzz + (ii-4)*nxx + jj] + P_pre[kk*nxx*nzz + (ii+4)*nxx + jj])
                  +   128.0f*(P_pre[kk*nxx*nzz + (ii-3)*nxx + jj] + P_pre[kk*nxx*nzz + (ii+3)*nxx + jj])
                  -  1008.0f*(P_pre[kk*nxx*nzz + (ii-2)*nxx + jj] + P_pre[kk*nxx*nzz + (ii+2)*nxx + jj])
                  +  8064.0f*(P_pre[kk*nxx*nzz + (ii-1)*nxx + jj] + P_pre[kk*nxx*nzz + (ii+1)*nxx + jj])
                  - 14350.0f*(P_pre[kk*nxx*nzz + ii*nxx + jj]))/(5040.0f*powf(dz,2));

            /* Calculating the future wave field through the previous fields */
            P_fut[kk*nxx*nzz + ii*nxx + jj] = (powf(dt,2)*powf(vp[kk*nxx*nzz + ii*nxx + jj],2)
                            * (d2_Px2 + d2_Py2 + d2_Pz2)) + 2.0f*P_pre[kk*nxx*nzz + ii*nxx + jj] 
                            -  P_pas[kk*nxx*nzz + ii*nxx + jj];
        }
    }
}

//
void FDM_8E2T_pvAcoustic3D(float *P, float *vx, float *vy, float *vz,float *K, float *rho, float dx, float dy, float dz, float dt, int nx, int ny, int nz)
{
    int index, ii, jj, kk;
    float dvx_dx, dvy_dy, dvz_dz;
    float dP_dx, dP_dy, dP_dz;
    float rhox, rhoy, rhoz;

    for (index = 0; index < nx*ny*nz; index++)
    {
        kk = floor(index/(nx*nz));          // indicador de matrizes (direção y)
        jj = index % nx;                    // indicador de colunas  (direção x)
        ii = floor((index % (nx*nz)) / nx); // indicador de linhas   (direção z)  

        if((ii > 3) && (ii < nz-3) && (kk > 3) && (kk < ny-4) && (jj > 3) && (jj < nx-3)) 
        { 
            dvx_dx = (75.0f*(vx[kk*nz*nx + ii*nx + (jj-4)] - vx[kk*nz*nx + ii*nx + (jj+3)]) +
                    1029.0f*(vx[kk*nz*nx + ii*nx + (jj+2)] - vx[kk*nz*nx + ii*nx + (jj-3)]) +
                    8575.0f*(vx[kk*nz*nx + ii*nx + (jj-2)] - vx[kk*nz*nx + ii*nx + (jj+1)]) +
                  128625.0f*(vx[kk*nz*nx + ii*nx + jj]     - vx[kk*nz*nx + ii*nx + (jj-1)]))/(107520.0f*dx);
            
            dvy_dy = (75.0f*(vy[(kk-4)*nz*nx + ii*nx + jj] - vy[(kk+3)*nz*nx + ii*nx + jj]) +
                    1029.0f*(vy[(kk+2)*nz*nx + ii*nx + jj] - vy[(kk-3)*nz*nx + ii*nx + jj]) +
                    8575.0f*(vy[(kk-2)*nz*nx + ii*nx + jj] - vy[(kk+1)*nz*nx + ii*nx + jj]) +
                  128625.0f*(vy[kk*nz*nx + ii*nx + jj]     - vy[(kk-1)*nz*nx + ii*nx + jj]))/(107520.0f*dy);
            
            dvz_dz = (75.0f*(vz[kk*nz*nx + (ii-4)*nx + jj] - vz[kk*nz*nx + (ii+3)*nx + jj]) +
                    1029.0f*(vz[kk*nz*nx + (ii+2)*nx + jj] - vz[kk*nz*nx + (ii-3)*nx + jj]) +
                    8575.0f*(vz[kk*nz*nx + (ii-2)*nx + jj] - vz[kk*nz*nx + (ii+1)*nx + jj]) +
                  128625.0f*(vz[kk*nz*nx + ii*nx + jj]     - vz[kk*nz*nx + (ii-1)*nx + jj]))/(107520.0f*dz);

            P[index] += dt*(- K[index]*(dvx_dx + dvy_dy + dvz_dz));
        }
    }

    for (index = 0; index < nx*ny*nz; index++)
    {
        kk = floor(index/(nx*nz));          // indicador de matrizes (direção y)
        jj = index % nx;                    // indicador de colunas  (direção x)
        ii = floor((index % (nx*nz)) / nx); // indicador de linhas   (direção z)  

        if ((ii >= 0) && (ii < nz) && (kk >= 0) && (kk < ny) && (jj >= 3) && (jj < nx-4))
        {
            dP_dx = (75.0f*(P[kk*nz*nx + ii*nx + (jj-3)] - P[kk*nz*nx + ii*nx + (jj+4)]) +
                   1029.0f*(P[kk*nz*nx + ii*nx + (jj+3)] - P[kk*nz*nx + ii*nx + (jj-2)]) +
                   8575.0f*(P[kk*nz*nx + ii*nx + (jj-1)] - P[kk*nz*nx + ii*nx + (jj+2)]) +
                 128625.0f*(P[kk*nz*nx + ii*nx + (jj+1)] - P[kk*nz*nx + ii*nx + jj]))/(107520.0f*dx); 

            rhox = (rho[kk*nz*nx + ii*nx + (jj+1)] + rho[kk*nz*nx + ii*nx + jj]) / 2.0f;

            vx[index] += dt*(- dP_dx / rhox);        
        }

        if ((ii >= 0) && (ii < nz) && (kk >= 3) && (kk < ny-4) && (jj >= 0) && (jj < nx))
        {
            dP_dy = (75.0f*(P[(kk-3)*nz*nx + ii*nx + jj] - P[(kk+4)*nz*nx + ii*nx + jj]) +
                   1029.0f*(P[(kk+3)*nz*nx + ii*nx + jj] - P[(kk-2)*nz*nx + ii*nx + jj]) +
                   8575.0f*(P[(kk-1)*nz*nx + ii*nx + jj] - P[(kk+2)*nz*nx + ii*nx + jj]) +
                 128625.0f*(P[(kk+1)*nz*nx + ii*nx + jj] - P[kk*nz*nx + ii*nx + jj]))/(107520.0f*dy);

            rhoy = (rho[(kk+1)*nz*nx + ii*nx + jj] + rho[kk*nz*nx + ii*nx + jj]) / 2.0f;

            vy[index] += dt*(- dP_dy / rhoy);        
        }

        if ((ii >= 3) && (ii < nz-3) && (kk >= 0) && (kk < ny) && (jj >= 0) && (jj < nx))
        {
            dP_dz = (75.0f*(P[kk*nz*nx + (ii-3)*nx + jj] - P[kk*nz*nx + (ii+4)*nx + jj]) +
                   1029.0f*(P[kk*nz*nx + (ii+3)*nx + jj] - P[kk*nz*nx + (ii-2)*nx + jj]) +
                   8575.0f*(P[kk*nz*nx + (ii-1)*nx + jj] - P[kk*nz*nx + (ii+2)*nx + jj]) +
                 128625.0f*(P[kk*nz*nx + (ii+1)*nx + jj] - P[kk*nz*nx + ii*nx + jj]))/(107520.0f*dz);

            rhoz = (rho[kk*nz*nx + (ii+1)*nx + jj] + rho[kk*nz*nx + ii*nx + jj]) / 2.0f;

            vz[index] += dt*(- dP_dz / rhoz);        
        }
    }
}

// Elástico Isotrópico 3D com matrizes lineares (Graves,1996) - Paulo Bastos, GISIS 22/07/2019 (GRAVES,1996)
void FDM_8E2T_elastic_Isotropic3D(float *Vx,float *Vy,float *Vz,float *Txx,float *Tyy,float *Tzz,float *Txy,float *Txz,float *Tyz, float *rho, float *M, float *L, int nx,int ny,int nz,float dt,float dx,float dy,float dz) 
{
    float d_Txx_dx, d_Txy_dy, d_Txz_dz, d_Vx_dz, d_Vx_dy, d_Vx_dx;
    float d_Txy_dx, d_Tyy_dy, d_Tyz_dz, d_Vy_dx, d_Vy_dz, d_Vy_dy;
    float d_Txz_dx, d_Tyz_dy, d_Tzz_dz, d_Vz_dx, d_Vz_dy, d_Vz_dz;

    int index, kk, jj, ii;

    float rhox, rhoy, rhoz; // Interpolated density 
    float M_xy, M_xz, M_yz; // Interpolated rigity

    for(index = 0; index < nx*ny*nz; index++) 
    {    
        kk = floor(index/(nx*nz));          // indicador de matrizes (direção y)
        jj = index % nx;                    // indicador de colunas  (direção x)
        ii = floor((index % (nx*nz)) / nx); // indicador de linhas   (direção z)  

        if((ii >= 3) && (ii < nz-4) && (kk >= 3) && (kk < ny-4) && (jj >= 3) && (jj < nx-4)) 
        {    
            d_Vx_dx = (75.0f*(Vx[kk*nz*nx + ii*nx + (jj-3)] - Vx[kk*nz*nx + ii*nx + (jj+4)]) +
                     1029.0f*(Vx[kk*nz*nx + ii*nx + (jj+3)] - Vx[kk*nz*nx + ii*nx + (jj-2)]) +
                     8575.0f*(Vx[kk*nz*nx + ii*nx + (jj-1)] - Vx[kk*nz*nx + ii*nx + (jj+2)]) +
                   128625.0f*(Vx[kk*nz*nx + ii*nx + (jj+1)] - Vx[kk*nz*nx + ii*nx + jj]))/(107520.0f*dx);

            d_Vy_dy = (75.0f*(Vy[(kk-3)*nz*nx + ii*nx + jj] - Vy[(kk+4)*nz*nx + ii*nx + jj]) +
                     1029.0f*(Vy[(kk+3)*nz*nx + ii*nx + jj] - Vy[(kk-2)*nz*nx + ii*nx + jj]) +
                     8575.0f*(Vy[(kk-1)*nz*nx + ii*nx + jj] - Vy[(kk+2)*nz*nx + ii*nx + jj]) +
                   128625.0f*(Vy[(kk+1)*nz*nx + ii*nx + jj] - Vy[kk*nz*nx + ii*nx + jj]))/(107520.0f*dy);

            d_Vz_dz = (75.0f*(Vz[kk*nz*nx + (ii-3)*nx + jj] - Vz[kk*nz*nx + (ii+4)*nx + jj]) +
                     1029.0f*(Vz[kk*nz*nx + (ii+3)*nx + jj] - Vz[kk*nz*nx + (ii-2)*nx + jj]) +
                     8575.0f*(Vz[kk*nz*nx + (ii-1)*nx + jj] - Vz[kk*nz*nx + (ii+2)*nx + jj]) +
                   128625.0f*(Vz[kk*nz*nx + (ii+1)*nx + jj] - Vz[kk*nz*nx + ii*nx + jj]))/(107520.0f*dz);

            Txx[index] += dt *(L[index] + 2*M[index])*d_Vx_dx +
                          dt * L[index]*(d_Vy_dy + d_Vz_dz);

            Tyy[index] += dt *(L[index] + 2*M[index])*d_Vy_dy +
                          dt * L[index]*(d_Vx_dx + d_Vz_dz);

            Tzz[index] += dt *(L[index] + 2*M[index])*d_Vz_dz +
                          dt * L[index]*(d_Vx_dx + d_Vy_dy);                    
        }

        if((ii >= 0) && (ii < nz) && (kk > 3) && (kk < ny-3) && (jj > 3) && (jj < nx-3)) 
        {
            d_Vx_dy = (75.0f*(Vx[(kk-4)*nz*nx + ii*nx + jj] - Vx[(kk+3)*nz*nx + ii*nx + jj]) +
                     1029.0f*(Vx[(kk+2)*nz*nx + ii*nx + jj] - Vx[(kk-3)*nz*nx + ii*nx + jj]) +
                     8575.0f*(Vx[(kk-2)*nz*nx + ii*nx + jj] - Vx[(kk+1)*nz*nx + ii*nx + jj]) +
                   128625.0f*(Vx[kk*nz*nx + ii*nx + jj]     - Vx[(kk-1)*nz*nx + ii*nx + jj]))/(107520.0f*dy);

            d_Vy_dx = (75.0f*(Vy[kk*nz*nx + ii*nx + (jj-4)] - Vy[kk*nz*nx + ii*nx + (jj+3)]) +
                     1029.0f*(Vy[kk*nz*nx + ii*nx + (jj+2)] - Vy[kk*nz*nx + ii*nx + (jj-3)]) +
                     8575.0f*(Vy[kk*nz*nx + ii*nx + (jj-2)] - Vy[kk*nz*nx + ii*nx + (jj+1)]) +
                   128625.0f*(Vy[kk*nz*nx + ii*nx + jj]     - Vy[kk*nz*nx + ii*nx + (jj-1)]))/(107520.0f*dx);

            M_xy = powf(0.25*(1/M[(kk+1)*nz*nx + ii*nx + (jj+1)] +
                              1/M[kk*nz*nx + ii*nx + (jj+1)] + 
                              1/M[(kk+1)*nz*nx + ii*nx + jj] +
                              1/M[kk*nz*nx + ii*nx + jj]),-1.0f);

            Txy[index] += dt*M_xy*(d_Vx_dy + d_Vy_dx);
        }

        if((ii > 3) && (ii < nz-3) && (kk >= 0) && (kk < ny) && (jj > 3) && (jj < nx-3)) 
        {
            d_Vx_dz = (75.0f*(Vx[kk*nz*nx + (ii-4)*nx + jj] - Vx[kk*nz*nx + (ii+3)*nx + jj]) +
                     1029.0f*(Vx[kk*nz*nx + (ii+2)*nx + jj] - Vx[kk*nz*nx + (ii-3)*nx + jj]) +
                     8575.0f*(Vx[kk*nz*nx + (ii-2)*nx + jj] - Vx[kk*nz*nx + (ii+1)*nx + jj]) +
                   128625.0f*(Vx[kk*nz*nx + ii*nx + jj]     - Vx[kk*nz*nx + (ii-1)*nx + jj]))/(107520.0f*dz);

            d_Vz_dx = (75.0f*(Vz[kk*nz*nx + ii*nx + (jj-4)] - Vz[kk*nz*nx + ii*nx + (jj+3)]) +
                     1029.0f*(Vz[kk*nz*nx + ii*nx + (jj+2)] - Vz[kk*nz*nx + ii*nx + (jj-3)]) +
                     8575.0f*(Vz[kk*nz*nx + ii*nx + (jj-2)] - Vz[kk*nz*nx + ii*nx + (jj+1)]) +
                   128625.0f*(Vz[kk*nz*nx + ii*nx + jj]     - Vz[kk*nz*nx + ii*nx + (jj-1)]))/(107520.0f*dx);

            M_xz = powf(0.25*(1/M[kk*nz*nx + (ii+1)*nx + (jj+1)] +
                              1/M[kk*nz*nx + ii*nx + (jj+1)] + 
                              1/M[kk*nz*nx + (ii+1)*nx + jj] +
                              1/M[kk*nz*nx + ii*nx + jj]),-1.0f);

            Txz[index] += dt*M_xz*(d_Vx_dz + d_Vz_dx);
        }
    
        if((ii > 3) && (ii < nz-3) && (kk > 3) && (kk < ny-3) && (jj >= 0) && (jj < nx)) 
        {
            d_Vy_dz = (75.0f*(Vy[kk*nz*nx + (ii-4)*nx + jj] - Vy[kk*nz*nx + (ii+3)*nx + jj]) +
                     1029.0f*(Vy[kk*nz*nx + (ii+2)*nx + jj] - Vy[kk*nz*nx + (ii-3)*nx + jj]) +
                     8575.0f*(Vy[kk*nz*nx + (ii-2)*nx + jj] - Vy[kk*nz*nx + (ii+1)*nx + jj]) +
                   128625.0f*(Vy[kk*nz*nx + ii*nx + jj]     - Vy[kk*nz*nx + (ii-1)*nx + jj]))/(107520.0f*dz);

            d_Vz_dy = (75.0f*(Vz[(kk-4)*nz*nx + ii*nx + jj] - Vz[(kk+3)*nz*nx + ii*nx + jj]) +
                     1029.0f*(Vz[(kk+2)*nz*nx + ii*nx + jj] - Vz[(kk-3)*nz*nx + ii*nx + jj]) +
                     8575.0f*(Vz[(kk-2)*nz*nx + ii*nx + jj] - Vz[(kk+1)*nz*nx + ii*nx + jj]) +
                   128625.0f*(Vz[kk*nz*nx + ii*nx + jj]     - Vz[(kk-1)*nz*nx + ii*nx + jj]))/(107520.0f*dy);

            M_yz = powf(0.25*(1/M[(kk+1)*nz*nx + (ii+1)*nx + jj] +
                              1/M[(kk+1)*nz*nx + ii*nx + jj] + 
                              1/M[kk*nz*nx + (ii+1)*nx + jj] +
                              1/M[kk*nz*nx + ii*nx + jj]),-1.0f);

            Tyz[index] += dt*M_yz*(d_Vy_dz + d_Vz_dy);
        }
    }

    for(index = 0; index < nx*ny*nz; index++) 
    {    
        kk = floor(index/(nx*nz));          // indicador de matrizes (direção y)
        jj = index % nx;                    // indicador de colunas  (direção x)
        ii = floor((index % (nx*nz)) / nx); // indicador de linhas   (direção z)  

        if((ii >= 3) && (ii < nz-4) && (kk >= 3) && (kk < ny-4) && (jj > 3) && (jj < nx-3)) 
        {
            d_Txx_dx = (75.0f*(Txx[kk*nz*nx + ii*nx + (jj-4)] - Txx[kk*nz*nx + ii*nx + (jj+3)]) +
                      1029.0f*(Txx[kk*nz*nx + ii*nx + (jj+2)] - Txx[kk*nz*nx + ii*nx + (jj-3)]) +
                      8575.0f*(Txx[kk*nz*nx + ii*nx + (jj-2)] - Txx[kk*nz*nx + ii*nx + (jj+1)]) +
                    128625.0f*(Txx[kk*nz*nx + ii*nx + jj]     - Txx[kk*nz*nx + ii*nx + (jj-1)]))/(107520.0f*dx);

            d_Txy_dy = (75.0f*(Txy[(kk-3)*nz*nx + ii*nx + jj] - Txy[(kk+4)*nz*nx + ii*nx + jj]) +
                      1029.0f*(Txy[(kk+3)*nz*nx + ii*nx + jj] - Txy[(kk-2)*nz*nx + ii*nx + jj]) +
                      8575.0f*(Txy[(kk-1)*nz*nx + ii*nx + jj] - Txy[(kk+2)*nz*nx + ii*nx + jj]) +
                    128625.0f*(Txy[(kk+1)*nz*nx + ii*nx + jj] - Txy[kk*nz*nx + ii*nx + jj]))/(107520.0f*dy);

            d_Txz_dz = (75.0f*(Txz[kk*nz*nx + (ii-3)*nx + jj] - Txz[kk*nz*nx + (ii+4)*nx + jj]) +
                      1029.0f*(Txz[kk*nz*nx + (ii+3)*nx + jj] - Txz[kk*nz*nx + (ii-2)*nx + jj]) +
                      8575.0f*(Txz[kk*nz*nx + (ii-1)*nx + jj] - Txz[kk*nz*nx + (ii+2)*nx + jj]) +
                    128625.0f*(Txz[kk*nz*nx + (ii+1)*nx + jj] - Txz[kk*nz*nx + ii*nx + jj]))/(107520.0f*dz);

            rhox = 0.5f * (rho[kk*nz*nx + ii*nx + (jj+1)] + rho[kk*nz*nx + ii*nx + jj]);

            Vx[index] += dt/rhox*(d_Txx_dx + d_Txy_dy + d_Txz_dz); 
        }
    
        if((ii >= 3) && (ii < nz-3) && (kk > 3) && (kk < ny-3) && (jj >= 3) && (jj < nx-4)) 
        {
            d_Txy_dx = (75.0f*(Txy[kk*nz*nx + ii*nx + (jj-3)] - Txy[kk*nz*nx + ii*nx + (jj+4)]) +
                      1029.0f*(Txy[kk*nz*nx + ii*nx + (jj+3)] - Txy[kk*nz*nx + ii*nx + (jj-2)]) +
                      8575.0f*(Txy[kk*nz*nx + ii*nx + (jj-1)] - Txy[kk*nz*nx + ii*nx + (jj+2)]) +
                    128625.0f*(Txy[kk*nz*nx + ii*nx + (jj+1)] - Txy[kk*nz*nx + ii*nx + jj]))/(107520.0f*dx);

            d_Tyy_dy = (75.0f*(Tyy[(kk-4)*nz*nx + ii*nx + jj] - Tyy[(kk+3)*nz*nx + ii*nx + jj]) +
                      1029.0f*(Tyy[(kk+2)*nz*nx + ii*nx + jj] - Tyy[(kk-3)*nz*nx + ii*nx + jj]) +
                      8575.0f*(Tyy[(kk-2)*nz*nx + ii*nx + jj] - Tyy[(kk+1)*nz*nx + ii*nx + jj]) +
                    128625.0f*(Tyy[kk*nz*nx + ii*nx + jj]     - Tyy[(kk-1)*nz*nx + ii*nx + jj]))/(107520.0f*dy);

            d_Tyz_dz = (75.0f*(Tyz[kk*nz*nx + (ii-3)*nx + jj] - Tyz[kk*nz*nx + (ii+4)*nx + jj]) +
                      1029.0f*(Tyz[kk*nz*nx + (ii+3)*nx + jj] - Tyz[kk*nz*nx + (ii-2)*nx + jj]) +
                      8575.0f*(Tyz[kk*nz*nx + (ii-1)*nx + jj] - Tyz[kk*nz*nx + (ii+2)*nx + jj]) +
                    128625.0f*(Tyz[kk*nz*nx + (ii+1)*nx + jj] - Tyz[kk*nz*nx + ii*nx + jj]))/(107520.0f*dz);

            rhoy = 0.5f * (rho[(kk+1)*nz*nx + ii*nx + jj] + rho[kk*nz*nx + ii*nx + jj]);

            Vy[index] += dt/rhoy*(d_Txy_dx + d_Tyy_dy + d_Tyz_dz); 
        }    

        if((ii > 3) && (ii < nz-3) && (kk >= 3) && (kk < ny-4) && (jj >= 3) && (jj < nx-4)) 
        {
            d_Txz_dx = (75.0f*(Txz[kk*nz*nx + ii*nx + (jj-3)] - Txz[kk*nz*nx + ii*nx + (jj+4)]) +
                      1029.0f*(Txz[kk*nz*nx + ii*nx + (jj+3)] - Txz[kk*nz*nx + ii*nx + (jj-2)]) +
                      8575.0f*(Txz[kk*nz*nx + ii*nx + (jj-1)] - Txz[kk*nz*nx + ii*nx + (jj+2)]) +
                    128625.0f*(Txz[kk*nz*nx + ii*nx + (jj+1)] - Txz[kk*nz*nx + ii*nx + jj]))/(107520.0f*dx);

            d_Tyz_dy = (75.0f*(Tyz[(kk-3)*nz*nx + ii*nx + jj] - Tyz[(kk+4)*nz*nx + ii*nx + jj]) +
                      1029.0f*(Tyz[(kk+3)*nz*nx + ii*nx + jj] - Tyz[(kk-2)*nz*nx + ii*nx + jj]) +
                      8575.0f*(Tyz[(kk-1)*nz*nx + ii*nx + jj] - Tyz[(kk+2)*nz*nx + ii*nx + jj]) +
                    128625.0f*(Tyz[(kk+1)*nz*nx + ii*nx + jj] - Tyz[kk*nz*nx + ii*nx + jj]))/(107520.0f*dy);

            d_Tzz_dz = (75.0f*(Tzz[kk*nz*nx + (ii-4)*nx + jj] - Tzz[kk*nz*nx + (ii+3)*nx + jj]) +
                      1029.0f*(Tzz[kk*nz*nx + (ii+2)*nx + jj] - Tzz[kk*nz*nx + (ii-3)*nx + jj]) +
                      8575.0f*(Tzz[kk*nz*nx + (ii-2)*nx + jj] - Tzz[kk*nz*nx + (ii+1)*nx + jj]) +
                    128625.0f*(Tzz[kk*nz*nx + ii*nx + jj]     - Tzz[kk*nz*nx + (ii-1)*nx + jj]))/(107520.0f*dz);

            rhoz = 0.5f * (rho[kk*nz*nx + (ii+1)*nx + jj] + rho[kk*nz*nx + ii*nx + jj]);

            Vz[index] += dt/rhoz*(d_Txz_dx + d_Tyz_dy + d_Tzz_dz); 
        }
    }    
}

// Elástico Isotrópico 3D com matrizes lineares (GitHub do Komatitsch,2009) - Paulo Bastos, GISIS 27/07/2020
void FDM_8E2T_PML_elastic_Isotropic3D(float *Vx,float *Vy,float *Vz,float *Txx,float *Tyy,float *Tzz,float *Txy,float *Txz,float *Tyz,float *rho,float *M,float *L,int nx,int ny,int nz,float dt,float dx,float dy,float dz) 
{
    float d_Txx_dx, d_Txy_dy, d_Txz_dz, d_Vx_dz, d_Vx_dy, d_Vx_dx;
    float d_Txy_dx, d_Tyy_dy, d_Tyz_dz, d_Vy_dx, d_Vy_dz, d_Vy_dy;
    float d_Txz_dx, d_Tyz_dy, d_Tzz_dz, d_Vz_dx, d_Vz_dy, d_Vz_dz;

    float Vx_x, Vx_y, Vx_z, Vy_x, Vy_y, Vy_z, Vz_x, Vz_y, Vz_z;
    float Txx_x,Txx_y,Txx_z,Tyy_x,Tyy_y,Tyy_z,Tzz_x,Tzz_y,Tzz_z;
    float Txy_x,Txy_y,Txy_z,Txz_x,Txz_y,Txz_z,Tyz_x,Tyz_y,Tyz_z;

    int index, kk, jj, ii;

    float rhox, rhoy, rhoz; // Interpolated density 
    float M_xy, M_xz, M_yz; // Interpolated rigity

    for(index = 0; index < nx*ny*nz; index++) 
    {    
        kk = floor(index/(nx*nz));          // indicador de matrizes (direção y)
        jj = index % nx;                    // indicador de colunas  (direção x)
        ii = floor((index % (nx*nz)) / nx); // indicador de linhas   (direção z)  

        if((ii > 3) && (ii < nz-3) && (kk > 3) && (kk < ny-3) && (jj >= 3) && (jj < nx-4)) 
        {    
            d_Vx_dx = (75.0f*(Vx[kk*nz*nx + ii*nx + (jj-3)] - Vx[kk*nz*nx + ii*nx + (jj+4)]) +
                     1029.0f*(Vx[kk*nz*nx + ii*nx + (jj+3)] - Vx[kk*nz*nx + ii*nx + (jj-2)]) +
                     8575.0f*(Vx[kk*nz*nx + ii*nx + (jj-1)] - Vx[kk*nz*nx + ii*nx + (jj+2)]) +
                   128625.0f*(Vx[kk*nz*nx + ii*nx + (jj+1)] - Vx[kk*nz*nx + ii*nx + jj]))/(107520.0f*dx);

            d_Vy_dy = (75.0f*(Vy[(kk-4)*nz*nx + ii*nx + jj] - Vy[(kk+3)*nz*nx + ii*nx + jj]) +
                     1029.0f*(Vy[(kk+2)*nz*nx + ii*nx + jj] - Vy[(kk-3)*nz*nx + ii*nx + jj]) +
                     8575.0f*(Vy[(kk-2)*nz*nx + ii*nx + jj] - Vy[(kk+1)*nz*nx + ii*nx + jj]) +
                   128625.0f*(Vy[kk*nz*nx + ii*nx + jj]     - Vy[(kk-1)*nz*nx + ii*nx + jj]))/(107520.0f*dy);

            d_Vz_dz = (75.0f*(Vz[kk*nz*nx + (ii-4)*nx + jj] - Vz[kk*nz*nx + (ii+3)*nx + jj]) +
                     1029.0f*(Vz[kk*nz*nx + (ii+2)*nx + jj] - Vz[kk*nz*nx + (ii-3)*nx + jj]) +
                     8575.0f*(Vz[kk*nz*nx + (ii-2)*nx + jj] - Vz[kk*nz*nx + (ii+1)*nx + jj]) +
                   128625.0f*(Vz[kk*nz*nx + ii*nx + jj]     - Vz[kk*nz*nx + (ii-1)*nx + jj]))/(107520.0f*dz);

            Txx_x = dt *(L[index] + 2*M[index]) * d_Vx_dx;
            Txx_y = dt * L[index] * d_Vy_dy;
            Txx_z = dt * L[index] * d_Vz_dz;

            Txx[index] += Txx_x + Txx_y + Txx_z;

            Tyy_x = dt * L[index] * d_Vx_dx;
            Tyy_y = dt *(L[index] + 2*M[index]) * d_Vy_dy;
            Tyy_z = dt * L[index] * d_Vz_dz;

            Tyy[index] += Tyy_x + Tyy_y + Tyy_z;
                
            Tzz_x = dt * L[index] * d_Vx_dx;
            Tzz_y = dt * L[index] * d_Vy_dy;
            Tzz_z = dt *(L[index] + 2*M[index]) * d_Vz_dz;

            Tzz[index] += Tzz_x + Tzz_y + Tzz_z;
        }

        if((ii >= 0) && (ii < nz) && (kk >= 3) && (kk < ny-4) && (jj > 3) && (jj < nx-3)) 
        {
            d_Vx_dy = (75.0f*(Vx[(kk-3)*nz*nx + ii*nx + jj] - Vx[(kk+4)*nz*nx + ii*nx + jj]) +
                     1029.0f*(Vx[(kk+3)*nz*nx + ii*nx + jj] - Vx[(kk-2)*nz*nx + ii*nx + jj]) +
                     8575.0f*(Vx[(kk-1)*nz*nx + ii*nx + jj] - Vx[(kk+2)*nz*nx + ii*nx + jj]) +
                   128625.0f*(Vx[(kk+1)*nz*nx + ii*nx + jj] - Vx[kk*nz*nx + ii*nx + jj]))/(107520.0f*dy);

            d_Vy_dx = (75.0f*(Vy[kk*nz*nx + ii*nx + (jj-4)] - Vy[kk*nz*nx + ii*nx + (jj+3)]) +
                     1029.0f*(Vy[kk*nz*nx + ii*nx + (jj+2)] - Vy[kk*nz*nx + ii*nx + (jj-3)]) +
                     8575.0f*(Vy[kk*nz*nx + ii*nx + (jj-2)] - Vy[kk*nz*nx + ii*nx + (jj+1)]) +
                   128625.0f*(Vy[kk*nz*nx + ii*nx + jj]     - Vy[kk*nz*nx + ii*nx + (jj-1)]))/(107520.0f*dx);

            // M_xy = powf(0.25*(1/M[(kk+1)*nz*nx + ii*nx + (jj+1)] +
            //                   1/M[kk*nz*nx + ii*nx + (jj+1)] + 
            //                   1/M[(kk+1)*nz*nx + ii*nx + jj] +
            //                   1/M[kk*nz*nx + ii*nx + jj]),-1.0f);

            M_xy = M[index];

            Txy_x = dt * M_xy * d_Vy_dx;
            Txy_y = dt * M_xy * d_Vx_dy;
            Txy_z = 0.0f;

            Txy[index] += Txy_x + Txy_y + Txy_z; 
        }

        if((ii >= 3) && (ii < nz-4) && (kk >= 0) && (kk < ny) && (jj > 3) && (jj < nx-3)) 
        {
            d_Vx_dz = (75.0f*(Vx[kk*nz*nx + (ii-3)*nx + jj] - Vx[kk*nz*nx + (ii+4)*nx + jj]) +
                     1029.0f*(Vx[kk*nz*nx + (ii+3)*nx + jj] - Vx[kk*nz*nx + (ii-2)*nx + jj]) +
                     8575.0f*(Vx[kk*nz*nx + (ii-1)*nx + jj] - Vx[kk*nz*nx + (ii+2)*nx + jj]) +
                   128625.0f*(Vx[kk*nz*nx + (ii+1)*nx + jj] - Vx[kk*nz*nx + ii*nx + jj]))/(107520.0f*dz);

            d_Vz_dx = (75.0f*(Vz[kk*nz*nx + ii*nx + (jj-4)] - Vz[kk*nz*nx + ii*nx + (jj+3)]) +
                     1029.0f*(Vz[kk*nz*nx + ii*nx + (jj+2)] - Vz[kk*nz*nx + ii*nx + (jj-3)]) +
                     8575.0f*(Vz[kk*nz*nx + ii*nx + (jj-2)] - Vz[kk*nz*nx + ii*nx + (jj+1)]) +
                   128625.0f*(Vz[kk*nz*nx + ii*nx + jj]     - Vz[kk*nz*nx + ii*nx + (jj-1)]))/(107520.0f*dx);

            // M_xz = powf(0.25*(1/M[kk*nz*nx + (ii+1)*nx + (jj+1)] +
            //                   1/M[kk*nz*nx + ii*nx + (jj+1)] + 
            //                   1/M[kk*nz*nx + (ii+1)*nx + jj] +
            //                   1/M[kk*nz*nx + ii*nx + jj]),-1.0f);

            M_xz = M[index];

            Txz_x = dt * M_xz * d_Vz_dx;
            Txz_y = 0.0f;
            Txz_z = dt * M_xz * d_Vx_dz;

            Txz[index] += Txz_x + Txz_y + Txz_z; 
        }
    
        if((ii >= 3) && (ii < nz-4) && (kk >= 3) && (kk < ny-4) && (jj >= 0) && (jj < nx)) 
        {
            d_Vy_dz = (75.0f*(Vy[kk*nz*nx + (ii-3)*nx + jj] - Vy[kk*nz*nx + (ii+4)*nx + jj]) +
                     1029.0f*(Vy[kk*nz*nx + (ii+3)*nx + jj] - Vy[kk*nz*nx + (ii-2)*nx + jj]) +
                     8575.0f*(Vy[kk*nz*nx + (ii-1)*nx + jj] - Vy[kk*nz*nx + (ii+2)*nx + jj]) +
                   128625.0f*(Vy[kk*nz*nx + (ii+1)*nx + jj] - Vy[kk*nz*nx + ii*nx + jj]))/(107520.0f*dz);

            d_Vz_dy = (75.0f*(Vz[(kk-3)*nz*nx + ii*nx + jj] - Vz[(kk+4)*nz*nx + ii*nx + jj]) +
                     1029.0f*(Vz[(kk+3)*nz*nx + ii*nx + jj] - Vz[(kk-2)*nz*nx + ii*nx + jj]) +
                     8575.0f*(Vz[(kk-1)*nz*nx + ii*nx + jj] - Vz[(kk+2)*nz*nx + ii*nx + jj]) +
                   128625.0f*(Vz[(kk+1)*nz*nx + ii*nx + jj] - Vz[kk*nz*nx + ii*nx + jj]))/(107520.0f*dy);

            // M_yz = powf(0.25*(1/M[(kk+1)*nz*nx + (ii+1)*nx + jj] +
            //                   1/M[(kk+1)*nz*nx + ii*nx + jj] + 
            //                   1/M[kk*nz*nx + (ii+1)*nx + jj] +
            //                   1/M[kk*nz*nx + ii*nx + jj]),-1.0f);

            M_yz = M[index];

            Tyz_x = 0.0f;
            Tyz_y = dt * M_yz * d_Vz_dy;
            Tyz_z = dt * M_yz * d_Vy_dz;

            Tyz[index] += Tyz_x + Tyz_y + Tyz_z; 
        }
    }

    for(index = 0; index < nx*ny*nz; index++) 
    {    
        kk = floor(index/(nx*nz));          // indicador de matrizes (direção y)
        jj = index % nx;                    // indicador de colunas  (direção x)
        ii = floor((index % (nx*nz)) / nx); // indicador de linhas   (direção z)  

        if((ii > 3) && (ii < nz-3) && (kk > 3) && (kk < ny-3) && (jj > 3) && (jj < nx-3)) 
        {
            d_Txx_dx = (75.0f*(Txx[kk*nz*nx + ii*nx + (jj-4)] - Txx[kk*nz*nx + ii*nx + (jj+3)]) +
                      1029.0f*(Txx[kk*nz*nx + ii*nx + (jj+2)] - Txx[kk*nz*nx + ii*nx + (jj-3)]) +
                      8575.0f*(Txx[kk*nz*nx + ii*nx + (jj-2)] - Txx[kk*nz*nx + ii*nx + (jj+1)]) +
                    128625.0f*(Txx[kk*nz*nx + ii*nx + jj]     - Txx[kk*nz*nx + ii*nx + (jj-1)]))/(107520.0f*dx);

            d_Txy_dy = (75.0f*(Txy[(kk-4)*nz*nx + ii*nx + jj] - Txy[(kk+3)*nz*nx + ii*nx + jj]) +
                      1029.0f*(Txy[(kk+2)*nz*nx + ii*nx + jj] - Txy[(kk-3)*nz*nx + ii*nx + jj]) +
                      8575.0f*(Txy[(kk-2)*nz*nx + ii*nx + jj] - Txy[(kk+1)*nz*nx + ii*nx + jj]) +
                    128625.0f*(Txy[kk*nz*nx + ii*nx + jj]     - Txy[(kk-1)*nz*nx + ii*nx + jj]))/(107520.0f*dy);

            d_Txz_dz = (75.0f*(Txz[kk*nz*nx + (ii-4)*nx + jj] - Txz[kk*nz*nx + (ii+3)*nx + jj]) +
                      1029.0f*(Txz[kk*nz*nx + (ii+2)*nx + jj] - Txz[kk*nz*nx + (ii-3)*nx + jj]) +
                      8575.0f*(Txz[kk*nz*nx + (ii-2)*nx + jj] - Txz[kk*nz*nx + (ii+1)*nx + jj]) +
                    128625.0f*(Txz[kk*nz*nx + ii*nx + jj]     - Txz[kk*nz*nx + (ii-1)*nx + jj]))/(107520.0f*dz);

            // rhox = 0.5f * (rho[kk*nz*nx + ii*nx + (jj+1)] + rho[kk*nz*nx + ii*nx + jj]);

            rhox = rho[index];

            Vx_x = dt/rhox * d_Txx_dx;
            Vx_y = dt/rhox * d_Txy_dy;
            Vx_z = dt/rhox * d_Txz_dz;    

            Vx[index] += Vx_x + Vx_y + Vx_z;
        }

        if((ii > 3) && (ii < nz-3) && (kk >= 3) && (kk < ny-4) && (jj >= 3) && (jj < nx-4)) 
        {
            d_Txy_dx = (75.0f*(Txy[kk*nz*nx + ii*nx + (jj-3)] - Txy[kk*nz*nx + ii*nx + (jj+4)]) +
                      1029.0f*(Txy[kk*nz*nx + ii*nx + (jj+3)] - Txy[kk*nz*nx + ii*nx + (jj-2)]) +
                      8575.0f*(Txy[kk*nz*nx + ii*nx + (jj-1)] - Txy[kk*nz*nx + ii*nx + (jj+2)]) +
                    128625.0f*(Txy[kk*nz*nx + ii*nx + (jj+1)] - Txy[kk*nz*nx + ii*nx + jj]))/(107520.0f*dx);

            d_Tyy_dy = (75.0f*(Tyy[(kk-3)*nz*nx + ii*nx + jj] - Tyy[(kk+4)*nz*nx + ii*nx + jj]) +
                      1029.0f*(Tyy[(kk+3)*nz*nx + ii*nx + jj] - Tyy[(kk-2)*nz*nx + ii*nx + jj]) +
                      8575.0f*(Tyy[(kk-1)*nz*nx + ii*nx + jj] - Tyy[(kk+2)*nz*nx + ii*nx + jj]) +
                    128625.0f*(Tyy[(kk+1)*nz*nx + ii*nx + jj] - Tyy[kk*nz*nx + ii*nx + jj]))/(107520.0f*dy);

            d_Tyz_dz = (75.0f*(Tyz[kk*nz*nx + (ii-4)*nx + jj] - Tyz[kk*nz*nx + (ii+3)*nx + jj]) +
                      1029.0f*(Tyz[kk*nz*nx + (ii+2)*nx + jj] - Tyz[kk*nz*nx + (ii-3)*nx + jj]) +
                      8575.0f*(Tyz[kk*nz*nx + (ii-2)*nx + jj] - Tyz[kk*nz*nx + (ii+1)*nx + jj]) +
                    128625.0f*(Tyz[kk*nz*nx + ii*nx + jj]     - Tyz[kk*nz*nx + (ii-1)*nx + jj]))/(107520.0f*dz);

            // rhoy = 0.5f * (rho[(kk+1)*nz*nx + ii*nx + jj] + rho[kk*nz*nx + ii*nx + jj]);

            rhoy = rho[index];

            Vy_x = dt/rhoy * d_Txy_dx;
            Vy_y = dt/rhoy * d_Tyy_dy;
            Vy_z = dt/rhoy * d_Tyz_dz;    

            Vy[index] += Vy_x + Vy_y + Vy_z;
        }    

        if((ii >= 3) && (ii < nz-4) && (kk > 3) && (kk < ny-3) && (jj >= 3) && (jj < nx-4)) 
        {
            d_Txz_dx = (75.0f*(Txz[kk*nz*nx + ii*nx + (jj-3)] - Txz[kk*nz*nx + ii*nx + (jj+4)]) +
                      1029.0f*(Txz[kk*nz*nx + ii*nx + (jj+3)] - Txz[kk*nz*nx + ii*nx + (jj-2)]) +
                      8575.0f*(Txz[kk*nz*nx + ii*nx + (jj-1)] - Txz[kk*nz*nx + ii*nx + (jj+2)]) +
                    128625.0f*(Txz[kk*nz*nx + ii*nx + (jj+1)] - Txz[kk*nz*nx + ii*nx + jj]))/(107520.0f*dx);

            d_Tyz_dy = (75.0f*(Tyz[(kk-4)*nz*nx + ii*nx + jj] - Tyz[(kk+3)*nz*nx + ii*nx + jj]) +
                      1029.0f*(Tyz[(kk+2)*nz*nx + ii*nx + jj] - Tyz[(kk-3)*nz*nx + ii*nx + jj]) +
                      8575.0f*(Tyz[(kk-2)*nz*nx + ii*nx + jj] - Tyz[(kk+1)*nz*nx + ii*nx + jj]) +
                    128625.0f*(Tyz[kk*nz*nx + ii*nx + jj]     - Tyz[(kk-1)*nz*nx + ii*nx + jj]))/(107520.0f*dy);

            d_Tzz_dz = (75.0f*(Tzz[kk*nz*nx + (ii-3)*nx + jj] - Tzz[kk*nz*nx + (ii+4)*nx + jj]) +
                      1029.0f*(Tzz[kk*nz*nx + (ii+3)*nx + jj] - Tzz[kk*nz*nx + (ii-2)*nx + jj]) +
                      8575.0f*(Tzz[kk*nz*nx + (ii-1)*nx + jj] - Tzz[kk*nz*nx + (ii+2)*nx + jj]) +
                    128625.0f*(Tzz[kk*nz*nx + (ii+1)*nx + jj] - Tzz[kk*nz*nx + ii*nx + jj]))/(107520.0f*dz);

            // rhoz = 0.5f * (rho[kk*nz*nx + (ii+1)*nx + jj] + rho[kk*nz*nx + ii*nx + jj]);

            rhoz = rho[index];

            Vz_x = dt/rhoz * d_Txz_dx;
            Vz_y = dt/rhoz * d_Tyz_dy;
            Vz_z = dt/rhoz * d_Tzz_dz;    

            Vz[index] += Vz_x + Vz_y + Vz_z;        
        }
    }
}

//
void cerjanCorners3D(float * prismX,float * prismY,float * prismZ,float * cube,float *factor,int n_abc,int nx,int ny,int nz) 
{
    int ii,jj,kk,index;
    float * slice = (float *) malloc(n_abc*n_abc*sizeof(float));

    /* Upper left absorbing corner */
    for(ii = n_abc-1; ii >= 0; ii--) 
    {
        slice[ii*n_abc + ii] = factor[ii];    
        for(kk = 1; kk <= n_abc; kk++) 
        {
            slice[(ii+kk)*n_abc + ii] = factor[ii];
            slice[ii*n_abc + ii + kk] = factor[ii];
        }    
    }

    /* Creating cube Upper left */
    for(ii = 0; ii < n_abc; ii++)
    {
        for(jj = 0; jj < n_abc; jj++)
        {
            cube[ii*n_abc*n_abc + jj*n_abc + jj] = slice[ii*n_abc + jj];
        }

        for(jj = 0; jj < n_abc-1; jj++)
        {
            for(kk = jj+1; kk < n_abc; kk++)
            {
                cube[ii*n_abc*n_abc + jj*n_abc + kk] = cube[ii*n_abc*n_abc + jj*n_abc + jj];
                cube[ii*n_abc*n_abc + kk*n_abc + jj] = cube[ii*n_abc*n_abc + jj*n_abc + jj];
            }
        }
    }

    /* Setting prism X Upper left */
    for (jj = 0; jj < nx; jj++)
    {
        for(index = 0; index < n_abc*n_abc; index++)
        {
            ii = index / n_abc; // z projection
            kk = index % n_abc; // y projection

            prismX[jj*n_abc*n_abc + ii*n_abc + kk] = slice[ii*n_abc + kk];
        }
    }

    /* Setting prism Y Upper left */
    for (kk = 0; kk < ny; kk++)
    {
        for(index = 0; index < n_abc*n_abc; index++)
        {
            ii = index / n_abc; // z projection
            jj = index % n_abc; // x projection

            prismY[kk*n_abc*n_abc + ii*n_abc + jj] = slice[ii*n_abc + jj];
        }
    }

    /* Setting prism Z Upper left */
    for (ii = 0; ii < nz; ii++)
    {
        for(index = 0; index < n_abc*n_abc; index++)
        {
            kk = index / n_abc; // y projection
            jj = index % n_abc; // x projection

            prismZ[ii*n_abc*n_abc + kk*n_abc + jj] = slice[kk*n_abc + jj];
        }
    }
}

//
void cerjanAcousticAbsorbingCondition3D(float * P_pre,float * P_fut,float *factor,float * prismX,float * prismY,float * prismZ,float * cubXYZ,int nxx,int nyy,int nzz,int n_abc)
{
    int index,kk,jj,ii;

    for(index = 0; index < nxx*nyy*nzz; index++)   /* Spatial loop */ 
    {
        kk = floor(index/(nxx*nzz));               // indicador de matrizes (direção y)
        jj = index % nxx;                          // indicador de colunas  (direção x)
        ii = floor((index % (nxx*nzz)) / nxx);     // indicador de linhas   (direção z)  

        // Absorvendo na direção X
        if((ii >= n_abc) && (ii < nzz-n_abc) && (jj >= 0) && (jj < n_abc) && (kk >= n_abc) && (kk < nyy-n_abc)) 
        {
            P_pre[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];
            P_fut[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];

            P_pre[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
            P_fut[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
        }

        // Absorvendo na direção Y
        if((ii >= n_abc) && (ii < nzz-n_abc) && (jj >= n_abc) && (jj < nxx-n_abc) && (kk >= 0) && (kk < n_abc)) 
        {
            P_pre[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];
            P_fut[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];

            P_pre[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
            P_fut[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
        }

        // Absorvendo na direção Z  
        if((ii >= 0) && (ii < n_abc) && (jj >= n_abc) && (jj < nxx-n_abc) && (kk >= n_abc) && (kk < nyy-n_abc)) 
        {
            P_pre[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];
            P_fut[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];

            // P_pre[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
            // P_fut[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
        }

        // Quinas no dimínio X
        if((ii >= 0) && (ii < n_abc) && (jj >= n_abc) && (jj < nxx-n_abc) && (kk >= 0) && (kk < n_abc)) 
        {
            P_pre[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            P_fut[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];

            P_pre[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            P_fut[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];

            P_pre[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            P_fut[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];

            P_pre[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            P_fut[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
        }        

        // Quinas no domínio Y
        if((ii >= 0) && (ii < n_abc) && (jj >= 0) && (jj < n_abc) && (kk >= n_abc) && (kk < nyy-n_abc)) 
        {
            P_pre[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            P_fut[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];

            P_pre[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            P_fut[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];

            P_pre[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            P_fut[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];

            P_pre[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            P_fut[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
        }        

        // Quinas no domínio Z
        if((ii >= n_abc) && (ii < nzz-n_abc) && (jj >= 0) && (jj < n_abc) && (kk >= 0) && (kk < n_abc)) 
        {
            P_pre[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            P_fut[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];

            P_pre[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            P_fut[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];

            P_pre[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            P_fut[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];

            P_pre[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            P_fut[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
        }        
    
        // Quinas tridimensionais
        if((ii >= 0) && (ii < n_abc) && (jj >= 0) && (jj < n_abc) && (kk >= 0) && (kk < n_abc))
        {
            P_pre[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            P_fut[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            P_pre[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            P_fut[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
             
            P_pre[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            P_fut[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            P_pre[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            P_fut[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            P_pre[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            P_fut[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            P_pre[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            P_fut[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
             
            P_pre[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            P_fut[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            P_pre[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            P_fut[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

        }
    }
}

//
void cerjanPVAcousticAbsorbingCondition3D(float * ps,float * vx,float * vy,float * vz,float *factor,float * prismX,float * prismY,float * prismZ,float * cubXYZ,int nxx,int nyy,int nzz,int n_abc)
{
    int index,kk,jj,ii;

    for(index = 0; index < nxx*nyy*nzz; index++)   /* Spatial loop */ 
    {
        kk = floor(index/(nxx*nzz));               // indicador de matrizes (direção y)
        jj = index % nxx;                          // indicador de colunas  (direção x)
        ii = floor((index % (nxx*nzz)) / nxx);     // indicador de linhas   (direção z)  

        // Absorvendo na direção X
        if((ii >= n_abc) && (ii < nzz-n_abc) && (jj >= 0) && (jj < n_abc) && (kk >= n_abc) && (kk < nyy-n_abc)) 
        {
            ps[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];
            vx[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];
            vy[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];
            vz[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];

            ps[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
            vx[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
            vy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
            vz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
        }

        // Absorvendo na direção Y
        if((ii >= n_abc) && (ii < nzz-n_abc) && (jj >= n_abc) && (jj < nxx-n_abc) && (kk >= 0) && (kk < n_abc)) 
        {
            ps[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];
            vx[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];
            vy[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];
            vz[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];

            ps[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
            vx[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
            vy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
            vz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
        }

        // Absorvendo na direção Z  
        if((ii >= 0) && (ii < n_abc) && (jj >= n_abc) && (jj < nxx-n_abc) && (kk >= n_abc) && (kk < nyy-n_abc)) 
        {
            ps[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];
            vx[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];
            vy[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];
            vz[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];

            ps[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
            vx[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
            vy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
            vz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
        }

        // Quinas no dimínio X
        if((ii >= 0) && (ii < n_abc) && (jj >= n_abc) && (jj < nxx-n_abc) && (kk >= 0) && (kk < n_abc)) 
        {
            ps[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            vx[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            vy[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            vz[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];

            ps[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            vx[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            vy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            vz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];

            ps[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            vx[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            vy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            vz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];

            ps[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            vx[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            vy[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            vz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
        }        

        // Quinas no domínio Y
        if((ii >= 0) && (ii < n_abc) && (jj >= 0) && (jj < n_abc) && (kk >= n_abc) && (kk < nyy-n_abc)) 
        {
            ps[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            vx[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            vy[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            vz[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];

            ps[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            vx[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            vy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            vz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];

            ps[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            vx[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            vy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            vz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];

            ps[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            vx[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            vy[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            vz[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
        }        

        // Quinas no domínio Z
        if((ii >= n_abc) && (ii < nzz-n_abc) && (jj >= 0) && (jj < n_abc) && (kk >= 0) && (kk < n_abc)) 
        {
            ps[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            vx[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            vy[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            vz[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];

            ps[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            vx[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            vy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            vz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];

            ps[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            vx[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            vy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            vz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];

            ps[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            vx[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            vy[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            vz[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
        }        
    
        // Quinas tridimensionais
        if((ii >= 0) && (ii < n_abc) && (jj >= 0) && (jj < n_abc) && (kk >= 0) && (kk < n_abc))
        {
            ps[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vx[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vy[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vz[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            ps[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vx[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
             
            ps[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vx[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            ps[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vx[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vy[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vz[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            ps[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vx[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            ps[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vx[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vy[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
             
            ps[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vx[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vy[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vz[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            ps[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vx[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vy[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            vz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
        }
    }
}

//
void cerjanElasticAbsorbingCondition3D(float *Vx,float *Vy,float *Vz,float *Txx,float *Tyy,float *Tzz,float *Txy,float *Txz,float *Tyz,float *factor,float *prismX,float *prismY,float *prismZ,float *cubXYZ,int nxx,int nyy,int nzz,int n_abc)
{
    int index,kk,jj,ii;

    for(index = 0; index < nxx*nyy*nzz; index++)   /* Spatial loop */ 
    {
        kk = floor(index/(nxx*nzz));               // indicador de matrizes (direção y)
        jj = index % nxx;                          // indicador de colunas  (direção x)
        ii = floor((index % (nxx*nzz)) / nxx);     // indicador de linhas   (direção z)  

        // Absorvendo na direção X
        if((ii >= n_abc) && (ii < nzz-n_abc) && (jj >= 0) && (jj < n_abc) && (kk >= n_abc) && (kk < nyy-n_abc)) 
        {
            Vx[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];
            Vy[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];
            Vz[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];
            Txx[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];
            Txy[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];
            Txz[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];
            Tyy[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];
            Tyz[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];
            Tzz[kk*nxx*nzz + ii*nxx + jj] *= factor[jj];

            Vx[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
            Vy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
            Vz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
            Txx[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
            Txy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
            Txz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
            Tyy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
            Tyz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
            Tzz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= factor[jj];
        }

        // Absorvendo na direção Y
        if((ii >= n_abc) && (ii < nzz-n_abc) && (jj >= n_abc) && (jj < nxx-n_abc) && (kk >= 0) && (kk < n_abc)) 
        {
            Vx[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Vy[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Vz[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Txx[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Txy[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Txz[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Tyy[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Tyz[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Tzz[kk*nxx*nzz + ii*nxx + jj] *= factor[kk];

            Vx[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Vy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Vz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Txx[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Txy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Txz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Tyy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Tyz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
            Tzz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= factor[kk];
        }

        // Absorvendo na direção Z  
        if((ii >= 0) && (ii < n_abc) && (jj >= n_abc) && (jj < nxx-n_abc) && (kk >= n_abc) && (kk < nyy-n_abc)) 
        {
            Vx[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];
            Vy[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];
            Vz[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];
            Txx[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];
            Txy[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];
            Txz[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];
            Tyy[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];
            Tyz[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];
            Tzz[kk*nxx*nzz + ii*nxx + jj] *= factor[ii];

            Vx[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
            Vy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
            Vz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
            Txx[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
            Txy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
            Txz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
            Tyy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
            Tyz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
            Tzz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= factor[ii];
        }

        // Quinas no dimínio X
        if((ii >= 0) && (ii < n_abc) && (jj >= n_abc) && (jj < nxx-n_abc) && (kk >= 0) && (kk < n_abc)) 
        {
            Vx[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Vy[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Vz[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Txx[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Txy[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Txz[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Tyy[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Tyz[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Tzz[kk*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];

            Vx[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Vy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Vz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Txx[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Txy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Txz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Tyy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Tyz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Tzz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];

            Vx[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Vy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Vz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Txx[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Txy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Txz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Tyy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Tyz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Tzz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];

            Vx[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Vy[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Vz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Txx[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Txy[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Txz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Tyy[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Tyz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
            Tzz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismX[(jj-n_abc)*n_abc*n_abc + ii*n_abc + kk];
        }        

        // Quinas no domínio Y
        if((ii >= 0) && (ii < n_abc) && (jj >= 0) && (jj < n_abc) && (kk >= n_abc) && (kk < nyy-n_abc)) 
        {
            Vx[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Vy[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Vz[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Txx[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Txy[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Txz[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Tyy[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Tyz[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Tzz[kk*nxx*nzz + ii*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];

            Vx[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Vy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Vz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Txx[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Txy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Txz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Tyy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Tyz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Tzz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];

            Vx[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Vy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Vz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Txx[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Txy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Txz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Tyy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Tyz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Tzz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];

            Vx[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Vy[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Vz[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Txx[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Txy[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Txz[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Tyy[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Tyz[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
            Tzz[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= prismY[(kk-n_abc)*n_abc*n_abc + ii*n_abc + jj];
        }        

        // Quinas no domínio Z
        if((ii >= n_abc) && (ii < nzz-n_abc) && (jj >= 0) && (jj < n_abc) && (kk >= 0) && (kk < n_abc)) 
        {
            Vx[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Vy[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Vz[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Txx[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Txy[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Txz[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Tyy[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Tyz[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Tzz[kk*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];

            Vx[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Vy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Vz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Txx[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Txy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Txz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Tyy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Tyz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Tzz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];

            Vx[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Vy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Vz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Txx[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Txy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Txz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Tyy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Tyz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Tzz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];

            Vx[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Vy[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Vz[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Txx[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Txy[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Txz[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Tyy[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Tyz[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
            Tzz[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= prismZ[(ii-n_abc)*n_abc*n_abc + kk*n_abc + jj];
        }        
    
        // Quinas tridimensionais
        if((ii >= 0) && (ii < n_abc) && (jj >= 0) && (jj < n_abc) && (kk >= 0) && (kk < n_abc))
        {
            Vx[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vy[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vz[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txx[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txy[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txz[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyy[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyz[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tzz[kk*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            Vx[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txx[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyy[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tzz[(nyy-kk-1)*nxx*nzz + ii*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
             
            Vx[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txx[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyy[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tzz[kk*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            Vx[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vy[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vz[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txx[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txy[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txz[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyy[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyz[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tzz[(nyy-kk-1)*nxx*nzz + ii*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            Vx[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txx[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyy[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tzz[kk*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            Vx[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vy[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txx[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txy[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyy[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tzz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + jj] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
             
            Vx[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vy[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vz[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txx[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txy[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txz[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyy[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyz[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tzz[kk*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];

            Vx[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vy[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Vz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txx[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txy[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Txz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyy[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tyz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
            Tzz[(nyy-kk-1)*nxx*nzz + (nzz-ii-1)*nxx + (nxx-jj-1)] *= cubXYZ[kk*n_abc*n_abc + ii*n_abc + jj];
        }
    }
}

//
void waveFieldUpdate(float * pas, float * pre, float * fut, int nPoints)
{
    int index;

    #pragma acc parallel loop present(pas[0:nPoints],pre[0:nPoints],fut[0:nPoints])
    for(index = 0; index < nPoints; ++index)
    {
        pas[index] = pre[index];         
        pre[index] = fut[index];         
    }
}

//
void getSquareSumField(float * waveField, float * waveFieldSum, int nxx, int nzz, int abc)
{
    int ii, jj, index;

    #pragma acc parallel loop present(waveField[0:nxx*nzz],waveFieldSum[0:nx*nz])
    for(index = 0; index < nxx*nzz; ++index)
    {
        ii = floor(index / nxx);  /* line indicator */
        jj = index % nxx;         /* column indicator */                              

        /* Collect the sum of wavefield to applies compensations to final image */
        if((ii >= abc) && (ii < nzz - abc) && (jj >= abc) && (jj < nxx - abc)) 
        {                            
            waveFieldSum[(ii-abc)*(nxx - 2*abc)  + (jj-abc)] += powf(waveField[ii*nxx + jj],2.0f);
        }
    }
}

//
void getWaveField(int kk,float * waveField,float * waveFieldToGet,int nxx,int nzz,int nt,int abc,int sampleInterval)
{
    int ii, jj, index;

    if(kk % sampleInterval == 0)
    {   
        #pragma acc parallel loop present(wvefield[0:nxx*nzz],wavefieldToGet[(nxx-2*abc)*(nzz-2*abc)*(nt/10)])                  
        for(index = 0; index < nxx*nzz; ++index)
        {
            ii = floor(index / nxx);  /* line indicator */
            jj = index % nxx;         /* column indicator */  

            /* Getting wave field to use in cross correlation */
            if((ii >= abc) && (ii < nzz - abc) && (jj >= abc) && (jj < nxx - abc)) 
            {
                waveFieldToGet[(kk/sampleInterval)*(nzz-2*abc)*(nxx-2*abc) + (ii-abc)*(nxx-2*abc) + (jj-abc)] = waveField[ii*nxx + jj];
            }
        }
    }
}

//
void crossCorrelation(int kk, float * image, float * directField, float * reverseField, int nxx, int nzz, int abc, int sampleInterval)
{
    int ii, jj, index;

    if(kk % sampleInterval == 0)
    {
        #pragma acc parallel loop present(image[0:nx*nz],directField[nx*nz*(nt/10)],reverseField[0:nxx*nzz])
        for(index = 0; index < nxx*nzz; ++index)
        {
            ii = floor(index / nxx);  /* line indicator */
            jj = index % nxx;         /* column indicator */                              

            /* Cross Correlation image condition */
            if((ii >= abc) && (ii < nzz - abc) && (jj >= abc) && (jj < nxx - abc)) 
            {                            
                image[(ii-abc)*(nxx-2*abc) + (jj-abc)] += directField[(kk/sampleInterval)*(nxx-2*abc)*(nzz-2*abc) + (ii-abc)*(nxx-2*abc) + (jj-abc)] * reverseField[ii*nxx + jj];
            }
        }
        
        // export_float32("formacao.bin",nx*nz,image);
    }   
}

//
float * DsumCompensation(float * image, float * directFieldSum, int nxx, int nzz, int abc)
{
    int ii, jj, index;
    int nx = nxx - 2*abc;
    int nz = nzz - 2*abc;
    
    float * imageCompensed = (float *) malloc(nx*nz*sizeof(float));

    for(index = 0; index < nx*nz; ++index)
    {
        ii = floor(index / nx);  /* line indicator */
        jj = index % nx;         /* column indicator */                              

        imageCompensed[ii*nx + jj] = image[ii*nx + jj] / directFieldSum[ii*nx + jj];
    }

    return imageCompensed;
}

//
float * RsumCompensation(float * image, float * reverseFieldSum, int nxx, int nzz, int abc)
{
    int ii, jj, index;
    int nx = nxx - 2*abc;
    int nz = nzz - 2*abc;
    
    float * imageCompensed = (float *) malloc(nx*nz*sizeof(float));

    for(index = 0; index < nx*nz; ++index)
    {
        ii = floor(index / nx);  /* line indicator */
        jj = index % nx;         /* column indicator */                              

        imageCompensed[ii*nx + jj] = image[ii*nx + jj] / reverseFieldSum[ii*nx + jj];
    }

    return imageCompensed;
}

//
float * sqrtDsumCompensation(float * image, float * directFieldSum, int nxx, int nzz, int abc)
{
    int ii, jj, index;
    int nx = nxx - 2*abc;
    int nz = nzz - 2*abc;
    
    float * imageCompensed = (float *) malloc(nx*nz*sizeof(float));

    for(index = 0; index < nx*nz; ++index)
    {
        ii = floor(index / nx);  /* line indicator */
        jj = index % nx;         /* column indicator */                              

        imageCompensed[ii*nx + jj] = image[ii*nx + jj] /  sqrt(directFieldSum[ii*nx + jj]);
    }

    return imageCompensed;
}

//
float * sqrtRsumCompensation(float * image, float * reverseFieldSum, int nxx, int nzz, int abc)
{
    int ii, jj, index;
    int nx = nxx - 2*abc;
    int nz = nzz - 2*abc;
    
    float * imageCompensed = (float *) malloc(nx*nz*sizeof(float));

    for(index = 0; index < nx*nz; ++index)
    {
        ii = floor(index / nx);  /* line indicator */
        jj = index % nx;         /* column indicator */                              

        imageCompensed[ii*nx + jj] = image[ii*nx + jj] / sqrt(reverseFieldSum[ii*nx + jj]);
    }

    return imageCompensed;
}

//
float * DRsumCompensation(float * image, float * directFieldSum,float * reverseFieldSum, int nxx, int nzz, int abc)
{
    int ii, jj, index;
    int nx = nxx - 2*abc;
    int nz = nzz - 2*abc;
    
    float * imageCompensed = (float *) malloc(nx*nz*sizeof(float));

    for(index = 0; index < nx*nz; ++index)
    {
        ii = floor(index / nx);  /* line indicator */
        jj = index % nx;         /* column indicator */                              

        imageCompensed[ii*nx + jj] = image[ii*nx + jj] / (directFieldSum[ii*nx + jj] * reverseFieldSum[ii*nx + jj]);
    }

    return imageCompensed;
}

//
float * laplaciano(float *image, int nx, int nz, float dh)
{
    int ii,jj,index;
    float *filtr_image = (float *) malloc(nx*nz*sizeof(float));
    float d2_U_dx2,d2_U_dz2;

    for(index = 0; index < nx*nz; index++)
    {
        ii = floor(index / nx);
        jj = index % nx;

        // d2_U_dx2 = (image[(ii-1)*nx + jj] - 2.0f*image[ii*nx + jj] + image[(ii+1)*nx + jj])/powf(dh,2.0f);
        d2_U_dx2 = 0.0f;

        d2_U_dz2 = (image[ii*nx + (jj-1)] - 2.0f*image[nx*ii + jj] + image[(ii+1)*nx + jj])/powf(dh,2.0f);

        filtr_image[ii*nx + jj] = d2_U_dx2 + d2_U_dz2; 
    }    

    return filtr_image;
}

# endif
