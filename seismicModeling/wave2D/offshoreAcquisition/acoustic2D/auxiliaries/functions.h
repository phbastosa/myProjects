# ifndef FUNCTIONS_H_DEFINED
# define FUNCTIONS_H_DEFINED

void setWaveField(float * p1, float * p2, float * p3, int n)
{
    for(int i = 0; i < n; i++)
    {
        p1[i] = 0.0f; p2[i] = 0.0f; p3[i] = 0.0f;
    } 
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

void exportVector(float * vector, int nPoints, char filename[])
{
    FILE * write = fopen((const char *) filename, "wb");
    fwrite(vector, sizeof(float), nPoints, write);
    fclose(write);
}

void readParameters(int *nx, int *nz, int *nt, float *dx, float *dz, float *dt, int *abc, int *nrec, int *nShots, int *nsrc, int *spread, char filename[])
{
    FILE * arq = fopen((const char *) filename,"r"); 
    if(arq != NULL) 
    {
        fscanf(arq,"%i",nx); fscanf(arq,"%i",nz); fscanf(arq,"%i",nt); 
        fscanf(arq,"%f",dx); fscanf(arq,"%f",dz); fscanf(arq,"%f",dt); 
        fscanf(arq,"%i",abc); fscanf(arq,"%i",nrec); fscanf(arq,"%i",nShots); 
        fscanf(arq,"%i",nsrc); fscanf(arq,"%i",spread); 
    } 
    fclose(arq);
}

float * getVelocities(int nxx, int nzz, float * vp)
{
    float v_max = vp[0];
    float v_min = vp[0];
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

void printStatus(float * vels, float dx, float dz, float dt, int nxx, int nzz, int abc, int nt)
{
    printf("Synthetic seismic modeling\n\n");
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

void modelingStatus(int shot, int time, int * xsrc, int n_shot, int * xrec, int spread, float dx, 
                    float dz, int nt, float * vels, float dt, int nxx, int nzz, int abc)
{
    if (time % (nt / 1000) == 0)
    {
        system("clear");        
        printStatus(vels,dx,dz,dt,nxx,nzz,abc,nt);    

        printf("------------------------------------------------------\n");
        printf("    Propagation progress: %.2f %%\n",(float) time/nt * 100.0f);
        printf("------------------------------------------------------\n\n");

        printf("Modeling status:\n");
        printf("   Shot position: %.1f meters\n",(xsrc[shot]-abc)*dx);
        printf("   Recivers position: %.1f - %.1f meters\n",(xrec[shot]-abc)*dx,(xrec[shot+spread-1]-abc)*dx);
        printf("   Total progress: %.2f %%\n",(float) shot/n_shot * 100.0f);
        printf("\nExported seismograms: %i of %i\n",shot,n_shot);
    }
}   

void ajustCoordinates(int *xrec, int *xsrc, int *topo, int nabc, int nxx, int nrec, int nsrc)
{
    for (int ii = 0; ii < nrec; ii++) xrec[ii] += nabc;
    for (int ii = 0; ii < nsrc; ii++) xsrc[ii] += nabc;    

    for (int ii = 0; ii < nxx; ii++) topo[ii] = nabc + 5;    
}

void FDM_8E2T_acoustic2D(int shot, int time, float * vp, float * P_pre, float * P_pas, float * P_fut,
                        float * source, int nsrc, int * topo, int * x_src, int nxx, int nzz, float dx, float dz, float dt)
{
    for(int index = 0; index < nxx*nzz; ++index) 
    {
        int ii = floor(index / nxx);  /* Line indicator */
        int jj = index % nxx;         /* Column indicator */  
        
        if((time < nsrc) && (index == 0))
        { 
            P_pre[topo[x_src[shot]]*nxx + x_src[shot]] = source[time] / (dx*dz); 
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

            P_fut[ii*nxx + jj] = powf(dt,2.0f)*powf(vp[ii*nxx + jj],2.0f)
            * (d2_Px2 + d2_Pz2) + 2.0f*P_pre[ii*nxx + jj] - P_pas[ii*nxx + jj];
        }
    }
}

void cerjanAbsorbingBoundaryCondition(float * P_pas,float * P_pre,float * P_fut,float * damp, int nxx, int nzz)
{
    for (int index = 0; index < nxx*nzz; index++)
    {
        P_pas[index] *= damp[index];
        P_pre[index] *= damp[index];
        P_fut[index] *= damp[index];
    }
}

void waveFieldUpdate(float * pas, float * pre, float * fut, int nPoints)
{
    for(int index = 0; index < nPoints; ++index)
    {
        pas[index] = pre[index];         
        pre[index] = fut[index];         
    }
}

void getSeismograms(float * seism, float * P_pre, int * xrec, int * topo, int nrec, int nxx, int shot, int time)
{
    for (int ii = 0; ii < nrec; ii++)
    {
        seism[time*nrec + ii] = P_pre[topo[xrec[shot + ii]]*nxx + xrec[shot + ii]];
    }
}

void getSnapshots(FILE * snap, float * snapshot, float * P_pre, float * vp, int nxx, int nzz, int nabc, int time, int nt, int nsnap, float parVel)
{
    if (time % (nt / nsnap) == 0)
    {
        for (int index = 0; index < nxx*nzz; index++)
        {
            int ii = floor(index / nxx);  
            int jj = index % nxx;

            if ((ii >= nabc) && (ii < nzz - nabc) && (jj >= nabc) && (jj < nxx - nabc))
            {
                snapshot[(ii-nabc)*(nxx-2*nabc) + (jj-nabc)] = P_pre[ii*nxx + jj] + parVel*vp[ii*nxx + jj];           
            }    
        }
        fwrite(snapshot,sizeof(float),(nxx-2*nabc)*(nzz-2*nabc),snap);
    }            
}

# endif