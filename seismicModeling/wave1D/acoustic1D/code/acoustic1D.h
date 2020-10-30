# ifndef ACOUSTIC1D_H
# define ACOUSTIC1D_H

void importVector(float * vector, int nPoints, char filename[])
{
    FILE * read = fopen((const char *) filename,"rb");
    fread(vector,sizeof(float),nPoints,read);
    fclose(read);
}

void exportVector(float * vector, int nPoints, char filename[])
{
    FILE * write = fopen((const char *) filename,"wb");
    fwrite(vector,sizeof(float),nPoints,write);
    fclose(write);
}

void readParameters(int * nz, int * nt, float * dz, float * dt, int * nabc, int * nsrc, char filename[])
{
    FILE * arq = fopen((const char *) filename,"r"); 
    if(arq != NULL) 
    {
        fscanf(arq,"%i",nz); fscanf(arq,"%i",nt); 
        fscanf(arq,"%f",dz); fscanf(arq,"%f",dt); 
        fscanf(arq,"%i",nabc); fscanf(arq,"%i",nsrc); 
    } 
    fclose(arq);
}

void setToZero(float * P_pas, float * P_pre, float * P_fut, int nPoints)
{
    for (int index = 0; index < nPoints; index++)
    {
        P_pas[index] = 0.0f;
        P_pre[index] = 0.0f;
        P_fut[index] = 0.0f;
    }
}

void updateWavefield(float * P_pas, float * P_pre, float * P_fut, int nPoints)
{
    for(int index = 0; index < nPoints; index++) 
    {
        P_pas[index] = P_pre[index]; 
        P_pre[index] = P_fut[index];    
    }
}

void acoustic1D_2E2T(float * P_pas, float * P_pre, float * P_fut, float * vp, float * damp, float * source, int timePointer, int nsrc, int nabc, int nz, float dz, float dt)
{
    for (int index = 0; index < nz; index ++)
    {
        if ((index == 0) && (timePointer < nsrc))
        {
            P_pre[nabc] += source[index] / dz;
        }    

        if ((index > 0) && (index < nz-1))
        {        
            float d2_P_dx2 = (P_pre[index+1] - 2.0f*P_pre[index] + P_pre[index-1])/powf(dz,2.0f); 
            float aux1 = powf(dt,2.0f)*powf(vp[index],2.0f);                  
            float aux2 = 2.0f*P_pre[index] - P_pas[index];

            P_fut[index] = (aux1*d2_P_dx2 + aux2) * damp[index] ;  
        }
    }
}

void acoustic1D_4E2T(float * P_pas, float * P_pre, float * P_fut, float * vp, float * damp, float * source, int timePointer, int nsrc, int nabc, int nz, float dz, float dt)
{
    for (int index = 0; index < nz; index ++)
    {
        if ((index == 0) && (timePointer < nsrc))
        {
            P_pre[nabc] += source[index] / dz;
        }    

        if ((index > 1) && (index < nz-2))
        {        
            float d2_P_dx2 = (-P_pre[index-2] + 16.0f*P_pre[index-1] - 30.0f*P_pre[index] + 16.0f*P_pre[index+1] - P_pre[index+2])/(12.0f*powf(dz,2.0f)); 
            float aux1 = powf(dt,2.0f)*powf(vp[index],2.0f);                  
            float aux2 = 2.0f*P_pre[index] - P_pas[index];                       

            P_fut[index] = (aux1*d2_P_dx2 + aux2) * damp[index];  
        }
    }
}

void acoustic1D_8E2T(float * P_pas, float * P_pre, float * P_fut, float * vp, float * damp, float * source, int timePointer, int nsrc, int nabc, int nz, float dz, float dt)
{
    for (int index = 0; index < nz; index ++)
    {
        if ((index == 0) && (timePointer < nsrc))
        {
            P_pre[nabc] += source[timePointer] / dz;
        }    

        if ((index > 3) && (index < nz-4))
        {        
            float d2_P_dx2 = (- 9.0f*(P_pre[index-4] + P_pre[index+4])
                          +   128.0f*(P_pre[index-3] + P_pre[index+3])
                          -  1008.0f*(P_pre[index-2] + P_pre[index+2])
                          +  8064.0f*(P_pre[index+1] + P_pre[index-1])
                          - 14350.0f*P_pre[index])/(5040.0f*powf(dz,2.0f));                       

            float aux1 = powf(dt,2.0f)*powf(vp[index],2.0f);                  
            float aux2 = 2.0f*P_pre[index] - P_pas[index];                       

            P_fut[index] = (aux1*d2_P_dx2 + aux2) * damp[index];  
        }
    }
}

void getSnapshots1D(float * snap, float * P_fut, int nzz, int nabc, int timePointer, int totalTime, int nSnap)
{
    if(timePointer % (totalTime / nSnap) == 0)
    {
        for (int ii = nabc; ii < nzz-nabc; ii++)
        {
            snap[(timePointer / (totalTime / nSnap))*(nzz-2*nabc) + (ii-nabc)] = P_fut[ii];
        }
    }
}

void getSeismogram1D(float * seism, float * P_fut, int nabc, int timePointer)
{
    seism[timePointer] = P_fut[nabc];
}

# endif