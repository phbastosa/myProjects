# ifndef ACOUSTIC1D_H
# define ACOUSTIC1D_H

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

void acoustic1D_2E2T(float * P_pas, float * P_pre,float * P_fut, float * vp, int nx, float dx, float dt)
{
    for(int ii = 1; ii < nx-1; ii++) 
    {    

        
        float d2_P_dx2 = (P_pre[ii+1] - 2.0f*P_pre[ii] + P_pre[ii-1])/powf(dx,2.0f); 
        float aux1 = powf(dt,2.0f)*powf(vp[ii],2.0f);                  
        float aux2 = 2.0f*P_pre[ii] - P_pas[ii];

        P_fut[ii] = aux1*d2_P_dx2 + aux2;  
    }
}

void acoustic1D_4E2T(float * P_pas, float * P_pre, float * P_fut, float * vp, int nx, float dx, float dt)
{
    for(int ii = 2; ii < nx-2; ii++) 
    {
        float d2_P_dx2 = (-P_pre[ii-2] + 16.0f*P_pre[ii-1] - 30.0f*P_pre[ii] + 16.0f*P_pre[ii + 1] - P_pre[ii + 2])/(12.0f*powf(dx,2.0f)); 
        float aux1 = powf(dt,2.0f)*powf(vp[ii],2.0f);                  
        float aux2 = 2.0f*P_pre[ii] - P_pas[ii];                       

        P_fut[ii] = aux1*d2_P_dx2 + aux2;  
    }
}

void acoustic1D_8E2T(float * P_pas, float * P_pre, float * P_fut, float * vp, int nx, float dx, float dt)
{
    for(int ii = 4; ii < nx-4; ii++) 
    {
        float d2_P_dx2 = (- 9.0f*(P_pre[ii-4] + P_pre[ii+4])
                      +   128.0f*(P_pre[ii-3] + P_pre[ii+3])
                      -  1008.0f*(P_pre[ii-2] + P_pre[ii+2])
                      +  8064.0f*(P_pre[ii+1] + P_pre[ii-1])
                      - 14350.0f*P_pre[ii])/(5040.0f*powf(dx,2.0f));                       

        float aux1 = powf(dt,2.0f)*powf(vp[ii],2.0f);                  
        float aux2 = 2.0f*P_pre[ii] - P_pas[ii];                       

        P_fut[ii] = aux1*d2_P_dx2 + aux2;  
    }
}

# endif