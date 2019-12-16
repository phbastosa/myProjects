#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void update_wavefield(float * P_pas, float * P_pre, float * P_fut, int length)
{
    int index;
    
    for(index = 0; index < length; index++) {
        P_pas[index] = P_pre[index]; 
        P_pre[index] = P_fut[index];    
    }
}

void acoustic_1D_2E2T(float * P_pas, float * P_pre,float * P_fut, float * Vp, int Nx, float dx, float dt)
{
    int ii;
    float aux1, aux2, d2_P_dx2;

    for(ii = 1; ii < Nx-1; ii++) 
    {
        d2_P_dx2 = (P_pre[ii+1] - 2.0*P_pre[ii] + P_pre[ii-1])/pow(dx,2.0); 
        aux1 = pow(dt,2.0)*pow(Vp[ii],2.0);                  
        aux2 = 2.0*P_pre[ii] - P_pas[ii];

        P_fut[ii] = aux1*d2_P_dx2 + aux2;  
    }
}

void acoustic_1D_4E2T(float * P_pas, float * P_pre,float * P_fut, float * Vp, int Nx, float dx, float dt)
{
    int ii;
    float aux1, aux2, d2_P_dx2;

    for(ii = 2; ii < Nx-2; ii++) 
    {
        d2_P_dx2 = (-P_pre[ii-2] + 16.0f*P_pre[ii-1] - 30.0f*P_pre[ii] + 16.0f*P_pre[ii + 1] - P_pre[ii + 2])/(12.0f*pow(dx,2.0)); 
        aux1 = pow(dt,2.0)*pow(Vp[ii],2.0);                  
        aux2 = 2.0*P_pre[ii] - P_pas[ii];                       

        P_fut[ii] = aux1*d2_P_dx2 + aux2;  
    }
}

void acoustic_1D_8E2T(float * P_pas, float * P_pre,float * P_fut, float * Vp, int Nx, float dx, float dt)
{
    int ii;
    float aux1, aux2, d2_P_dx2;

    for(ii = 4; ii < Nx-4; ii++) 
    {
        d2_P_dx2 = (- 9.0f*(P_pre[ii-4] + P_pre[ii+4])
                +   128.0f*(P_pre[ii-3] + P_pre[ii+3])
                -  1008.0f*(P_pre[ii-2] + P_pre[ii+2])
                +  8064.0f*(P_pre[ii+1] + P_pre[ii-1])
                - 14350.0f*P_pre[ii])/(5040.0f*powf(dx,2));                       

        aux1 = pow(dt,2.0)*pow(Vp[ii],2.0);                  
        aux2 = 2.0*P_pre[ii] - P_pas[ii];                       

        P_fut[ii] = aux1*d2_P_dx2 + aux2;  
    }
}

void acoustic_2D_2E2T(int shot_pointer, int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int z_src,
                      int *x_src, int n_shots, int nxx,int nzz, float dx, float dz, float dt)
{
    int ii,jj,index;
    float d2_P_dx2, d2_P_dz2;

    for(index = 0; index < nxx*nzz; ++index) /* Spatial loop */ 
    {
        ii = floor(index / nxx);  /* Line indicator */
        jj = index % nxx;         /* Column indicator */  

        if(time_pointer < nsrc) P_pre[z_src*nxx + x_src[shot_pointer]] = source[time_pointer]; /* Applying source*/

        if((ii > 0) && (ii < nzz-1) && (jj > 0) && (jj < nxx-1)) 
        {
            /* Second derivative of the pressure with respect to x, 
            using second-order finite difference discretization */
            d2_P_dx2 = (P_pre[ii*nxx + (jj-1)] - 2.0f*P_pre[ii*nxx + jj] + P_pre[ii*nxx + (jj+1)])/powf(dx,2);

            /* Second derivative of the pressure with respect to z, 
            using second-order finite difference discretization */    
            d2_P_dz2 = (P_pre[(ii-1)*nxx + jj] - 2.0f*P_pre[ii*nxx + jj] + P_pre[(ii+1)*nxx + jj])/powf(dz,2);

            /* Calculating the future wave field through the previous fields */
            P_fut[ii*nxx + jj] = (powf(dt,2)*powf(vp[ii*nxx + jj],2)
                               * (d2_P_dx2 + d2_P_dz2)) + 2.0f*P_pre[ii*nxx + jj] 
                               -  P_pas[ii*nxx + jj];
        }
    }
}

void acoustic_2D_4E2T(int shot_pointer, int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int z_src,
                      int *x_src, int n_shots, int nxx,int nzz, float dx, float dz, float dt)
{
    int ii,jj,index;
    float d2_P_dx2, d2_P_dz2;

    for(index = 0; index < nxx*nzz; ++index) /* Spatial loop */ 
    {
        ii = floor(index / nxx);  /* Line indicator */
        jj = index % nxx;         /* Column indicator */  

        if(time_pointer < nsrc) P_pre[z_src*nxx + x_src[shot_pointer]] = source[time_pointer]; /* Applying source*/

        if((ii > 1) && (ii < nzz-2) && (jj > 1) && (jj < nxx-2)) 
        {
            /* Second derivative of the pressure with respect to x, 
            using fourth-order finite difference discretization */
            d2_P_dx2 = (- P_pre[ii*nxx + (jj-2)] + 16.0f*P_pre[ii*nxx + (jj-1)] 
                -   30.0f*P_pre[ii*nxx + jj] + 16.0f*P_pre[ii*nxx + (jj+1)] 
                -         P_pre[ii*nxx + (jj+2)])/(12.0f*powf(dx,2));

            /* Second derivative of the pressure with respect to z, 
            using fourth-order finite difference discretization */    
            d2_P_dz2 = (- P_pre[(ii-2)*nxx + jj] + 16.0f*P_pre[(ii-1)*nxx + jj] 
                -   30.0f*P_pre[ii*nxx + jj] + 16.0f*P_pre[(ii+1)*nxx + jj] 
                -         P_pre[(ii+2)*nxx + jj])/(12.0f*powf(dz,2));

            /* Calculating the future wave field through the previous fields */
            P_fut[ii*nxx + jj] = (powf(dt,2)*powf(vp[ii*nxx + jj],2)
                               * (d2_P_dx2 + d2_P_dz2)) + 2.0f*P_pre[ii*nxx + jj] 
                               -  P_pas[ii*nxx + jj];
        }
    }
}

void acoustic_2D_8E2T(int shot_pointer, int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int z_src,
                      int *x_src, int n_shots, int nxx,int nzz, float dx, float dz, float dt)
{
    int ii,jj,index;
    float d2_P_dx2, d2_P_dz2;

    for(index = 0; index < nxx*nzz; ++index) /* Spatial loop */ 
    {
        ii = floor(index / nxx);  /* Line indicator */
        jj = index % nxx;         /* Column indicator */  

        if(time_pointer < nsrc) P_pre[z_src*nxx + x_src[shot_pointer]] = source[time_pointer]; /* Applying source*/

        if((ii > 3) && (ii < nzz-4) && (jj > 3) && (jj < nxx-4)) 
        {
            /* Second derivative of the pressure with respect to x, 
            using eighth-order finite difference discretization */
            d2_P_dx2 = (- 9.0f*(P_pre[ii*nxx + (jj-4)] + P_pre[ii*nxx + (jj+4)])
                +   128.0f*(P_pre[ii*nxx + (jj-3)] + P_pre[ii*nxx + (jj+3)])
                -  1008.0f*(P_pre[ii*nxx + (jj-2)] + P_pre[ii*nxx + (jj+2)])
                +  8064.0f*(P_pre[ii*nxx + (jj+1)] + P_pre[ii*nxx + (jj-1)])
                - 14350.0f*(P_pre[ii*nxx + jj]))/(5040.0f*powf(dx,2));

            /* Second derivative of the pressure with respect to z, 
            using eighth-order finite difference discretization */    
            d2_P_dz2 = (- 9.0f*(P_pre[(ii-4)*nxx + jj] + P_pre[(ii+4)*nxx + jj])
                +   128.0f*(P_pre[(ii-3)*nxx + jj] + P_pre[(ii+3)*nxx + jj])
                -  1008.0f*(P_pre[(ii-2)*nxx + jj] + P_pre[(ii+2)*nxx + jj])
                +  8064.0f*(P_pre[(ii-1)*nxx + jj] + P_pre[(ii+1)*nxx + jj])
                - 14350.0f*(P_pre[ii*nxx + jj]))/(5040.0f*powf(dz,2));

            /* Calculating the future wave field through the previous fields */
            P_fut[ii*nxx + jj] = (powf(dt,2)*powf(vp[ii*nxx + jj],2)
                            * (d2_P_dx2 + d2_P_dz2)) + 2.0f*P_pre[ii*nxx + jj] 
                            -  P_pas[ii*nxx + jj];
        }
    }
}

void acoustic_3D_2E2T(int shot_pointer, int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int *z_src,
                      int *y_src, int *x_src, int n_shots, int nxx, int nyy, int nzz, 
                      float dx, float dy, float dz, float dt)
{
    int index, ii, jj, kk;
    float d2_P_dx2, d2_P_dy2, d2_P_dz2;

    for(index = 0; index < nxx*nyy*nzz; index++) {    

        kk = floor(index/(nxx*nyy));           /* Matrix indicator (z direction) */
        jj = index % nyy;                      /* Column indicator (x direction) */
        ii = floor((index % (nxx*nyy)) / nxx); /* Line indicator   (y direction) */  

        if(time_pointer < nsrc) P_pre[z_src[shot_pointer]*nyy*nxx + y_src[shot_pointer]*nxx + x_src[shot_pointer]] = source[time_pointer]; /* Applying source*/

        if((ii > 0) && (ii < nyy-1) && (jj > 0) && (jj < nxx-1) && (kk > 0) && (kk < nzz-1)) 
        {
            /* Second derivative of the pressure with respect to x, 
            using second-order finite difference discretization */
            d2_P_dx2 = P_pre[kk*nyy*nxx + ii*nxx + (jj-1)] - 2.0f*P_pre[kk*nyy*nxx + ii*nxx + jj] + P_pre[kk*nyy*nxx + ii*nxx + (jj+1)]/powf(dx,2.0);

            /* Second derivative of the pressure with respect to y, 
            using second-order finite difference discretization */
            d2_P_dy2 = P_pre[kk*nyy*nxx + (ii-1)*nxx + jj] - 2.0f*P_pre[kk*nyy*nxx + ii*nxx + jj] + P_pre[kk*nyy*nxx + (ii+1)*nxx + jj]/powf(dy,2.0);

            /* Second derivative of the pressure with respect to z, 
            using second-order finite difference discretization */    
            d2_P_dz2 = P_pre[(kk-1)*nyy*nxx + ii*nxx + jj] - 2.0f*P_pre[kk*nyy*nxx + ii*nxx + jj] + P_pre[(kk+1)*nyy*nxx + ii*nxx + jj]/powf(dz,2.0);

            /* Calculating the future wave field through the previous fields */
            P_fut[ii*nxx + jj] = (powf(dt,2)*powf(vp[kk*nxx*nyy + ii*nxx + jj],2)
                            * (d2_P_dx2 + d2_P_dy2 + d2_P_dz2)) + 2.0f*P_pre[kk*nxx*nyy + ii*nxx + jj] 
                            -  P_pas[kk*nxx*nyy + ii*nxx + jj];
        }
    }
}

void acoustic_3D_4E2T(int shot_pointer, int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int *z_src,
                      int *y_src, int *x_src, int n_shots, int nxx, int nyy, int nzz, 
                      float dx, float dy, float dz, float dt)
{
    int index, ii, jj, kk;
    float d2_P_dx2, d2_P_dy2, d2_P_dz2;

    for(index = 0; index < nxx*nyy*nzz; index++) {    

        kk = floor(index/(nxx*nyy));           /* Matrix indicator (z direction) */
        jj = index % nyy;                      /* Column indicator (x direction) */
        ii = floor((index % (nxx*nyy)) / nxx); /* Line indicator   (y direction) */  

        if(time_pointer < nsrc) P_pre[z_src[shot_pointer]*nyy*nxx + y_src[shot_pointer]*nxx + x_src[shot_pointer]] = source[time_pointer]; /* Applying source*/

        if((ii > 1) && (ii < nyy-2) && (jj > 1) && (jj < nxx-2) && (kk > 1) && (kk < nzz-2)) 
        {
            /* Second derivative of the pressure with respect to x, 
            using fourth-order finite difference discretization */
            d2_P_dx2 = (- P_pre[kk*nxx*nyy + ii*nxx + (jj-2)] + 16.0f*P_pre[kk*nxx*nyy + ii*nxx + (jj-1)] 
                -   30.0f*P_pre[kk*nxx*nyy + ii*nxx + jj] + 16.0f*P_pre[kk*nxx*nyy + ii*nxx + (jj+1)] 
                -         P_pre[kk*nxx*nyy + ii*nxx + (jj+2)])/(12.0f*powf(dx,2));
            
            /* Second derivative of the pressure with respect to y, 
            using fourth-order finite difference discretization */
            d2_P_dy2 = (- P_pre[kk*nxx*nyy + (ii-2)*nxx + jj] + 16.0f*P_pre[kk*nxx*nyy + (ii-1)*nxx + jj] 
                -   30.0f*P_pre[kk*nxx*nyy + ii*nxx + jj] + 16.0f*P_pre[kk*nxx*nyy + (ii+1)*nxx + jj] 
                -         P_pre[kk*nxx*nyy + (ii+2)*nxx + jj])/(12.0f*powf(dy,2));
            
            /* Second derivative of the pressure with respect to z, 
            using fourth-order finite difference discretization */    
            d2_P_dz2 = (- P_pre[(kk-2)*nxx*nyy + ii*nxx + jj] + 16.0f*P_pre[(kk-1)*nxx*nyy + ii*nxx + jj] 
                -   30.0f*P_pre[kk*nxx*nyy + ii*nxx + jj] + 16.0f*P_pre[(kk+1)*nxx*nyy + ii*nxx + jj] 
                -         P_pre[(kk+2)*nxx*nyy + ii*nxx + jj])/(12.0f*powf(dz,2));

            /* Calculating the future wave field through the previous fields */
            P_fut[ii*nxx + jj] = (powf(dt,2)*powf(vp[kk*nxx*nyy + ii*nxx + jj],2)
                            * (d2_P_dx2 + d2_P_dy2 + d2_P_dz2)) + 2.0f*P_pre[kk*nxx*nyy + ii*nxx + jj] 
                            -  P_pas[kk*nxx*nyy + ii*nxx + jj];
        }
    }
}

void acoustic_3D_8E2T(int shot_pointer, int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int *z_src,
                      int *y_src, int *x_src, int n_shots, int nxx, int nyy, int nzz, 
                      float dx, float dy, float dz, float dt)
{
    int index, ii, jj, kk;
    float d2_P_dx2, d2_P_dy2, d2_P_dz2;

    for(index = 0; index < nxx*nyy*nzz; index++) {    

        kk = floor(index/(nxx*nyy));           /* Matrix indicator (z direction) */
        jj = index % nyy;                      /* Column indicator (x direction) */
        ii = floor((index % (nxx*nyy)) / nxx); /* Line indicator   (y direction) */  

        if(time_pointer < nsrc) P_pre[z_src[shot_pointer]*nyy*nxx + y_src[shot_pointer]*nxx + x_src[shot_pointer]] = source[time_pointer]; /* Applying source*/

        if((ii > 3) && (ii < nyy-4) && (jj > 3) && (jj < nxx-4) && (kk > 3) && (kk < nzz-4)) 
        {
            /* Second derivative of the pressure with respect to x, 
            using fourth-order finite difference discretization */
            d2_P_dx2 = (- 9.0f*(P_pre[kk*nxx*nyy + ii*nxx + (jj-4)] + P_pre[kk*nxx*nyy + ii*nxx + (jj+4)])
                    +   128.0f*(P_pre[kk*nxx*nyy + ii*nxx + (jj-3)] + P_pre[kk*nxx*nyy + ii*nxx + (jj+3)])
                    -  1008.0f*(P_pre[kk*nxx*nyy + ii*nxx + (jj-2)] + P_pre[kk*nxx*nyy + ii*nxx + (jj+2)])
                    +  8064.0f*(P_pre[kk*nxx*nyy + ii*nxx + (jj+1)] + P_pre[kk*nxx*nyy + ii*nxx + (jj-1)])
                    - 14350.0f*(P_pre[kk*nxx*nyy + ii*nxx + jj]))/(5040.0f*powf(dx,2));

            /* Second derivative of the pressure with respect to y, 
            using fourth-order finite difference discretization */
            d2_P_dy2 = (- 9.0f*(P_pre[kk*nxx*nyy + (ii-4)*nxx + jj] + P_pre[kk*nxx*nyy + (ii+4)*nxx + jj])
                    +   128.0f*(P_pre[kk*nxx*nyy + (ii-3)*nxx + jj] + P_pre[kk*nxx*nyy + (ii+3)*nxx + jj])
                    -  1008.0f*(P_pre[kk*nxx*nyy + (ii-2)*nxx + jj] + P_pre[kk*nxx*nyy + (ii+2)*nxx + jj])
                    +  8064.0f*(P_pre[kk*nxx*nyy + (ii+1)*nxx + jj] + P_pre[kk*nxx*nyy + (ii-1)*nxx + jj])
                    - 14350.0f*(P_pre[kk*nxx*nyy + ii*nxx + jj]))/(5040.0f*powf(dy,2));

            /* Second derivative of the pressure with respect to z, 
            using fourth-order finite difference discretization */    
            d2_P_dz2 = (- 9.0f*(P_pre[(kk-4)*nxx*nyy + ii*nxx + jj] + P_pre[(kk+4)*nxx*nyy + ii*nxx + jj])
                    +   128.0f*(P_pre[(kk-3)*nxx*nyy + ii*nxx + jj] + P_pre[(kk+3)*nxx*nyy + ii*nxx + jj])
                    -  1008.0f*(P_pre[(kk-2)*nxx*nyy + ii*nxx + jj] + P_pre[(kk+2)*nxx*nyy + ii*nxx + jj])
                    +  8064.0f*(P_pre[(kk+1)*nxx*nyy + ii*nxx + jj] + P_pre[(kk-1)*nxx*nyy + ii*nxx + jj])
                    - 14350.0f*(P_pre[kk*nxx*nyy + ii*nxx + jj]))/(5040.0f*powf(dy,2));

            /* Calculating the future wave field through the previous fields */
            P_fut[ii*nxx + jj] = (powf(dt,2)*powf(vp[kk*nxx*nyy + ii*nxx + jj],2)
                            * (d2_P_dx2 + d2_P_dy2 + d2_P_dz2)) + 2.0f*P_pre[kk*nxx*nyy + ii*nxx + jj] 
                            -  P_pas[kk*nxx*nyy + ii*nxx + jj];
        }
    }
}
