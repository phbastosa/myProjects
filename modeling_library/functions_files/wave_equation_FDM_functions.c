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

void acoustic_2D_2E2T(int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int z_src,
                      int x_src, int nxx,int nzz, float dx, float dz, float dt)
{
    int ii,jj,index;
    float d2_P_dx2, d2_P_dz2;

    for(index = 0; index < nxx*nzz; ++index) /* Spatial loop */ 
    {
        ii = floor(index / nxx);  /* Line indicator */
        jj = index % nxx;         /* Column indicator */  

        if(time_pointer < nsrc) P_pre[z_src*nxx + x_src] = source[time_pointer]; /* Applying source*/

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

void acoustic_2D_4E2T(int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int z_src,
                      int x_src, int nxx,int nzz, float dx, float dz, float dt)
{
    int ii,jj,index;
    float d2_P_dx2, d2_P_dz2;

    for(index = 0; index < nxx*nzz; ++index) /* Spatial loop */ 
    {
        ii = floor(index / nxx);  /* Line indicator */
        jj = index % nxx;         /* Column indicator */  

        if(time_pointer < nsrc) P_pre[z_src*nxx + x_src] = source[time_pointer]; /* Applying source*/

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

void acoustic_2D_8E2T(int time_pointer, float *vp, float *P_pre,
                      float *P_pas, float *P_fut, float *source, int nsrc, int z_src,
                      int x_src, int nxx,int nzz, float dx, float dz, float dt)
{
    int ii,jj,index;
    float d2_P_dx2, d2_P_dz2;

    for(index = 0; index < nxx*nzz; ++index) /* Spatial loop */ 
    {
        ii = floor(index / nxx);  /* Line indicator */
        jj = index % nxx;         /* Column indicator */  

        if(time_pointer < nsrc) P_pre[z_src*nxx + x_src] = source[time_pointer]; /* Applying source*/

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

void elastic_isotropic_2D_wave_8E2T_tension_stencil(float *U, float *V, float *Txx, float *Tzz, float *Txz,
    float *rho, float *M, float *L, int nxx, int nzz, float dt, float ds) {

    int index,ii,jj;
    float rho_int,L_int,M_int;
    float d_Txx_dx, d_Txz_dz, d_Txz_dx, d_Tzz_dz;
    float d_U_dx, d_V_dz, d_U_dz, d_V_dx;

    for(index = 0; index < nxx*nzz; index++) {              

        ii = (int) index / nxx;      /* Line indicator */
        jj = (int) index % nxx;      /* Column indicator */ 

        if((ii > 3) && (ii < nzz-3) && (jj > 3) && (jj < nxx-3)) { // vacuum (estencil na tensão)
         
            d_Txx_dx = (75*(Txx[(jj-3) + ii*nxx] - Txx[(jj+4) + ii*nxx]) +
                      1029*(Txx[(jj+3) + ii*nxx] - Txx[(jj-2) + ii*nxx]) +
                      8575*(Txx[(jj-1) + ii*nxx] - Txx[(jj+2) + ii*nxx]) +
                    128625*(Txx[(jj+1) + ii*nxx] - Txx[jj + ii*nxx]))/ds;

            d_Txz_dz = (75*(Txz[jj + (ii-4)*nxx] - Txz[jj + (ii+3)*nxx]) +
                      1029*(Txz[jj + (ii+2)*nxx] - Txz[jj + (ii-3)*nxx]) + 
                      8575*(Txz[jj + (ii-2)*nxx] - Txz[jj + (ii+1)*nxx]) +
                    128625*(Txz[jj + ii*nxx]     - Txz[jj + (ii-1)*nxx]))/ds;

            rho_int = 0.5*(rho[jj + ii*nxx] + rho[(jj+1) + ii*nxx]);

            U[jj + ii*nxx] += dt/(rho[jj + ii*nxx]*107520)*(d_Txx_dx + d_Txz_dz);  
        }
    }

    for(index = 0; index < nxx*nzz; index++) {              

        ii = (int) index / nxx;      /* Line indicator */
        jj = (int) index % nxx;      /* Column indicator */ 
      
        if((ii >= 3) && (ii < nzz-4) && (jj >= 3) && (jj < nxx-4)) { // vacuum

        
            d_Txz_dx = (75*(Txz[(jj-4) + ii*nxx] - Txz[(jj+3) + ii*nxx]) +
                      1029*(Txz[(jj+2) + ii*nxx] - Txz[(jj-3) + ii*nxx]) +
                      8575*(Txz[(jj-2) + ii*nxx] - Txz[(jj+1) + ii*nxx]) +
                    128625*(Txz[jj + ii*nxx]     - Txz[(jj-1) + ii*nxx]))/ds;

            d_Tzz_dz = (75*(Tzz[jj + (ii-3)*nxx] - Tzz[jj + (ii+4)*nxx]) + 
                      1029*(Tzz[jj + (ii+3)*nxx] - Tzz[jj + (ii-2)*nxx]) +
                      8575*(Tzz[jj + (ii-1)*nxx] - Tzz[jj + (ii+2)*nxx]) +
                    128625*(Tzz[jj + (ii+1)*nxx] - Tzz[jj + ii*nxx]))/ds;

            rho_int = 0.5*(rho[jj + ii*nxx] + rho[jj + (ii+1)*nxx]);

            V[jj + ii*nxx] += dt/(rho_int*107520)*(d_Txz_dx + d_Tzz_dz); 
        }
    }
    
    for(index = 0; index < nxx*nzz; index++) {              
        
        ii = (int) index / nxx;      /* Line indicator */
        jj = (int) index % nxx;      /* Column indicator */

        if((ii > 3) && (ii < nzz-3) && (jj >= 3) && (jj < nxx-4)) { //vacuum

            d_U_dx = (75*(U[(jj-4) + ii*nxx] - U[(jj+3) + ii*nxx]) + 
                    1029*(U[(jj+2) + ii*nxx] - U[(jj-3) + ii*nxx]) +
                    8575*(U[(jj-2) + ii*nxx] - U[(jj+1) + ii*nxx]) + 
                  128625*(U[jj + ii*nxx]     - U[(jj-1) + ii*nxx]))/ds;

            d_V_dz = (75*(V[jj + (ii-4)*nxx] - V[jj + (ii+3)*nxx]) +   
                    1029*(V[jj + (ii+2)*nxx] - V[jj + (ii-3)*nxx]) +
                    8575*(V[jj + (ii-2)*nxx] - V[jj + (ii+1)*nxx]) +
                  128625*(V[jj + ii*nxx]     - V[jj + (ii-1)*nxx]))/ds;     

            L_int = powf(0.25*(1/L[jj + (ii+1)*nxx] + 1/L[jj + ii*nxx] + 1/L[(jj+1) + ii*nxx] + 1/L[(jj+1) + (ii+1)*nxx]),(-1));
            M_int = powf(0.25*(1/M[jj + (ii+1)*nxx] + 1/M[jj + ii*nxx] + 1/M[(jj+1) + ii*nxx] + 1/M[(jj+1) + (ii+1)*nxx]),(-1));

            Txx[jj + ii*nxx] += (L_int + 2*M_int)*dt/107520*d_U_dx + 
                                 L_int*dt/107520*d_V_dz;   
        
            Tzz[jj+ ii*nxx] += (L_int + 2*M_int)*dt/107520*d_V_dz +  
                                L_int*dt/107520*d_U_dx;
        }
    }

    for(index = 0; index < nxx*nzz; index++) {              
        
        ii = (int) index / nxx;      /* Line indicator */
        jj = (int) index % nxx;      /* Column indicator */

        if((ii >= 3) && (ii < nzz-4) && (jj > 3) && (jj < nxx-3)) {  // vacuum

            d_U_dz = (75*(U[jj + (ii-3)*nxx] - U[jj + (ii+4)*nxx]) +
                    1029*(U[jj + (ii+3)*nxx] - U[jj + (ii-2)*nxx]) +
                    8575*(U[jj + (ii-1)*nxx] - U[jj + (ii+2)*nxx]) +
                  128625*(U[jj + (ii+1)*nxx] - U[jj + ii*nxx]))/ds;

            d_V_dx = (75*(V[(jj-3) + ii*nxx] - V[(jj+4) + ii*nxx]) +
                    1029*(V[(jj+3) + ii*nxx] - V[(jj-2) + ii*nxx]) +
                    8575*(V[(jj-1) + ii*nxx] - V[(jj+2) + ii*nxx]) +
                  128625*(V[(jj+1) + ii*nxx] - V[jj + ii*nxx]))/ds;

            Txz[jj + ii*nxx] += M[jj + ii*nxx]*dt/107520*(d_U_dz + d_V_dx);            
        }      
    }
}

void elastic_isotropic_2D_wave_8E2T_velocity_stencil(float *U, float *V, float *Txx, float *Tzz, float *Txz,
                                float *rho, float *M, float *L, int nxx, int nzz, float dt,
                                float ds) {

    int index,ii,jj;
    float rho_int,L_int,M_int;
    float d_Txx_dx, d_Txz_dz, d_Txz_dx, d_Tzz_dz;
    float d_U_dx, d_V_dz, d_U_dz, d_V_dx; 

    #pragma acc parallel loop independent present(U[0:nxx*nzz],rho[0:nxx*nzz],Txx[0:nxx*nzz],Txz[0:nxx*nzz]) copyin(Tzz[0:nzz*nxx]) copyout(U[0:nxx*nzz])
    for(index = 0; index < nxx*nzz; index++) {              

        ii = (int) index / nxx;      // indicador de linhas   (direção y)
        jj = (int) index % nxx;      // indicador de colunas  (direção x) 

        if((ii > 3) && (ii < nzz-3) && (jj >= 3) && (jj < nxx-4)) {
        
            d_Txx_dx = (75*(Txx[(jj-3) + ii*nxx] - Txx[(jj+4) + ii*nxx]) +
                      1029*(Txx[(jj+3) + ii*nxx] - Txx[(jj-2) + ii*nxx]) +
                      8575*(Txx[(jj-1) + ii*nxx] - Txx[(jj+2) + ii*nxx]) +
                    128625*(Txx[(jj+1) + ii*nxx] - Txx[jj + ii*nxx]))/ds;

            d_Txz_dz = (75*(Txz[jj + (ii-4)*nxx] - Txz[jj + (ii+3)*nxx]) +
                      1029*(Txz[jj + (ii+2)*nxx] - Txz[jj + (ii-3)*nxx]) + 
                      8575*(Txz[jj + (ii-2)*nxx] - Txz[jj + (ii+1)*nxx]) +
                    128625*(Txz[jj + ii*nxx]     - Txz[jj + (ii-1)*nxx]))/ds;

            U[jj + ii*nxx] += dt/(rho[jj + ii*nxx]*107520)*(d_Txx_dx + d_Txz_dz);  
        }
    }

    #pragma acc parallel loop independent present(V[0:nxx*nzz],Txz[0:nxx*nzz],Tzz[0:nxx*nzz],rho[0:nxx*nzz]) copyout(V[0:nxx*nzz]) 
    for(index = 0; index < nxx*nzz; index++) {              

        ii = (int) index / nxx;      // indicador de linhas  (direção y)
        jj = (int) index % nxx;      // indicador de colunas (direção x) 

        if((ii >= 3) && (ii < nzz-4) && (jj > 3) && (jj < nxx-3)) {
        
            d_Txz_dx = (75*(Txz[(jj-4) + ii*nxx] - Txz[(jj+3) + ii*nxx]) +
                      1029*(Txz[(jj+2) + ii*nxx] - Txz[(jj-3) + ii*nxx]) +
                      8575*(Txz[(jj-2) + ii*nxx] - Txz[(jj+1) + ii*nxx]) +
                    128625*(Txz[jj + ii*nxx]     - Txz[(jj-1) + ii*nxx]))/ds;

            d_Tzz_dz = (75*(Tzz[jj + (ii-3)*nxx] - Tzz[jj + (ii+4)*nxx]) + 
                      1029*(Tzz[jj + (ii+3)*nxx] - Tzz[jj + (ii-2)*nxx]) +
                      8575*(Tzz[jj + (ii-1)*nxx] - Tzz[jj + (ii+2)*nxx]) +
                    128625*(Tzz[jj + (ii+1)*nxx] - Tzz[jj + ii*nxx]))/ds;

            rho_int = 0.25*(rho[jj + ii*nxx] + rho[jj + (ii+1)*nxx] + rho[(jj+1) + (ii+1)*nxx] + rho[(jj+1) + ii*nxx]);

            V[jj + ii*nxx] += dt/(rho_int*107520)*(d_Txz_dx + d_Tzz_dz); 
        }
    }
    
    #pragma acc parallel loop independent present(U[0:nxx*nzz],V[0:nxx*nzz],Txx[0:nxx*nzz],Tzz[0:nxx*nzz],L[0:nxx*nzz],M[0:nxx*nzz]) copyout(Txx[0:nxx*nzz],Tzz[0:nxx*nzz])
    for(index = 0; index < nxx*nzz; index++) {              
        
        ii = (int) index / nxx;      // indicador de linhas  (direção y)
        jj = (int) index % nxx;      // indicador de colunas (direção x)

        if((ii > 3) && (ii < nzz-3) && (jj > 3) && (jj < nxx-3)) {

            L_int = 0.5*(L[jj + (ii+1)*nxx] + L[jj + ii*nxx]);
            M_int = 0.5*(M[jj + (ii+1)*nxx] + M[jj + ii*nxx]);

            d_U_dx = (75*(U[(jj-4) + ii*nxx] - U[(jj+3) + ii*nxx]) + 
                    1029*(U[(jj+2) + ii*nxx] - U[(jj-3) + ii*nxx]) +
                    8575*(U[(jj-2) + ii*nxx] - U[(jj+1) + ii*nxx]) + 
                  128625*(U[jj + ii*nxx]     - U[(jj-1) + ii*nxx]))/ds;

            d_V_dz = (75*(V[jj + (ii-4)*nxx] - V[jj + (ii+3)*nxx]) +   
                    1029*(V[jj + (ii+2)*nxx] - V[jj + (ii-3)*nxx]) +
                    8575*(V[jj + (ii-2)*nxx] - V[jj + (ii+1)*nxx]) +
                  128625*(V[jj + ii*nxx]     - V[jj + (ii-1)*nxx]))/ds;     

            Txx[jj + ii*nxx] += (L_int + 2*M_int)*dt/107520*d_U_dx + 
                                 L_int*dt/107520*d_V_dz;   
        
            Tzz[jj+ ii*nxx] += (L_int + 2*M_int)*dt/107520*d_V_dz +  
                                L_int*dt/107520*d_U_dx;
        }
    }

    #pragma acc parallel loop independent present(Txz[0:nzz*nxx],M[0:nzz*nxx],U[0:nzz*nxx],V[0:nzz*nxx]) copyout(Txz[0:nzz*nxx])
    for(index = 0; index < nxx*nzz; index++) {              
        
        ii = (int) index / nxx;      // indicador de linhas  (direção y)
        jj = (int) index % nxx;      // indicador de colunas (direção x)

        if((ii >= 3) && (ii < nzz-4) && (jj >= 3) && (jj < nxx-4)) {

            M_int = 0.5*(M[(jj+1) + ii*nxx] + M[jj + ii*nxx]); 

            d_U_dz = (75*(U[jj + (ii-3)*nxx] - U[jj + (ii+4)*nxx]) +
                    1029*(U[jj + (ii+3)*nxx] - U[jj + (ii-2)*nxx]) +
                    8575*(U[jj + (ii-1)*nxx] - U[jj + (ii+2)*nxx]) +
                  128625*(U[jj + (ii+1)*nxx] - U[jj + ii*nxx]))/ds;

            d_V_dx = (75*(V[(jj-3) + ii*nxx] - V[(jj+4) + ii*nxx]) +
                    1029*(V[(jj+3) + ii*nxx] - V[(jj-2) + ii*nxx]) +
                    8575*(V[(jj-1) + ii*nxx] - V[(jj+2) + ii*nxx]) +
                  128625*(V[(jj+1) + ii*nxx] - V[jj + ii*nxx]))/ds;

            Txz[jj + ii*nxx] += M_int*dt/107520*(d_U_dz + d_V_dx);            
        }
    }
}
