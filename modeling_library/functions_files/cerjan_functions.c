# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "../header_files/cerjan_functions.h"

void Cerjan_1D_acoustic_attenuation(float * P_pre, float * P_fut, float * factor, int abc_layer, int Nx)
{
    int index;

    for (index = 0; index < abc_layer; ++index) 
    {
        P_pre[index] = P_pre[index] * factor[index];
        P_fut[index] = P_fut[index] * factor[index];
        P_pre[Nx-index-1] = P_pre[Nx-index-1] * factor[index]; 
        P_fut[Nx-index-1] = P_fut[Nx-index-1] * factor[index]; 
    }
}

void Cerjan_2D_corners(float *up_left, float *up_right, float *down_left, 
                       float *down_right, float *factor, int n_abc) 
{
    int ii,jj,kk;

    /* Upper left absorbing corner*/
    for(ii = n_abc-1; ii >= 0; ii++) 
    {
        up_left[ii*n_abc + ii] = factor[ii];    
        for(kk = 1; kk <= n_abc; kk++) 
        {
            up_left[(ii+kk)*n_abc + ii] = factor[ii];
            up_left[ii*n_abc + ii + kk] = factor[ii];
        }    
    }

    /* Upper right absorbing corner */
    for(ii = 0; ii < n_abc; ii++) 
    {
        for(jj = 0; jj < n_abc; jj++) {
            up_right[ii*n_abc + jj] = up_left[ii*n_abc + n_abc-jj-1];
        }
    }

    /* Botton left absorbing corner */
    for(ii = 0; ii < n_abc; ii++) 
    {
        for(jj = 0; jj < n_abc; jj++) 
        {
            down_left[ii*n_abc + jj] = up_left[(n_abc-ii-1)*n_abc + jj];
        }
    }

    /* Botton right absorbing corner */ 
    for(ii = 0; ii < n_abc; ii++) 
    {
        for(jj = 0; jj < n_abc; jj++) 
        {
            down_right[ii*n_abc + jj] = down_left[ii*n_abc + n_abc-jj-1];
        }
    }
}

void Cerjan_2D_acoustic_attenuation(float *P_pre, float *P_fut, float *factor, 
                                    int nxx, int nzz, int abc_layer, float *down_left, 
                                    float *down_right, float *up_left, float *up_right)
{
    int ii,jj,index;
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
            P_pre[ii*nxx + jj] = P_pre[ii*nxx + jj]*up_left[ii*abc_layer + jj];
            P_fut[ii*nxx + jj] = P_fut[ii*nxx + jj]*up_left[ii*abc_layer + jj];
        }

        /* Upper right absorbing corner */
        if((ii >= 0) && (ii < abc_layer) && (jj >= nxx-abc_layer) && (jj < nxx)) 
        {
            P_pre[ii*nxx + jj] = P_pre[ii*nxx + jj]*up_right[ii*abc_layer + (jj-nxx+abc_layer)];
            P_fut[ii*nxx + jj] = P_fut[ii*nxx + jj]*up_right[ii*abc_layer + (jj-nxx+abc_layer)];
        }

        /* Bottom left absorbing corner */
        if((ii >= nzz-abc_layer) && (ii < nzz) && (jj >= 0) && (jj < abc_layer)) 
        {
            P_pre[ii*nxx + jj] = P_pre[ii*nxx + jj]*down_left[(ii-nzz+abc_layer)*abc_layer + jj];
            P_fut[ii*nxx + jj] = P_fut[ii*nxx + jj]*down_left[(ii-nzz+abc_layer)*abc_layer + jj];
        }

        /* Botton right absorbing corner */ 
        if((ii >= nzz-abc_layer) && (ii < nzz) && (jj >= nxx-abc_layer) && (jj < nxx)) 
        {
            P_pre[ii*nxx + jj] = P_pre[ii*nxx + jj]*down_right[(ii-nzz+abc_layer)*abc_layer + (jj-nxx+abc_layer)];
            P_fut[ii*nxx + jj] = P_fut[ii*nxx + jj]*down_right[(ii-nzz+abc_layer)*abc_layer + (jj-nxx+abc_layer)];
        }
    }
}

void Cerjan_2D_elastic_attenuation(float *U, float *V, float *Txx, float *Tzz, 
                                   float *Txz, int nxx, int nzz, float *factor,
                                   float *up_left, float *up_right, float *down_left,
                                   float *down_right, int abc_layer) 
{
    int index,ii,jj,kk;

    for(index = 0; index < nxx*nzz; index++) { 
        
        ii = floor(index/nxx);  /* Line indicator */
        jj = index % nxx;       /* Column indicator */

        /* Left and right absorbing layers */
        if((ii >= abc_layer) && (ii < nzz-abc_layer) && (jj >= 0) && (jj < abc_layer)) 
        {
            U[ii*nxx+jj] = U[ii*nxx+jj]*factor[jj];
            U[ii*nxx+(nxx-jj)] = U[ii*nxx+(nxx-jj)]*factor[jj];

            V[ii*nxx+jj] = V[ii*nxx+jj]*factor[jj];
            V[ii*nxx+(nxx-jj)] = V[ii*nxx+(nxx-jj)]*factor[jj];

            Txx[ii*nxx+jj] = Txx[ii*nxx+jj]*factor[jj];
            Txx[ii*nxx+(nxx-jj)] = Txx[ii*nxx+(nxx-jj)]*factor[jj];

            Tzz[ii*nxx+jj] = Tzz[ii*nxx+jj]*factor[jj];
            Tzz[ii*nxx+(nxx-jj)] = Tzz[ii*nxx+(nxx-jj)]*factor[jj];

            Txz[ii*nxx+jj] = Txz[ii*nxx+jj]*factor[jj];
            Txz[ii*nxx+(nxx-jj)] = Txz[ii*nxx+(nxx-jj)]*factor[jj];
        }

        /* Up and bottom absorbing layers */
        if((jj >= abc_layer) && (jj < nxx-abc_layer) && (ii >= 0) && (ii < abc_layer)) {
            U[ii*nxx+jj] = U[ii*nxx+jj]*factor[ii];
            U[(nzz-ii)*nxx-jj] = U[(nzz-ii)*nxx-jj]*factor[ii];

            V[ii*nxx+jj] = V[ii*nxx+jj]*factor[ii];
            V[(nzz-ii)*nxx-jj] = V[(nzz-ii)*nxx-jj]*factor[ii];
            
            Txx[ii*nxx+jj] = Txx[ii*nxx+jj]*factor[ii];
            Txx[(nzz-ii)*nxx-jj] = Txx[(nzz-ii)*nxx-jj]*factor[ii];
            
            Tzz[ii*nxx+jj] = Tzz[ii*nxx+jj]*factor[ii];
            Tzz[(nzz-ii)*nxx-jj] = Tzz[(nzz-ii)*nxx-jj]*factor[ii];
            
            Txz[ii*nxx+jj] = Txz[ii*nxx+jj]*factor[ii];            
            Txz[(nzz-ii)*nxx-jj] = Txz[(nzz-ii)*nxx-jj]*factor[ii];
        }
           
        /* Upper left absorbing corner */ 
        if((ii >= 0) && (ii < abc_layer) && (jj >= 0) && (jj < abc_layer)) {
            U[ii*nxx + jj] = U[ii*nxx + jj]*up_left[ii*abc_layer + jj];
            V[ii*nxx + jj] = V[ii*nxx + jj]*up_left[ii*abc_layer + jj];
            Txx[ii*nxx + jj] = Txx[ii*nxx + jj]*up_left[ii*abc_layer + jj];
            Tzz[ii*nxx + jj] = Tzz[ii*nxx + jj]*up_left[ii*abc_layer + jj];
            Txz[ii*nxx + jj] = Txz[ii*nxx + jj]*up_left[ii*abc_layer + jj];
        }

        /* Upper right absorbing corner */
        if((ii >= 0) && (ii < abc_layer) && (jj >= nxx-abc_layer) && (jj < nxx)) {
            U[ii*nxx + jj] = U[ii*nxx + jj]*up_right[ii*abc_layer + (jj-nxx+abc_layer)];
            V[ii*nxx + jj] = V[ii*nxx + jj]*up_right[ii*abc_layer + (jj-nxx+abc_layer)];
            Txx[ii*nxx + jj] = Txx[ii*nxx + jj]*up_right[ii*abc_layer + (jj-nxx+abc_layer)];
            Tzz[ii*nxx + jj] = Tzz[ii*nxx + jj]*up_right[ii*abc_layer + (jj-nxx+abc_layer)];
            Txz[ii*nxx + jj] = Txz[ii*nxx + jj]*up_right[ii*abc_layer + (jj-nxx+abc_layer)];
        }

        /* Bottom left absorbing corner */
        if((ii >= nzz-abc_layer) && (ii < nzz) && (jj >= 0) && (jj < abc_layer)) {
            U[ii*nxx + jj] = U[ii*nxx + jj]*down_left[(ii-nzz+abc_layer)*abc_layer + jj];
            Txx[ii*nxx + jj] = Txx[ii*nxx + jj]*down_left[(ii-nzz+abc_layer)*abc_layer + jj];
            Tzz[ii*nxx + jj] = Tzz[ii*nxx + jj]*down_left[(ii-nzz+abc_layer)*abc_layer + jj];
            V[ii*nxx + jj] = V[ii*nxx + jj]*down_left[(ii-nzz+abc_layer)*abc_layer + jj];
            Txz[ii*nxx + jj] = Txz[ii*nxx + jj]*down_left[(ii-nzz+abc_layer)*abc_layer + jj];
        }

        /* Botton right absorbing corner */
        if((ii >= nzz-abc_layer) && (ii < nzz) && (jj >= nxx-abc_layer) && (jj < nxx)) {
            U[ii*nxx + jj] = U[ii*nxx + jj]*down_right[(ii-nzz+abc_layer)*abc_layer + (jj-nxx+abc_layer)];
            V[ii*nxx + jj] = V[ii*nxx + jj]*down_right[(ii-nzz+abc_layer)*abc_layer + (jj-nxx+abc_layer)];
            Txx[ii*nxx + jj] = Txx[ii*nxx + jj]*down_right[(ii-nzz+abc_layer)*abc_layer + (jj-nxx+abc_layer)];
            Tzz[ii*nxx + jj] = Tzz[ii*nxx + jj]*down_right[(ii-nzz+abc_layer)*abc_layer + (jj-nxx+abc_layer)];
            Txz[ii*nxx + jj] = Txz[ii*nxx + jj]*down_right[(ii-nzz+abc_layer)*abc_layer + (jj-nxx+abc_layer)];
        }
    }
}

float * factor_attenuation(float parameter, int n_bondary)
{
    int index;
    float * factor = (float *) malloc(n_bondary*sizeof(float));
    
    for(index = 0; index < n_bondary; ++index)
    {
        factor[index] = exp(-pow(parameter*(n_bondary-index),2));
    }

    return (factor);        
}
