# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include "acoustic1D.h"

int main(int argc, char **argv) 
{
    float totalTime;  /* Total time */
    time_t ti, tf;    /* Initial and final time */
    ti = time(NULL);  /* Starting the time */

    int nx = 500;
    int nt = 1001;         
    int sx = 250;      /* Shot position in the modeling */   
    int nsrc = 500;    /* Total samples of the source */
    int nabc = 100;    /* Samples used in the attenuation layer */
    
    float dx = 10.0f;
    float dt = 0.001f;
    float fcut = 30.0f; 

    float * pas = (float *) malloc(nx*sizeof(float));   /* Future wavefield */
    float * pre = (float *) malloc(nx*sizeof(float));   /* Present wavefield */
    float * fut = (float *) malloc(nx*sizeof(float));   /* Past wavefield */
    float * vel = (float *) malloc(nx*sizeof(float));   /* Compressional wave velocities */ 

    for (int ii = 0; ii < nx; ii++) vel[ii] = 2000.0f;

    setToZero(pas,pre,fut,nx);

    for(int kk = 0; kk < nt; kk++) 
    {
        if(kk % 200 == 0) printf("Tempo = %.2f segundos\n",dt*k);

        if(kk < nsrc) pre[sx] =  /* applying the source in modeling */ 

        acoustic1D_8E2T(pas,pre,fut,vel,nx,dx,dt);    

        updateWavefield(pas,pre,fut,nx);
    }

    tf = time(NULL);                /* Initializing the final time */
    totalTime = difftime(tf, ti);   /* Calculating the difference between the initial and final times */
    printf("\nExecution time: \033[31m%.0fs\n\033[m", totalTime);

    return 0;
}
