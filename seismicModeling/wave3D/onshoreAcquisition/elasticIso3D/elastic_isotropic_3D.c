# include <time.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include "../../auxiliaryCodes/generalFunctions/functions.h"

int main(int argc, char **argv) 
{
    float total_time;
    time_t t_0, t_f;
    t_0 = time(NULL);

    int nx = 350;        /* Tamanho do modelo na direção de x */
    int ny = 350;        /* Tamanho do modelo na direção de y */
    int nz = 400;        /* Tamanho do modelo na direção de z */
    int nt = 1501;       /* Total de amostras no tempo */

    int borda = 50;

    float dx = 15.0f;   /* Discretização espacial da malha em x */
    float dy = 15.0f;   /* Discretização espacial da malha em x */ 
    float dz = 10.0f;   /* Discretização espacial da malha em z */
    float dt = 0.002f;  /* Discretização temporal da malha */
    
    int * sx = (int *) malloc(sizeof(int));     /* Posição da fonte em x no modelo */
    int * sy = (int *) malloc(sizeof(int));     /* Posição da fonte em y no modelo */
    int * sz = (int *) malloc(sizeof(int));     /* Posição da fonte em z no modelo */
    
    sx[0] = (int) nx/2;
    sy[0] = (int) ny/2;    
    sz[0] = (int) 2*borda;    

    int nsrc = 300;     /* Total de amostras em tempo da fonte */
    float freq = 30.0f; /* Frequência de corte da fonte */
            
    float *vp = (float *) malloc(nx*ny*nz*sizeof(float));
    float *vs = (float *) malloc(nx*ny*nz*sizeof(float));
    float *rho = (float *) malloc(nx*ny*nz*sizeof(float));
    
    float *M = (float *) malloc(nx*ny*nz*sizeof(float));
    float *L = (float *) malloc(nx*ny*nz*sizeof(float));
        
    float *Txx = (float *) malloc(nx*ny*nz*sizeof(float));
    float *Tyy = (float *) malloc(nx*ny*nz*sizeof(float));
    float *Tzz = (float *) malloc(nx*ny*nz*sizeof(float));
    float *Txz = (float *) malloc(nx*ny*nz*sizeof(float));
    float *Tyz = (float *) malloc(nx*ny*nz*sizeof(float));
    float *Txy = (float *) malloc(nx*ny*nz*sizeof(float));
    
    float *Vx = (float *) malloc(nx*ny*nz*sizeof(float));
    float *Vy = (float *) malloc(nx*ny*nz*sizeof(float));
    float *Vz = (float *) malloc(nx*ny*nz*sizeof(float));

    int ii,jj,kk,index;    
    for(index = 0; index < nx*ny*nz; index++) 
    {
        kk = floor(index/(nx*nz));        // matrizes (direção y)
        jj = index % nx;                  // colunas  (direção x)
        ii = floor((index % (nx*nz))/nx); // linhas   (direção z)  

        if (ii < 150)
        {
            vp[index] = 1500.0f;
            vs[index] = 0.0f;
            rho[index] = 1000.0f;
        }
        else if (ii < 200)
        {
            vp[index] = 1700.0f;
            vs[index] = vp[index]/1.7;
            rho[index] = 310*powf(vp[index],0.25);
        }
        else if (ii < 250)
        {
            vp[index] = 1900.0f;
            vs[index] = vp[index]/1.7;
            rho[index] = 310*powf(vp[index],0.25);
        }
        else if (ii < 300)
        {
            vp[index] = 2300.0f;
            vs[index] = vp[index]/1.7;
            rho[index] = 310*powf(vp[index],0.25);
        }        
        else 
        {
            vp[index] = 2500.0f;
            vs[index] = vp[index]/1.7;
            rho[index] = 310*powf(vp[index],0.25);
        }

        // vp[index] = 1500.0f;
        // vs[index] = 0.0f;
        // rho[index] = 1000.0f;
        M[index] = rho[index]*pow(vs[index],2.0f);
        L[index] = rho[index]*pow(vp[index],2.0f) - 2.0f*M[index];
    }

    float * seismVx = (float *) malloc(nx*nt*sizeof(float));
    float * seismVy = (float *) malloc(nx*nt*sizeof(float));
    float * seismVz = (float *) malloc(nx*nt*sizeof(float));
    float * seismPs = (float *) malloc(nx*nt*sizeof(float));

    float * pressure = (float *) malloc(nx*nz*sizeof(float));

    float * source = generateStaggeredSource(freq,dt,nsrc);

    // Cerjan boundary condition
    float * factor = factorAttenuation(borda,0.0055f);
    float * prismX = (float *) malloc(borda*borda*nx*sizeof(float));
    float * prismY = (float *) malloc(borda*borda*ny*sizeof(float)); 
    float * prismZ = (float *) malloc(borda*borda*nz*sizeof(float));
    float * cubXYZ = (float *) malloc(borda*borda*borda*sizeof(float));
    cerjanCorners3D(prismX,prismY,prismZ,cubXYZ,factor,borda,nx,ny,nz);

    memSet(Txx,nx*ny*nz);
    memSet(Tyy,nx*ny*nz);
    memSet(Tzz,nx*ny*nz);
    memSet(Txz,nx*ny*nz);
    memSet(Tyz,nx*ny*nz);
    memSet(Txy,nx*ny*nz);
   
    memSet(Vx,nx*ny*nz);
    memSet(Vy,nx*ny*nz);
    memSet(Vz,nx*ny*nz);

    FILE * ps = fopen("snapPs.bin","ab");
    
    FILE * vx = fopen("snapVx.bin","ab");
    FILE * vy = fopen("snapVy.bin","ab");
    FILE * vz = fopen("snapVz.bin","ab");

    for (kk = 0; kk < nt; kk++) 
    {
        if(kk % (nt/10) == 0)
            printf("Tempo = %.2f segundos\n",dt*kk);

        if(kk < nsrc)
        {
            Txx[sy[0]*nx*nz + sz[0]*nx + sx[0]] += source[kk];
            Tyy[sy[0]*nx*nz + sz[0]*nx + sx[0]] += source[kk];
            Tzz[sy[0]*nx*nz + sz[0]*nx + sx[0]] += source[kk];
        }

        // FDM_8E2T_PML_elastic_Isotropic3D(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,rho,M,L,nx,ny,nz,dt,dx,dy,dz);
        FDM_8E2T_elastic_Isotropic3D(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,rho,M,L,nx,ny,nz,dt,dx,dy,dz);
    
        cerjanElasticAbsorbingCondition3D(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,factor,prismX,prismY,prismZ,cubXYZ,nx,ny,nz,borda);

        for (index = 0; index < nx*nz; index++)
        {
            jj = floor(index / nx);  /* Line indicator */
            ii = index % nx;         /* Column indicator */  

            pressure[jj*nx + ii] = (Txx[sy[0]*nx*nz + jj*nx + ii] + 
                                    Tyy[sy[0]*nx*nz + jj*nx + ii] + 
                                    Tzz[sy[0]*nx*nz + jj*nx + ii])/3; 
        }

        // Snapshots
        if ((kk % 10 == 0) && (kk >= nsrc/2))
        {
            for (index = 0; index < nx*nz; index++)
            {
                jj = floor(index / nx);  /* Line indicator */
                ii = index % nx;         /* Column indicator */  

                // Plano XY
                // fwrite(&Vx[jj*nx*nz + sz[0]*nx + ii],sizeof(float),1,vx);
                // fwrite(&Vy[jj*nx*nz + sz[0]*nx + ii],sizeof(float),1,vy);
                // fwrite(&Vz[jj*nx*nz + sz[0]*nx + ii],sizeof(float),1,vz);
                
                // Plano YZ
                // fwrite(&Vx[ii*nx*nz + jj*nx + sx[0]],sizeof(float),1,vx);
                // fwrite(&Vy[ii*nx*nz + jj*nx + sx[0]],sizeof(float),1,vy);
                // fwrite(&Vz[ii*nx*nz + jj*nx + sx[0]],sizeof(float),1,vz);
                
                // Plano XZ
                fwrite(&Vx[sy[0]*nx*nz + jj*nx + ii],sizeof(float),1,vx);
                fwrite(&Vy[sy[0]*nx*nz + jj*nx + ii],sizeof(float),1,vy);
                fwrite(&Vz[sy[0]*nx*nz + jj*nx + ii],sizeof(float),1,vz);
                fwrite(&pressure[jj*nx + ii],sizeof(float),1,ps);
            }
        }
 
        // Sismograma plano XZ
        for (index = 0; index < nx; index++)
        {
            seismVx[kk*nx + index] = Vx[sy[0]*nx*nz + 2*borda*nx + index]; 
            seismVy[kk*nx + index] = Vy[sy[0]*nx*nz + 2*borda*nx + index]; 
            seismVz[kk*nx + index] = Vz[sy[0]*nx*nz + 2*borda*nx + index]; 
            seismPs[kk*nx + index] = pressure[(2*borda)*nx + index];
        }
    }
    
    exportVector(seismVx,nx*nt,"seismVx.bin");    
    exportVector(seismVy,nx*nt,"seismVy.bin");    
    exportVector(seismVz,nx*nt,"seismVz.bin");    
    exportVector(seismPs,nx*nt,"seismPs.bin");    

    t_f = time(NULL);
    total_time = difftime(t_f, t_0);
    printf("\nExecution time: \033[31m%.0fs\n\033[m", total_time);    
    return 0;
} 