# include <stdio.h>
# include <stdlib.h>
# include <math.h>

/* Function that read binary files. Made by: Felipe Capucci GISIS-UFF */
float * import_float32(char * file, int n_points) 
{
    FILE * fp;
	size_t result;

	fp = fopen(file, "rb");
	if (fp == NULL) {fputs ("File error\n", stderr); exit (1);}

	float * buffer = (float *) malloc(n_points*sizeof(float));

	if (buffer == NULL) {fputs ("Memory error\n",stderr); exit (2);}

	result = fread(buffer,sizeof(float),n_points,fp);
	if (result != n_points) {fputs ("Reading error\n",stderr); exit (3);}

	fclose (fp);
	return (buffer);
}

/* Function that write binary files. Made by: Felipe Capucci GISIS-UFF */
void export_float32(char * file, int n_points, float * vector) 
{
    FILE * fp;
	fp = fopen(file, "wb");
    	if(fp != NULL) fwrite((char *)vector, n_points*sizeof(float), 1, fp);
	fclose(fp);
}

void exporting_2D_snapshots(int time_pointer, float *snaps, float *WaveField, float *vp,
                         char snaps_file[], int nsrc, int nxx, int nzz)
{
    int ii, jj, index;

    /* Exporting the snapshots */
    if((time_pointer % 100 == 0) && (time_pointer > nsrc/2)) 
    {
        for(index = 0; index < nxx*nzz; ++index) 
        {                    
            /* Current calculated wave field whith the velocities model */
            snaps[index] = WaveField[index] + 0.000001*vp[index]; 
        }
        
        FILE *snap = fopen(snaps_file,"ab");
        for (ii = 0; ii < nxx; ii++) {
            for(jj = 0; jj < nzz; jj++) {
                fwrite(&WaveField[jj*nxx + ii],sizeof(float),1,snap);        
            }
        }
        fclose(snap);            
    }
}

void exporting_pointer_seismogram(char * file_name, int nx, int nt, float * seismogram)
{
    int ii,jj;
    
    FILE * seism = fopen(file_name,"wb");
    for (ii = 0; ii < nx; ii++) 
    {
        for (jj = 0; jj < nt; jj++) 
        {
            fwrite(&seismogram[jj*nx + ii],sizeof(float),1,seism);
        }
    }
    fclose(seism);
}

int describe_2D_model_stability(float * vp, float * vs, int nx, int nz, float f_cut, int nt, float dt, float ds)
{
    int index,ii,jj;
    float v_max,v_min;
    float pi = 4*atan(1);
    float sug_ds, sug_dt;

    for(index = 0; index < nx*nz; ++index) 
    {
        ii = floor(index / nx);  /* line indicator */
        jj = index % nx;         /* column indicator */  

        /* Finding the highest velocity in the model */
        v_max = vp[0];
        if(v_max < vp[ii*nx + jj]) v_max = vp[ii*nx + jj]; 
        
        /* Finding the slowest velocity in the model */
        v_min = vs[0];
        if(v_min > vs[ii*nx + jj]) v_min = vs[ii*nx + jj]; 
    }

    /* Nyquist limit numerical stability */
    if(((v_max*dt)/ds <= sqrt(4/(3*powf(pi,2)))) && (ds <= v_min/(3.3*f_cut)) && (dt <= (v_min/(3.3*f_cut))/(4*v_max)))  
    {    
        printf("\nThe modeling will be stable according to Nyquist limit numerical stability!\n");        
    } 
    else 
    { 
        printf("\nThe modeling won't be stable according to Nyquist limit numerical stability!\n\n"); 
        sug_ds = v_min/(3.3*f_cut);
        sug_dt = (v_min/(3.3*f_cut))/(4*v_max); 
        printf("Sugestions: \n ds <= %f \n dt <= %f\n\n",sug_ds,sug_dt);
        return 1;       
    }

    printf("\nParameters:\n");
    printf("   Highest velocity in the model = %.1f m/s\n",v_max);
    printf("   Slowest velocity in the model = %.1f m/s\n",v_min);
    printf("   Spatial discratization = %.1f meters\n",ds);
    printf("   Temporal discratization = %.4f seconds\n",dt);
    printf("   Source cutoff frequency = %.1f Hz\n",f_cut);
    printf("   Total modeling time = %.2f seconds\n",dt*nt);
    printf("   Horizontal length of model = %.0f meters\n",ds*nx);
    printf("   Vertical length of model = %.0f meters\n\n",ds*nz);
}