# include <stdio.h>
# include <stdlib.h>
# include "../header_files/io_functions.h"

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

void exporting_2D_snapshots(int time_pointer, float *snaps, float *P_pre, float *vp,
                         char snaps_file[], int nsrc, int nxx, int nzz)
{
    int ii, jj, index;

    /* Exporting the snapshots */
    if((time_pointer % 100 == 0) && (time_pointer > nsrc/2)) 
    {
        for(index = 0; index < nxx*nzz; ++index) 
        {                    
            /* Current calculated wave field whith the velocities model */
            snaps[index] = P_pre[index] + 0.000001*vp[index]; 
        }
        export_float32(snaps_file,nxx*nzz,snaps);            
    }
}
