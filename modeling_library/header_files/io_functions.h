#ifndef IO_FUNCTIONS_H_DEFINED
#define IO_FUNCTIONS_H_DEFINED

float * import_float32(char * file, int n_points); 
void export_float32(char * file, int n_points, float * vector); 
void exporting_2D_snapshots(int time_pointer, float *snaps, float *P_pre, float *vp,
                            char snaps_file[], int nsrc, int nxx, int nzz);

#endif // IO_FUNCTIONS_H_DEFINED