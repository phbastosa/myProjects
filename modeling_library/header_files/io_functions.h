#ifndef IO_FUNCTIONS_H_DEFINED
#define IO_FUNCTIONS_H_DEFINED

float * import_float32(char * file, int n_points); 
void export_float32(char * file, int n_points, float * vector); 
void exporting_2D_snapshots(int time_pointer, float *snaps, float *WaveField, float *vp,
                            char snaps_file[], int nsrc, int nxx, int nzz);
void exporting_pointer_seismogram(char * file_name, int nx, int nt, float * seismogram);
int describe_2D_model_stability(float * vp, float * vs, int nx, int nz, float f_cut, int nt, float dt, float ds);

#endif // IO_FUNCTIONS_H_DEFINED