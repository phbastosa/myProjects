# include <stdio.h>
# include <stdlib.h>

int main(int argc, char **argv) 
{
    int ii,jj;
    int nz, nx;
    
    nx = atoi(argv[2]);
    nz = atoi(argv[3]);

    float (*vp)[nx] = malloc(sizeof(float[nz][nx]));
    float (*va)[nx] = malloc(sizeof(float[nz][nx]));

    FILE *read = fopen(argv[1],"rb");
    for(ii = 0; ii < nx; ii++)
    {
        for(jj = 0; jj < nz; jj++)
        {
            fread(&vp[jj][ii], sizeof(float), 1, read);
        }
    }
    fclose(read);

    for(ii = 0; ii < nz; ii++)
    {
        for(jj = 0; jj < nx; jj++)
        {
            va[ii][jj] = 1/vp[ii][jj];
        }
    }

    FILE *write = fopen(argv[4],"wb");
    for(ii = 0; ii < nx; ii++)
    {
        for(jj = 0; jj < nz; jj++)
        {
            fwrite(&va[jj][ii], sizeof(float), 1, write);
        }
    }
    fclose(write);

    return 0;
}