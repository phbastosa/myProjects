# include "utils.hpp"
# include "cubic.hpp"
# include "linear.hpp"

void build_function(float * f, float * x, float * y, float * z, int nx, int ny, int nz)
{
    float xi =-10.0f;
    float yi =-10.0f;
    float zi =-10.0f;

    float xf = 10.0f;
    float yf = 10.0f;
    float zf = 10.0f;

    linspace(x, xi, xf, nx);
    linspace(y, yi, yf, ny);
    linspace(z, zi, zf, nz);

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
                f[i + j*nx + k*nx*ny] = sinf(x[i]) + sinf(y[j]) + sinf(z[k]);
        }
    }
}

void perform_linear_interpolation(float * result, float * F, float * xp, float * yp, float * zp, float * x, float * y, float * z, int skipx, int skipy, int skipz, int nx, int ny, int nz)
{
    float P[2][2][2];

    int rNx = (int)((nx-1)/skipx) + 1;
    int rNy = (int)((nz-1)/skipy) + 1;

    for (int idx = skipx + 1; idx < nx - skipx - 1; idx++)
    {   
        for (int idy = skipy + 1; idy < ny - skipy - 1; idy++)
        {   
            for (int idz = skipz + 1; idz < nz - skipz - 1; idz++)
            {   
                int ipx = (int)(idx/skipx);
                int ipy = (int)(idy/skipy);
                int ipz = (int)(idz/skipz);

                float dx = (x[idx] - xp[ipx]) / (xp[ipx+1] - xp[ipx]);  
                float dy = (y[idy] - yp[ipy]) / (yp[ipy+1] - yp[ipy]);  
                float dz = (z[idz] - zp[ipz]) / (zp[ipz+1] - zp[ipz]);  

                for (int pIdx = 0; pIdx < 2; pIdx++)
                {
                    for (int pIdy = 0; pIdy < 2; pIdy++)
                    {
                        for (int pIdz = 0; pIdz < 2; pIdz++)
                            P[pIdx][pIdy][pIdz] = F[(ipx + pIdx) + (ipy + pIdy)*rNx + (ipz + pIdz)*rNx*rNy];                
                    }
                }

                result[idx + idy*nx + idz*nx*ny] = linear3d(P, dx, dy, dz);   
            }
        }
    }
}

void perform_cubic_interpolation(float * result, float * F, float * xp, float * yp, float * zp, float * x, float * y, float * z, int skipx, int skipy, int skipz, int nx, int ny, int nz)
{
    float P[4][4][4];

    int rNx = (int)((nx-1)/skipx) + 1;
    int rNy = (int)((nz-1)/skipz) + 1;

    for (int idx = skipx + 1; idx < nx - skipx - 1; idx++)
    {   
        for (int idy = skipy + 1; idy < ny - skipy - 1; idy++)
        {   
            for (int idz = skipz + 1; idz < nz - skipz - 1; idz++)
            {   
                int ipx = (int)(idx/skipx);
                int ipy = (int)(idy/skipy);
                int ipz = (int)(idz/skipz);

                float dx = (x[idx] - xp[ipx]) / (xp[ipx+1] - xp[ipx]);  
                float dy = (y[idy] - yp[ipy]) / (yp[ipy+1] - yp[ipy]);  
                float dz = (z[idz] - zp[ipz]) / (zp[ipz+1] - zp[ipz]);  

                for (int pIdx = 0; pIdx < 4; pIdx++)
                {
                    for (int pIdy = 0; pIdy < 4; pIdy++)
                    {
                        for (int pIdz = 0; pIdz < 4; pIdz++)
                            P[pIdx][pIdy][pIdz] = F[(ipx + pIdx - 1) + (ipy + pIdy - 1)*rNx + (ipz + pIdz - 1)*rNx*rNy];                
                    }
                }
                
                result[idx + idy*nx + idz*nx*ny] = cubic3d(P, dx, dy, dz);   
            }
        }
    }
}

int main()
{
    int original_Nx = 1001;
    int original_Ny = 1001;
    int original_Nz = 1001;
    
    int rescaled_Nx = 51;
    int rescaled_Ny = 51;
    int rescaled_Nz = 51;

    int original_N = original_Nx*original_Ny*original_Nz;
    int rescaled_N = rescaled_Nx*rescaled_Ny*rescaled_Nz;

    float * original_x = new float[original_Nx]();
    float * original_y = new float[original_Nx]();
    float * original_z = new float[original_Nx]();

    float * rescaled_x = new float[rescaled_Nx]();
    float * rescaled_y = new float[rescaled_Nx]();
    float * rescaled_z = new float[rescaled_Nx]();

    float * original_F = new float[original_N]();
    float * rescaled_F = new float[rescaled_N]();

    float * result = new float[original_N]();

    build_function(original_F, original_x, original_y, original_z, original_Nx, original_Ny, original_Nz);
    export_array("original3d.bin", original_F, original_N);

    build_function(rescaled_F, rescaled_x, rescaled_y, rescaled_z, rescaled_Nx, rescaled_Ny, rescaled_Nz);
    export_array("rescaled3d.bin", rescaled_F, rescaled_N);

    int skipx = (int)((original_Nx - 1)/(rescaled_Nx - 1));
    int skipy = (int)((original_Ny - 1)/(rescaled_Ny - 1));
    int skipz = (int)((original_Nz - 1)/(rescaled_Nz - 1));

    auto tli = std::chrono::system_clock::now();    
    perform_linear_interpolation(result, rescaled_F, rescaled_x, rescaled_y, rescaled_z, original_x, original_y, original_z, skipx, skipy, skipz, original_Nx, original_Ny, original_Nz);
    auto tlf = std::chrono::system_clock::now();
    export_array("linear3d.bin", result, original_N);
    std::chrono::duration<double> tl = tlf - tli;
    std::cout << "3D linear interpolation runtime: " << tl.count() << " s." << std::endl;

    auto tci = std::chrono::system_clock::now();    
    perform_cubic_interpolation(result, rescaled_F, rescaled_x, rescaled_y, rescaled_z, original_x, original_y, original_z, skipx, skipy, skipz, original_Nx, original_Ny, original_Nz);
    auto tcf = std::chrono::system_clock::now();
    export_array("cubic3d.bin", result, original_N);
    std::chrono::duration<double> tc = tcf - tci;
    std::cout << "3D cubic interpolation runtime: " << tc.count() << " s." << std::endl;

    return 0;
}