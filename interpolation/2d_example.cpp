# include "utils.hpp"
# include "cubic.hpp"
# include "linear.hpp"

void build_function(float * f, float * x, float * y, int Nx, int Ny)
{
    float xi =-10.0f;
    float xf = 10.0f;
    float yi =-10.0f;
    float yf = 10.0f;

    linspace(x, xi, xf, Nx);
    linspace(y, yi, yf, Ny);

    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
            f[i + j*Nx] = sinf(x[i]) + sinf(y[j]);
    }
}

void perform_linear_interpolation(float * result, float * F, float * xp, float * yp, float * x, float * y, int skipx, int skipy, int nx, int ny)
{
    float P[2][2];

    int rNx = (int)((nx - 1)/skipx) + 1;

    for (int idx = skipx + 1; idx < nx - skipx - 1; idx++)
    {   
        for (int idy = skipy + 1; idy < ny - skipy - 1; idy++)
        {   
            int ipx = (int)(idx/skipx);
            int ipy = (int)(idy/skipy);

            float dx = (x[idx] - xp[ipx]) / (xp[ipx+1] - xp[ipx]);  
            float dy = (y[idy] - yp[ipy]) / (yp[ipy+1] - yp[ipy]);  

            for (int pIdx = 0; pIdx < 2; pIdx++)
            {
                for (int pIdy = 0; pIdy < 2; pIdy++)    
                    P[pIdx][pIdy] = F[(ipx + pIdx) + (ipy + pIdy)*rNx];                
            }
            
            result[idx + idy*nx] = linear2d(P, dx, dy);   
        }
    }
}

void perform_cubic_interpolation(float * result, float * F, float * xp, float * yp, float * x, float * y, int skipx, int skipy, int nx, int ny)
{
    float P[4][4];

    int rNx = (int)((nx - 1)/skipy) + 1;

    for (int idx = skipx + 1; idx < nx - skipx - 1; idx++)
    {   
        for (int idy = skipy + 1; idy < ny - skipy - 1; idy++)
        {   
            int ipx = (int)(idx/skipx);
            int ipy = (int)(idy/skipy);

            float dx = (x[idx] - xp[ipx]) / (xp[ipx+1] - xp[ipx]);  
            float dy = (y[idy] - yp[ipy]) / (yp[ipy+1] - yp[ipy]);  

            for (int pIdx = 0; pIdx < 4; pIdx++)
            {
                for (int pIdy = 0; pIdy < 4; pIdy++)    
                    P[pIdx][pIdy] = F[(ipx + pIdx - 1) + (ipy + pIdy - 1)*rNx];                
            }

            result[idx + idy*nx] = cubic2d(P, dx, dy);   
        }
    }
}

int main()
{
    int original_Nx = 1001;
    int original_Ny = 1001;
    
    int rescaled_Nx = 51;
    int rescaled_Ny = 51;

    int original_N = original_Nx*original_Ny;
    int rescaled_N = rescaled_Nx*rescaled_Ny;

    float * original_x = new float[original_N]();
    float * original_y = new float[original_N]();

    float * rescaled_x = new float[rescaled_N]();
    float * rescaled_y = new float[rescaled_N]();

    float * original_F = new float[original_N]();
    float * rescaled_F = new float[rescaled_N]();

    float * result = new float[original_N]();

    build_function(original_F, original_x, original_y, original_Nx, original_Ny);
    export_array("original2d.bin", original_F, original_N);

    build_function(rescaled_F, rescaled_x, rescaled_y, rescaled_Nx, rescaled_Ny);
    export_array("rescaled2d.bin", rescaled_F, rescaled_N);

    int skipx = (int)((original_Nx - 1)/(rescaled_Nx - 1));
    int skipy = (int)((original_Ny - 1)/(rescaled_Ny - 1));

    auto tli = std::chrono::system_clock::now();    
    perform_linear_interpolation(result, rescaled_F, rescaled_x, rescaled_y, original_x, original_y, skipx, skipy, original_Nx, original_Ny);
    auto tlf = std::chrono::system_clock::now();
    export_array("linear2d.bin", result, original_N);
    std::chrono::duration<double> tl = tlf - tli;
    std::cout << "2D linear interpolation runtime: " << tl.count() << " s." << std::endl;

    auto tci = std::chrono::system_clock::now();    
    perform_cubic_interpolation(result, rescaled_F, rescaled_x, rescaled_y, original_x, original_y, skipx, skipy, original_Nx, original_Ny);
    auto tcf = std::chrono::system_clock::now();
    export_array("cubic2d.bin", result, original_N);
    std::chrono::duration<double> tc = tcf - tci;
    std::cout << "2D cubic interpolation runtime: " << tc.count() << " s." << std::endl;

    return 0;
}