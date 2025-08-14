# include "utils.hpp"
# include "cubic.hpp"
# include "linear.hpp"

void build_function(float * f, float * x, int nx)
{
    float xi =-10.0f;
    float xf = 10.0f;

    linspace(x, xi, xf, nx);

    for (int i = 0; i < nx; i++)
        f[i] = sinf(x[i]);
}

void perform_linear_interpolation(float * result, float * F, float * xp, float * x, int skipx, int nx)
{
    float P[2];

    for (int idx = skipx + 1; idx < nx - skipx - 1; idx++)
    {   
        int ipx = (int)(idx/skipx);

        float dx = (x[idx] - xp[ipx]) / (xp[ipx + 1] - xp[ipx]);  

        for (int pIdx = 0; pIdx < 2; pIdx++)
            P[pIdx] = F[ipx + pIdx];                

        result[idx] = linear1d(P, dx);   
    }
}

void perform_cubic_interpolation(float * result, float * F, float * xp, float * x, int skipx, int nx)
{
    float P[4];

    for (int idx = skipx + 1; idx < nx - skipx - 1; idx++)
    {   
        int ipx = (int)(idx/skipx);

        float dx = (x[idx] - xp[ipx]) / (xp[ipx + 1] - xp[ipx]);  
        
        for (int pIdx = 0; pIdx < 4; pIdx++)
            P[pIdx] = F[ipx + pIdx - 1];                

        result[idx] = cubic1d(P, dx);
    }
}

int main()
{
    int original_Nx = 1001;
    int rescaled_Nx = 51;

    float * original_x = new float[original_Nx]();
    float * rescaled_x = new float[rescaled_Nx]();

    float * original_F = new float[original_Nx]();
    float * rescaled_F = new float[rescaled_Nx]();

    float * result = new float[original_Nx]();

    build_function(original_F, original_x, original_Nx);
    export_array("original1d.bin", original_F, original_Nx);

    build_function(rescaled_F, rescaled_x, rescaled_Nx);
    export_array("rescaled1d.bin", rescaled_F, rescaled_Nx);

    int skipx = (int)((original_Nx - 1)/(rescaled_Nx - 1));

    auto tli = std::chrono::system_clock::now();    
    perform_linear_interpolation(result, rescaled_F, rescaled_x, original_x, skipx, original_Nx);
    auto tlf = std::chrono::system_clock::now();
    export_array("linear1d.bin", result, original_Nx);
    std::chrono::duration<double> tl = tlf - tli;
    std::cout << "1D linear interpolation runtime: " << tl.count() << " s." << std::endl;

    auto tci = std::chrono::system_clock::now();    
    perform_cubic_interpolation(result, rescaled_F, rescaled_x, original_x, skipx, original_Nx);
    auto tcf = std::chrono::system_clock::now();
    export_array("cubic1d.bin", result, original_Nx);
    std::chrono::duration<double> tc = tcf - tci;
    std::cout << "1D cubic interpolation runtime: " << tc.count() << " s." << std::endl;

    return 0;
}