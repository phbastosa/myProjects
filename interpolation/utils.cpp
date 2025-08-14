# include "utils.hpp" 

void linspace(float * x, float xi, float xf, int n)
{
    float spacing = (xf - xi)/(n - 1);
    for (int i = 0; i < n; i++)
        x[i] = xi + i*spacing;
}

void export_array(std::string path, float * array, int n)
{
    std::ofstream file(path, std::ios::out);
    file.write((char *) array, n*sizeof(float));
    file.close();
}
