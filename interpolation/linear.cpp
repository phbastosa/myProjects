# include "linear.hpp"

float linear1d(float P[2], float dx)
{
    return P[0] + dx*(P[1] - P[0]);
}

float linear2d(float P[2][2], float dx, float dy)
{    
    float p[2];
    p[0] = linear1d(P[0], dy);
    p[1] = linear1d(P[1], dy);
    return linear1d(p, dx);
}

float linear3d(float P[2][2][2], float dx, float dy, float dz)
{    
    float p[2];
    p[0] = linear2d(P[0], dy, dz);
    p[1] = linear2d(P[1], dy, dz);
    return linear1d(p, dx);
}
