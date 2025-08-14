# include "cubic.hpp"

float cubic1d(float P[4], float dx)
{
    return P[1] + 0.5f*dx*(P[2] - P[0] + dx*(2.0f*P[0] - 5.0f*P[1] + 4.0f*P[2] - P[3] + dx*(3.0f*(P[1] - P[2]) + P[3] - P[0])));
}

float cubic2d(float P[4][4], float dx, float dy)
{    
    float p[4];
    p[0] = cubic1d(P[0], dy);
    p[1] = cubic1d(P[1], dy);
    p[2] = cubic1d(P[2], dy);
    p[3] = cubic1d(P[3], dy);    
    return cubic1d(p, dx);
}

float cubic3d(float P[4][4][4], float dx, float dy, float dz)
{    
    float p[4];
    p[0] = cubic2d(P[0], dy, dz);
    p[1] = cubic2d(P[1], dy, dz);
    p[2] = cubic2d(P[2], dy, dz);
    p[3] = cubic2d(P[3], dy, dz);
    return cubic1d(p, dx);
}
