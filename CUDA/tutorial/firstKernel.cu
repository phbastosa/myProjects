# include <cuda.h>
# include <stdio.h>
# include <cuda_runtime.h>

__global__ void addInt(int *a, int *b)
{
    a[0] += b[0];
}

int main(int argc, char **argv)
{
    int a = 5;
    int b = 9;

    int *d_a, *d_b;

    cudaMalloc(&d_a,sizeof(int));
    cudaMalloc(&d_b,sizeof(int));

    cudaMemcpy(d_a,&a,sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_b,&b,sizeof(int),cudaMemcpyHostToDevice);

    addInt<<<1,1>>>(d_a,d_b);

    cudaMemcpy(&a,d_a,sizeof(int),cudaMemcpyDeviceToHost);

    cudaFree(d_a);
    cudaFree(d_b);

    printf("The awnser is %i\n",a);

    return 0;
}