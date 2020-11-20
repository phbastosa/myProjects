# include <cuda.h>
# include <stdio.h>
# include <stdlib.h>
# include <cuda_runtime.h>

__global__ void addInts(int * a, int * b, int n)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(index < n) a[index] += b[index];
}

int main(int argc, char **argv)
{
    int n = 100;

    int * h_a = (int *) malloc(n*sizeof(int)); 
    int * h_b = (int *) malloc(n*sizeof(int)); 

    for (int i = 0; i < n; i++)
    {
        h_a[i] = rand() % 1000;
        h_b[i] = rand() % 1000;
    }

    int * d_a;
    int * d_b;

    cudaMalloc(&d_a,n*sizeof(int));
    cudaMalloc(&d_b,n*sizeof(int));
    
    cudaMemcpy(d_a,h_a,n*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_b,h_b,n*sizeof(int),cudaMemcpyHostToDevice);

    addInts<<<1,n>>>(d_a,d_b,n);

    cudaMemcpy(h_a,d_a,n*sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(h_b,d_b,n*sizeof(int),cudaMemcpyDeviceToHost);

    printf("Results of summation:\n");
    for (int i = n-10; i < n; i++) 
        printf("%d + %d = %d\n",h_a[i]-h_b[i],h_b[i],h_a[i]);

    delete[] h_a;
    delete[] h_b;

    cudaFree(d_a);
    cudaFree(d_b);

    return 0;
}