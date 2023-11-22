#include <stdio.h>
#include <papi.h>

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
static void HandleError( cudaError_t err,
    const char *file,
    int line ) {
    if (err != cudaSuccess) {
    printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
    file, line );
    }
}

__device__ unsigned int decomposeCUDA(int nVertex, int i, int j = 0){
    return i + j * nVertex;
}

__global__ void updateEdge(int threadSize, int nVertex, unsigned int* apsp){
    int idx = threadIdx.x + blockIdx.x * threadSize;
    if(idx > nVertex * nVertex)
        return;
    int i = idx % nVertex; 
    int j = (idx / nVertex) % nVertex;
    for(int k = 0; k < nVertex; ++k){
        if(apsp[decomposeCUDA(nVertex,i,j)] > apsp[decomposeCUDA(nVertex,i,k)] + apsp[decomposeCUDA(nVertex,k,j)]){
            apsp[decomposeCUDA(nVertex,i,j)] = apsp[decomposeCUDA(nVertex,i,k)] + apsp[decomposeCUDA(nVertex,k,j)];
        }
    }
}

unsigned int decomposeGPU(int nVertex, int i, int j = 0){
    return i + j * nVertex;
}

void pathDouble(unsigned int nVertex, unsigned int* apsp){
    unsigned long long allFrom = PAPI_get_real_usec();
    unsigned int* apsp_cuda;

    HANDLE_ERROR(cudaMalloc((void**)&apsp_cuda, sizeof(unsigned int) * nVertex * nVertex));
    
    HANDLE_ERROR(cudaMemcpy(apsp_cuda, apsp, sizeof(unsigned int) * nVertex * nVertex, cudaMemcpyHostToDevice));
    
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    unsigned long long computingSize = nVertex * nVertex ;
    unsigned long long threadSize = prop.maxThreadsPerBlock;
    unsigned long long blockSize = (computingSize + threadSize - 1)/threadSize;
    if(threadSize > computingSize){
        threadSize = computingSize;
    }
    int time = 0;
    int currentDist = 1;
    while(currentDist < nVertex - 1){
        ++time;
        currentDist *= 2;
    }
    unsigned long long allTo = PAPI_get_real_usec();
    for(int i = 0; i < time; ++i)
        updateEdge<<<blockSize, threadSize>>>(threadSize, nVertex, apsp_cuda);

    unsigned long long freeFrom = PAPI_get_real_usec();
    HANDLE_ERROR(cudaMemcpy(apsp, apsp_cuda, sizeof(unsigned int) * nVertex * nVertex, cudaMemcpyDeviceToHost));
    
    HANDLE_ERROR(cudaFree(apsp_cuda));
    
    unsigned long long freeTo = PAPI_get_real_usec();
    printf("%d %d %d %d %lu %lu\n", threadSize, blockSize, computingSize, time, allTo - allFrom, freeTo - freeFrom);
}