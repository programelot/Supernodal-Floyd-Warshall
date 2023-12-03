// Single Source Shortest Path Algorithm testing Application //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include <stdio.h>
#include "Graph/Graph.hpp"
#include "Graph/Converter.hpp"
#include "Graph/CSRGraph.hpp"
#include "Matrix/MtxReader.hpp"
#include "Algorithm/SSSP.hpp"
#include <time.h>

int main(int argc, char* argv[]){
    if(argc != 2){
        printf("Usage : APSP MatrixFile \n");
        return 0;
    }
    clock_t start, finish;

    Graph g = MtxReader::Instance().Read(argv[1]);
    CSRGraph csr = Converter::Instance().ToCSR(g);
    printf("Graph prepared\n");
    weight_t* result = nullptr;
    start = clock();
    SSSP(0, csr, &result);
    finish = clock();

    printf("%fsec\n", (double)(finish - start)/CLOCKS_PER_SEC);

    // for(int i = 0; i < N ; ++i){
    //     printf("%f ", result[i]);
    // }
    printf("\n");
    delete[] result;
    return 0;
}