// Single Source Shortest Path Algorithm testing Application //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //


#include <stdio.h>
#include "Graph/Graph.hpp"
#include "Graph/Converter.hpp"
#include "Graph/CSRGraph.hpp"
#include "Algorithm/APSP.hpp"
#include <time.h>

int main(){
    constexpr size_t N = 512;
    constexpr size_t E_Min = 2;
    constexpr size_t E_Max = 5;
    Graph g(N);
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < rand()%(E_Max - E_Min) + E_Min; ++j){
            size_t to = rand()%N;
            weight_t v = (rand()%1000)/100;
            g.connect(i, to, v);
            // printf("%d=>%d %f\n", i, to, v);
        }
    }
    CSRGraph csr = Converter::Instance().ToCSR(g);
    weight_t* result = nullptr;
    APSP(csr, &result);

    for(int i = 0; i < N ; ++i){
        for(int j = 0; j < N; ++j){
            printf("%f ", result[i * N + j]);
        }
    }
    printf("\n");
    delete[] result;
    return 0;
}