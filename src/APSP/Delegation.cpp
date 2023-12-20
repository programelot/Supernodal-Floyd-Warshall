// Deligated APSP algorithm library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

// Parallelism : O(|V|) //

// Work complexity : O(|V|) * Delegated algorithm //
// Example 
/// Bellman-ford    : O(|E||V|^2)       //
/// Dijkstra        : O(|E||V|log|V|)   //

// Time complexity with PRAM : The same with Delegated algorith //
// Example 
/// Bellman-ford    : O(|E||V|)       //
/// Dijkstra        : O(|E|log|V|)   //

#include "Common/Type.hpp"
#include "Graph/CSRGraph.hpp"
#include "Algorithm/APSP.hpp"
#include "Algorithm/SSSP.hpp"
#include <vector>
#include <omp.h>

void APSP(const CSRGraph& input_graph, weight_t* distance){
    dataSize_t size = input_graph.Size();
    #pragma omp parallel
    {
        #pragma omp for
        for(dataSize_t i = 0; i < size; ++i){
            weight_t* SSSP_result = new weight_t[size];
            SSSP(i, input_graph, SSSP_result);
            for(dataSize_t j = 0; j < size; ++j){
                distance[size * i + j] = SSSP_result[j];
            }
            delete[] SSSP_result;
        }
    }
}