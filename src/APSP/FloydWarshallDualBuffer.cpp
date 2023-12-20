// FloydWarshall algorithm library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

// Parallelism : O(|V|^2) //
// Work complexity : O(|V|^3) //
// Time complexity with PRAM : O(|V|) //

#include "Common/Type.hpp"
#include "Graph/CSRGraph.hpp"
#include "Algorithm/APSP.hpp"
#include <omp.h>

void APSP(const CSRGraph& input_graph, weight_t* distance){
    
    dataSize_t size = input_graph.Size();
    dataSize_t* rowPtr = input_graph.RowPtr();
    dataSize_t* colIdx = input_graph.ColIdx();
    weight_t* value = input_graph.Value();
    weight_t* result[2] = {new weight_t[size * size], new weight_t[size * size]};

    for(dataSize_t i = 0; i < size ; ++i){
        for(dataSize_t j = 0; j < size; ++j){
            result[0][i * size + j] = kWeightInf;
        }
    }

    for(dataSize_t i = 0; i < size; ++i){
        result[0][i * size + i] = 0;
    }

    for(dataSize_t i = 0; i < size; ++i){
        for(dataSize_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
            result[0][i * size + colIdx[j]] = value[j];
        }
    }
    
    // int nested = omp_get_nested();
    // omp_set_nested(1);
    for(dataSize_t k = 0; k < size; ++k){
        #pragma omp parallel
        {
            #pragma omp for
            for(dataSize_t i = 0; i < size; ++i){
                // #pragma omp parallel
                // {
                //     #pragma omp for
                        for(dataSize_t j = 0; j < size; ++j){
                            weight_t oldValue = result[k % 2][i * size + j];
                            weight_t newValue = result[k % 2][i * size + k] + result[k % 2][k * size + j];
                            result[(k + 1) % 2][i * size + j] = (oldValue > newValue) ? newValue : oldValue;
                        }
                // }
            }
        }
    }
    // omp_set_nested(nested);

    for(dataSize_t i = 0; i < size; ++i){
        for(dataSize_t j = 0; j < size; ++j){
            distance[i * size + j] = result[size % 2][i * size + j];
        }
    }
    delete[] result[(size + 1) % 2];
}
