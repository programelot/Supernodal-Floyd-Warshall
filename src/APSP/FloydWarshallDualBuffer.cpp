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
    weight_t* result[2];
    result[size % 2] = distance;
    result[(size + 1) % 2] = new weight_t[size * size];

    for(dataSize_t i = 0; i < size ; ++i){
        for(dataSize_t j = 0; j < size; ++j){
            for(dataSize_t k = 0; k < 2; ++k){
                result[k][i * size + j] = kWeightInf;
            }
        }
    }

    for(dataSize_t i = 0; i < size; ++i){
        for(dataSize_t k = 0; k < 2; ++k){
            result[k][i * size + i] = 0;
        }
    }

    for(dataSize_t i = 0; i < size; ++i){
        for(dataSize_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
            for(dataSize_t k = 0; k < 2; ++k){
                result[k][i * size + colIdx[j]] = value[j];
            }
        }
    }
    
    // int nested = omp_get_nested();
    // omp_set_nested(1);
    for(dataSize_t k = 0; k < size; ++k){
        int now = k % 2;
        int next = (k + 1) % 2;
        #pragma omp parallel
        {
            #pragma omp for
            for(dataSize_t i = 0; i < size; ++i){
                if((result[now][i * size + k] == kWeightInf))
                    continue;
                // #pragma omp parallel
                // {
                //     #pragma omp for
                        for(dataSize_t j = 0; j < size; ++j){
                            if(result[now][k * size + j] == kWeightInf)
                                continue;
                            weight_t oldValue = result[now][i * size + j];
                            weight_t newValue = result[now][i * size + k] + result[now][k * size + j];
                            // result[next][i * size + j] = (oldValue > newValue) ? newValue : oldValue;
                            if(newValue > oldValue)
                                newValue = oldValue;
                            oldValue = result[next][i * size + j];
                            if(oldValue > newValue)
                                result[next][i * size + j] = newValue;
                        }
                // }
            }
        }
    }
    // omp_set_nested(nested);

    delete[] result[(size + 1) % 2];
}
