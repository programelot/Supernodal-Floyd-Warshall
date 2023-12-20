// FloydWarshall algorithm library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

// Parallelism : O(1) //
// Work complexity : O(|V|^3) //
// Time complexity with PRAM : O(|V|^3) //

#include "Common/Type.hpp"
#include "Graph/CSRGraph.hpp"
#include "Algorithm/APSP.hpp"

void APSP(const CSRGraph& input_graph, weight_t* distance){
    
    dataSize_t size = input_graph.Size();
    dataSize_t* rowPtr = input_graph.RowPtr();
    dataSize_t* colIdx = input_graph.ColIdx();
    weight_t* value = input_graph.Value();

    for(dataSize_t i = 0; i < size ; ++i){
        for(dataSize_t j = 0; j < size; ++j){
            distance[i * size + j] = kWeightInf;
        }
    }

    for(dataSize_t i = 0; i < size; ++i){
        distance[i * size + i] = 0;
    }

    for(dataSize_t i = 0; i < size; ++i){
        for(dataSize_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
            distance[i * size + colIdx[j]] = value[j];
        }
    }

    for(dataSize_t k = 0; k < size; ++k){
        for(dataSize_t i = 0; i < size; ++i){
            for(dataSize_t j = 0; j < size; ++j){
                weight_t oldValue = distance[i * size + j];
                weight_t newValue = distance[i * size + k] + distance[k * size + j];
                if(oldValue > newValue)
                    distance[i * size + j] = newValue;
            }   
        }
    }
}