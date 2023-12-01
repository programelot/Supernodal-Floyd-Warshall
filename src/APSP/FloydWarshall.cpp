// FloydWarshall algorithm library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include "Common/Type.hpp"
#include "Graph/CSRGraph.hpp"
#include "Algorithm/APSP.hpp"

void APSP(const CSRGraph& input_graph, weight_t** distance){
    
    size_t size = input_graph.Size();
    size_t* rowPtr = input_graph.RowPtr();
    size_t* colIdx = input_graph.ColIdx();
    weight_t* value = input_graph.Value();
    *distance = new weight_t[size * size];
    weight_t* result = *distance;

    for(size_t i = 0; i < size * size; ++i){
        result[i] = kWeightInf;
    }

    for(size_t i = 0; i < size; ++i){
        result[i * size + i] = 0;
    }

    for(size_t i = 0; i < size; ++i){
        for(size_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
            result[i * size + colIdx[j]] = value[j];
        }
    }

    for(size_t k = 0; k < size; ++k){
       for(size_t i = 0; i < size; ++i){
            for(size_t j = 0; j < size; ++j){
                weight_t oldValue = result[i * size + j];
                weight_t newValue = result[i * size + k] + result[k * size + j];
                if(oldValue > newValue)
                    result[i * size + j] = newValue;
            }   
        }
    }
}