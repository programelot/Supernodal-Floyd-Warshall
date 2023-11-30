// Bllman-ford algorithm library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include "Graph/CSRGraph.hpp"
#include "Algorithm/SSSP.hpp"
#include "Common/Type.hpp"

namespace{
    weight_t min(weight_t a, weight_t b){
        return a > b? b : a;
    }
}

void SSSP(int src, const CSRGraph& input_graph, weight_t** distance){

    size_t size = input_graph.Size();
    size_t* rowPtr = input_graph.RowPtr();
    size_t* colIdx = input_graph.ColIdx();
    weight_t* value = input_graph.Value();

    *distance = new weight_t[size];
    weight_t* result = *distance;
    //Initialize
    for(int i = 0; i < size; ++i){
        result[i] = kWeightInf;
    }
    result[src] = 0;
    for(size_t i = rowPtr[src]; i < rowPtr[src + 1]; ++i){
        result[colIdx[i]] = value[i];
    }
    for(size_t repeat = 0; repeat < size - 1; ++repeat){
        for(size_t i = 0; i < size; ++i){
             //If i is not reachable, pass
            if(result[i] == kWeightInf){ 
                continue;
            }
            for(size_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
                result[colIdx[j]] = min(result[colIdx[j]], result[i] + value[j]);
            }
        }
    }
}