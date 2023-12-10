// Bllman-ford algorithm library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //
// Time complexity : O(|E||V|) //

#include "Graph/CSRGraph.hpp"
#include "Algorithm/SSSP.hpp"
#include "Common/Type.hpp"

namespace{
    weight_t min(weight_t a, weight_t b){
        return a > b? b : a;
    }
}

void SSSP(size_t src, const CSRGraph& input_graph, weight_t** distance){

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
    // bool quickExit = false;
    for(size_t repeat = 0; repeat < size - 1; ++repeat){
        // if(quickExit) break;//If it is stable, exit it fast
        // quickExit = true;
        for(size_t i = 0; i < size; ++i){
             //If i is not reachable, pass
            if(result[i] == kWeightInf){ 
                continue;
            }
            for(size_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
                weight_t oldVal = result[colIdx[j]];
                weight_t newVal = result[i] + value[j];
                if(newVal < oldVal){
                    result[colIdx[j]] = newVal;
                    // quickExit = false;
                }
            }
        }
    }
}