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

void SSSP(dataSize_t src, const CSRGraph& input_graph, weight_t* distance){

    dataSize_t size = input_graph.Size();
    dataSize_t* rowPtr = input_graph.RowPtr();
    dataSize_t* colIdx = input_graph.ColIdx();
    weight_t* value = input_graph.Value();

    //Initialize
    for(int i = 0; i < size; ++i){
        distance[i] = kWeightInf;
    }
    distance[src] = 0;
    bool quickExit = false;
    for(dataSize_t repeat = 0; repeat < size - 1; ++repeat){
        if(quickExit) break;//If it is stable, exit it fast
        quickExit = true;
        for(dataSize_t i = 0; i < size; ++i){
             //If i is not reachable, pass
            if(distance[i] == kWeightInf){ 
                continue;
            }
            for(dataSize_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
                weight_t oldVal = distance[colIdx[j]];
                weight_t newVal = distance[i] + value[j];
                if(newVal < oldVal){
                    distance[colIdx[j]] = newVal;
                    quickExit = false;
                }
            }
        }
    }
}