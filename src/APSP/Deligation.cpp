// Djikstra algorithm library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include "Common/Type.hpp"
#include "Graph/CSRGraph.hpp"
#include "Algorithm/APSP.hpp"
#include "Algorithm/SSSP.hpp"

void APSP(const CSRGraph& input_graph, weight_t** distance){
    size_t size = input_graph.Size();
    *distance = new weight_t[size * size];
    weight_t* result = *distance;
    for(size_t i = 0; i < size; ++i){
        weight_t* SSSP_result;
        SSSP(i, input_graph, &SSSP_result);
        for(size_t j = 0; j < size; ++j){
            result[size * i + j] = SSSP_result[j];
        }
        delete[] SSSP_result;
    }
}