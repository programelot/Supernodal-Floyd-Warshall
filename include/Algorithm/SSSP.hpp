// Single source shortest path library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#ifndef INCLUDE_SSSP_HPP
#define INCLUDE_SSSP_HPP

#include "Common/Type.hpp"

//Standard SSSP algorithm
//input_graph is an input graph for SSSP algorithm.
//input_graph is used by read-only.
//distance will be used to get the result.
//It will be assigned inside of SSSP function.
void SSSP(dataSize_t src, const CSRGraph& input_graph, weight_t* distance);

#endif