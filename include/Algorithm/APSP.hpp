// All pair shortest path library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#ifndef INCLUDE_APSP_HPP
#define INCLUDE_APSP_HPP

#include "Common/Type.hpp"

//Standard APSP algorithm
//input_graph is an input graph for APSP algorithm.
//input_graph is used by read-only.
//distance will be used to get the result.
//It will be assigned inside of APSP function.
void APSP(const CSRGraph& input_graph, weight_t* distance);

#endif