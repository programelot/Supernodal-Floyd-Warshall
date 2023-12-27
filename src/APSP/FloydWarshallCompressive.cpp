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

    //Data structures for fast edge iteration with small number of edges.
    //Store the index of edges that exists.
    //It is not guarantee to be ordered.
    dataSize_t* inEdge = new dataSize_t[size * size];
    dataSize_t* inEdgeNum = new dataSize_t[size];
    dataSize_t* outEdge = new dataSize_t[size * size];
    dataSize_t* outEdgeNum = new dataSize_t[size];

    for(dataSize_t i = 0; i < size ; ++i){
        for(dataSize_t j = 0; j < size; ++j){
            distance[i * size + j] = kWeightInf;
        }
    }

    for(dataSize_t i = 0; i < size; ++i){
        distance[i * size + i] = 0;
        inEdgeNum[i] = 0;
        outEdgeNum[i] = 0;
    }

    for(dataSize_t i = 0; i < size; ++i){
        for(dataSize_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
            dataSize_t to = colIdx[j];
            distance[i * size + to] = value[j];
            outEdge[i * size + outEdgeNum[i]] = to;
            ++outEdgeNum[i];
            inEdge[to * size + inEdgeNum[to]] = i;
            ++inEdgeNum[to];
        }
    }

    for(dataSize_t k = 0; k < size; ++k){
        for(dataSize_t i = 0; i < inEdgeNum[k]; ++i){
            for(dataSize_t j = 0; j < outEdgeNum[k]; ++j){
                dataSize_t from = inEdge[k *size + i];
                dataSize_t to = outEdge[k * size + j];
                weight_t oldValue = distance[from * size + to];
                weight_t newValue = distance[from * size + k] + distance[k * size + to];
                if(oldValue > newValue){
                    if(oldValue == kWeightInf){
                        outEdge[from * size + outEdgeNum[from]] = to;
                        ++outEdgeNum[from];
                        inEdge[to * size + inEdgeNum[to]] = from;
                        ++inEdgeNum[to];
                    }
                    distance[from * size + to] = newValue;
                }
            }   
        }
    }

    delete[] inEdge;
    delete[] inEdgeNum;
    delete[] outEdge;
    delete[] outEdgeNum;
}