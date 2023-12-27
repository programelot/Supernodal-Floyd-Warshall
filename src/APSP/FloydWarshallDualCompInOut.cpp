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
#include <stdlib.h>

void APSP(const CSRGraph& input_graph, weight_t* distance){
    
    dataSize_t size = input_graph.Size();
    dataSize_t* rowPtr = input_graph.RowPtr();
    dataSize_t* colIdx = input_graph.ColIdx();
    weight_t* value = input_graph.Value();
    weight_t* result[2];
    result[0] = new weight_t[size * size]; //To read
    result[1] = distance; //To write
    
    dataSize_t* updateEdge = new dataSize_t[size * size];
    dataSize_t* updateSize = new dataSize_t[size];

    dataSize_t* inEdge = new dataSize_t[size * size];
    dataSize_t* inEdgeNum = new dataSize_t[size];
    dataSize_t* outEdge = new dataSize_t[size * size];
    dataSize_t* outEdgeNum = new dataSize_t[size];

    omp_lock_t* locks = new omp_lock_t[size];
    
    for(dataSize_t i = 0; i < size ; ++i){
        for(dataSize_t j = 0; j < size; ++j){
            for(dataSize_t k = 0; k < 2; ++k){
                result[k][i * size + j] = kWeightInf;
            }
        }
    }

    for(dataSize_t i = 0; i < size; ++i){
        omp_init_lock(&locks[i]);
        inEdgeNum[i] = 0;
        outEdgeNum[i] = 0;
        updateSize[i] = 0;
        for(dataSize_t k = 0; k < 2; ++k){
            result[k][i * size + i] = 0;
        }
    }

    for(dataSize_t i = 0; i < size; ++i){
        for(dataSize_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
            dataSize_t to = colIdx[j];
            inEdge[to * size + inEdgeNum[to]] = i;
            ++inEdgeNum[to];
            outEdge[i * size + outEdgeNum[i]] = to;
            ++outEdgeNum[i];
            for(dataSize_t k = 0; k < 2; ++k){
                result[k][i * size + to] = value[j];
            }
        }
    }
    
    for(dataSize_t k = 0; k < size; ++k){
        //Update values
        #pragma omp parallel
        {
            #pragma omp for
            for(dataSize_t i = 0; i < size; ++i){
                for(dataSize_t j = 0; j < updateSize[i]; ++j){
                    dataSize_t to = updateEdge[i * size + j];
                    if(result[0][i * size + to] > result[1][i * size + to])
                        result[0][i * size + to] = result[1][i * size + to];
                }
                updateSize[i] = 0;
            }
        }
        dataSize_t inEdgeNumOrg = inEdgeNum[k];
        dataSize_t outEdgeNumOrg = outEdgeNum[k];
        #pragma omp parallel
        {
            #pragma omp for
            for(dataSize_t i = 0; i < inEdgeNumOrg; ++i){
                for(dataSize_t j = 0; j < outEdgeNumOrg; ++j){
                    dataSize_t from = inEdge[k * size + i];
                    dataSize_t to = outEdge[k * size + j];
                    weight_t oldValue = result[0][from * size + to];
                    weight_t newValue = result[0][from * size + k] + result[0][k * size + to];
                    // result[next][i * size + to] = (oldValue > newValue) ? newValue : oldValue;
                    if(newValue > oldValue){
                        newValue = oldValue;
                    }
                    else if((oldValue == kWeightInf) && (newValue < kWeightInf)){
                        //Critical section
                        {
                            omp_set_lock(&locks[to]);
                            inEdge[to * size + inEdgeNum[to]] = from;
                            ++inEdgeNum[to];
                            omp_unset_lock(&locks[to]);
                        }
                        outEdge[from * size + outEdgeNum[from]] = to;
                        ++outEdgeNum[from];
                    }
                    oldValue = result[1][from * size + to];
                    if(oldValue > newValue){
                        updateEdge[from * size + updateSize[from]] = to;
                        ++updateSize[from];
                        result[1][from * size + to] = newValue;
                    }
                }
            }
        }
    }
    for(dataSize_t i = 0; i < size; ++i){
        omp_destroy_lock(&locks[i]);
    }
    delete[] result[0];
    delete[] locks;
    delete[] inEdge;
    delete[] inEdgeNum;
    delete[] outEdge;
    delete[] outEdgeNum;
    delete[] updateEdge;
    delete[] updateSize;
}
