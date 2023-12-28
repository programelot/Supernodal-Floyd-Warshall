// FloydWarshall algorithm library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

// Parallelism : O(|V|^2/|B|^2) //
// Work complexity : O(|V|^3) //
// Time complexity with PRAM : O(|V||B|^2) where |B| is the size of block //

// Example 
/// |B| = |V|^(1/4)
// Parallelism : O(|V|^(3/2)) //
// Work complexity : O(|V|^3) //
// Time complexity with PRAM : O(|V|^(3/2)) //

/// |B| = |V|^(1/2)
// Parallelism : O(|V|) //
// Work complexity : O(|V|^3) //
// Time complexity with PRAM : O(|V|^2) //

/// |B| = |V| (Same with Floydwarshall)
// Parallelism : O(1) //
// Work complexity : O(|V|^3) //
// Time complexity with PRAM : O(|V|^3) //

#include "Common/Type.hpp"
#include "Graph/CSRGraph.hpp"
#include "Algorithm/APSP.hpp"
#include <cmath>
#include <thread>

namespace{
    void UpdateBlock(weight_t* distance, dataSize_t size, 
                    dataSize_t bI,dataSize_t bJ, dataSize_t bK,
                    dataSize_t blockSize, dataSize_t numBlocks,
                    dataSize_t* inEdge, dataSize_t* inEdgeNum,
                    dataSize_t* outEdge, dataSize_t* outEdgeNum){
        dataSize_t fromK = bK * blockSize;
        dataSize_t toK = (bK + 1) * blockSize;
        toK = toK < size ? toK : size;
        for(dataSize_t k = fromK; k < toK; ++ k){
            dataSize_t inEdgeNumOrg = inEdgeNum[k * numBlocks + bI];
            dataSize_t outEdgeNumOrg = outEdgeNum[k * numBlocks + bJ];
            dataSize_t fromBase = k * size + bI * blockSize;
            dataSize_t toBase = k * size + bJ * blockSize;
            for(dataSize_t i = 0; i < inEdgeNumOrg; ++i){
                for(dataSize_t j = 0; j < outEdgeNumOrg; ++j){
                    dataSize_t from = inEdge[fromBase + i];
                    dataSize_t to = outEdge[toBase + j];
                    weight_t oldValue = distance[from * size + to];
                    weight_t newValue = distance[from * size + k] + distance[k * size + to];
                    if(oldValue > newValue){
                        if(oldValue == kWeightInf){
                            outEdge[from * size + bJ * blockSize + outEdgeNum[from * numBlocks + bJ]] = to;
                            ++outEdgeNum[from * numBlocks + bJ];
                            inEdge[to * size + bI * blockSize  + inEdgeNum[to * numBlocks + bI]] = from;
                            ++inEdgeNum[to * numBlocks + bI];
                        }
                        distance[from * size + to] = newValue;
                    }
                }   
            }
        }
    }
}

void APSP(const CSRGraph& input_graph, weight_t* distance){
    
    dataSize_t size = input_graph.Size();
    dataSize_t* rowPtr = input_graph.RowPtr();
    dataSize_t* colIdx = input_graph.ColIdx();
    weight_t* value = input_graph.Value();

    //Each block may have different size of workload.
    dataSize_t blockSize = 16;
    {
        //Assign twice of thread for some less workloaded blocks
        dataSize_t maxThreadNum = std::thread::hardware_concurrency() * 2;
        if(maxThreadNum > 1){
            blockSize = size/maxThreadNum;
        }
        if(blockSize < 16) blockSize = 16;
    }
    dataSize_t numBlocks = (size + blockSize - 1)/blockSize;

    dataSize_t* inEdge = new dataSize_t[size * size];
    dataSize_t* inEdgeNum = new dataSize_t[size * numBlocks];
    dataSize_t* outEdge = new dataSize_t[size * size];
    dataSize_t* outEdgeNum = new dataSize_t[size * numBlocks];

    for(dataSize_t i = 0; i < size; ++i){
        for(dataSize_t j = 0; j < numBlocks; ++j){
            inEdgeNum[i * numBlocks + j] = 0;
            outEdgeNum[i * numBlocks + j] = 0;
        }
    }
    
    for(dataSize_t i = 0; i < size * size; ++i){
        distance[i] = kWeightInf;
    }

    for(dataSize_t i = 0; i < size; ++i){
        distance[i * size + i] = 0;
    }

    for(dataSize_t i = 0; i < size; ++i){
        for(dataSize_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
            dataSize_t to = colIdx[j];
            dataSize_t inBlockIdx = i/blockSize;
            dataSize_t outBlockIdx = to/blockSize;

            distance[i * size + to] = value[j];
            
            inEdge[to * size + inBlockIdx * blockSize + inEdgeNum[to * numBlocks + inBlockIdx]]= i;
            ++inEdgeNum[to * numBlocks + inBlockIdx];
            outEdge[i * size + outBlockIdx * blockSize + outEdgeNum[i * numBlocks + outBlockIdx]]= to;
            ++outEdgeNum[i * numBlocks + outBlockIdx];
        }
    }

    for(dataSize_t blk = 0; blk < numBlocks; ++blk){
        //Diagonal Update
        UpdateBlock(distance, size, blk, blk, blk,
                    blockSize, numBlocks,
                    inEdge, inEdgeNum, outEdge, outEdgeNum);
        
        //Panel Update
        #pragma omp parallel
        {
            #pragma omp for
            for(dataSize_t blkE = 0; blkE < numBlocks - 1; ++blkE){
                dataSize_t blkEU = blkE; 
                if(blkEU >= blk) ++blkEU;
                UpdateBlock(distance, size, blkEU, blk, blk,
                                        blockSize, numBlocks,
                                        inEdge, inEdgeNum, outEdge, outEdgeNum);

                UpdateBlock(distance, size, blk, blkEU, blk,
                                        blockSize, numBlocks,
                                        inEdge, inEdgeNum, outEdge, outEdgeNum);
            }
        }

        //MinPlus Outer Product
        
        #pragma omp parallel
        {
            #pragma omp for
            for(dataSize_t blkI = 0; blkI < numBlocks - 1; ++blkI){
                dataSize_t blkIU = blkI; 
                if(blkIU >= blk) ++blkIU;
                for(dataSize_t blkJ = 0; blkJ < numBlocks - 1; ++blkJ){
                    dataSize_t blkJU = blkJ; 
                    if(blkJU >= blk) ++blkJU;
                    UpdateBlock(distance, size, blkIU, blkJU, blk,
                                            blockSize, numBlocks, 
                                            inEdge, inEdgeNum, outEdge, outEdgeNum);
                }
            }
        }
    }

    delete[] inEdge;
    delete[] inEdgeNum;
    delete[] outEdge;
    delete[] outEdgeNum;
}