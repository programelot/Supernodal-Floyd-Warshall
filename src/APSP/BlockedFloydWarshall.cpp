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

// |B| defined by paralleism
// Parallelism : P = O(|V|^2/|B|^2) //
// Work complexity : O(|V|^3) //
// Time complexity with PRAM : O(|V||B|^2) where |B| is the size of block //

#include "Common/Type.hpp"
#include "Graph/CSRGraph.hpp"
#include "Algorithm/APSP.hpp"
#include <cmath>
#include <thread>

namespace{
    void UpdateBlock(weight_t* distance, dataSize_t size, 
                    dataSize_t bI,dataSize_t bJ, dataSize_t bK,
                    dataSize_t blockSize){
        dataSize_t fromI = bI * blockSize;
        dataSize_t fromJ = bJ * blockSize;
        dataSize_t fromK = bK * blockSize;
        dataSize_t toI = (bI + 1) * blockSize;
        dataSize_t toJ = (bJ + 1) * blockSize;
        dataSize_t toK = (bK + 1) * blockSize;
        toI = toI < size ? toI : size;
        toJ = toJ < size ? toJ : size;
        toK = toK < size ? toK : size;
        for(dataSize_t k = fromK; k < toK; ++ k){
            for(dataSize_t i = fromI; i < toI; ++i){
                if((distance[i * size + k] == kWeightInf))
                    continue;
                for(dataSize_t j = fromJ; j < toJ; ++j){
                    if(distance[k * size + j] == kWeightInf)
                        continue;
                    weight_t oldValue = distance[i * size + j];
                    weight_t newValue = distance[i * size + k] + distance[k * size + j];
                    if(oldValue > newValue)
                        distance[i * size + j] = newValue;
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

    // dataSize_t blockSize = sqrt(sqrt(size));
    // dataSize_t blockSize = sqrt(size);
    // dataSize_t blockSize = size;

    //Block size to use maximum number of thread at Panel update.
    dataSize_t blockSize = 16;
    {
        dataSize_t maxThreadNum = std::thread::hardware_concurrency();
        if(maxThreadNum > 1){
            blockSize = size/maxThreadNum;
        }
        if(blockSize < 16) blockSize = 16;
    }
    dataSize_t numBlocks = (size + blockSize - 1)/blockSize;

    for(dataSize_t i = 0; i < size * size; ++i){
        distance[i] = kWeightInf;
    }

    for(dataSize_t i = 0; i < size; ++i){
        distance[i * size + i] = 0;
    }

    for(dataSize_t i = 0; i < size; ++i){
        for(dataSize_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
            distance[i * size + colIdx[j]] = value[j];
        }
    }
    for(dataSize_t blk = 0; blk < numBlocks; ++blk){
        //Diagonal Update
        UpdateBlock(distance, size, blk, blk, blk, blockSize);
        
        //Panel Update
        #pragma omp parallel
        {
            #pragma omp for
            for(dataSize_t blkE = 0; blkE < numBlocks - 1; ++blkE){
                dataSize_t blkEU = blkE; 
                if(blkEU >= blk) ++blkEU;
                UpdateBlock(distance, size, blkEU, blk, blk, blockSize);
                UpdateBlock(distance, size, blk,  blkEU, blk, blockSize);
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
                    UpdateBlock(distance, size, blkIU, blkJU, blk, blockSize);
                }
            }
        }
    }
}