// Supernodal floydWarshall algorithm library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

// Work complexity : O(|V|^3) //
// Work complexity : O() //
// Time complexity with PRAM : O() //

#include "Common/Type.hpp"
#include "Graph/CSRGraph.hpp"
#include "Algorithm/APSP.hpp"
#include "ETree/ETree.hpp"
#include <omp.h>

namespace{

    void UpdateBlock(weight_t* distance, dataSize_t size, dataSize_t bI, dataSize_t bJ, dataSize_t bK, const std::vector<Supernode>& supernodes){
        dataSize_t fromI = supernodes[bI].from;
        dataSize_t toI = supernodes[bI].to;
        dataSize_t fromJ = supernodes[bJ].from;
        dataSize_t toJ = supernodes[bJ].to;
        dataSize_t fromK = supernodes[bK].from;
        dataSize_t toK = supernodes[bK].to;
        for(dataSize_t k = fromK; k < toK; ++ k){
            for(dataSize_t i = fromI; i < toI; ++i){
                weight_t inValue = distance[i * size + k];
                if(inValue == kWeightInf)
                    continue;
                for(dataSize_t j = fromJ; j < toJ; ++j){
                    weight_t outValue = distance[k * size + j];
                    if(outValue == kWeightInf)
                        continue;
                    weight_t oldValue = distance[i * size + j];
                    weight_t newValue = inValue + outValue;
                    if(oldValue > newValue)
                        distance[i * size + j] = newValue;
                }   
            }
        }
    }

    void UpdateBlock(weight_t* distance, dataSize_t size, dataSize_t bI, dataSize_t bJ, dataSize_t bK, const std::vector<Supernode>& supernodes, omp_lock_t* lock){
        omp_set_lock(lock);
        UpdateBlock(distance, size, bI, bJ, bK, supernodes);
        omp_unset_lock(lock);
    }
}

void APSP(const CSRGraph& input_graph, weight_t* distance){
    
    dataSize_t size = input_graph.Size();
    dataSize_t* rowPtr = input_graph.RowPtr();
    dataSize_t* colIdx = input_graph.ColIdx();
    weight_t* value = input_graph.Value();
    ETree eTree(input_graph, 0, 16, 100);

    weight_t* permDistance = new weight_t[size * size];

    for(dataSize_t i = 0; i < size * size; ++i){
        permDistance[i] = kWeightInf;
    }

    for(dataSize_t i = 0; i < size; ++i){
        permDistance[i * size + i] = 0;
    }

    const CSRGraph& permGraph = eTree.PermGraph();
    dataSize_t* rowPtrPerm = permGraph.RowPtr();
    dataSize_t* colIdxPerm = permGraph.ColIdx();
    weight_t* valuePerm = permGraph.Value();

    for(dataSize_t i = 0; i < size; ++i){
        for(dataSize_t j = rowPtrPerm[i]; j < rowPtrPerm[i + 1]; ++j){
            permDistance[i * size + colIdxPerm[j]] = valuePerm[j];
        }
    }

    int nested = omp_get_nested();
    omp_set_nested(1);

    const std::vector<Supernode>& supernodes = eTree.Supernodes();
    const std::vector<dataSize_t>& depthIndex = eTree.DepthIndex();

    dataSize_t numSupernodes = supernodes.size();

    omp_lock_t* locks = new omp_lock_t[numSupernodes * numSupernodes];

    for(dataSize_t i = 0; i < numSupernodes * numSupernodes; ++i){
        omp_init_lock(&locks[i]);
    }

    std::vector<std::vector<dataSize_t>> ancestors; // Ancestors
    std::vector<std::vector<dataSize_t>> descendants; // Descendants
    ancestors.resize(numSupernodes);
    descendants.resize(numSupernodes);
    //Find ancestor and descendants
    #pragma omp parallel
    {
        #pragma omp for
        for(dataSize_t i = 0; i < numSupernodes; ++i){
            //Search area to be searched
            //Find ancestors
            {
                dataSize_t parent = supernodes[i].parent;
                while(parent != size){
                    ancestors[i].emplace_back(parent);
                    parent = supernodes[parent].parent;
                }
            }
            ///Find descendants
            {
                dataSize_t search = 0;
                dataSize_t left = supernodes[i].left;
                dataSize_t right = supernodes[i].right;
                if(left != size){ // && right != size
                    descendants[i].emplace_back(left);
                    descendants[i].emplace_back(right);
                }
                while(search < descendants[i].size()){
                    left = supernodes[descendants[i][search]].left;
                    right = supernodes[descendants[i][search]].right;
                    if(left != size){// && right != size
                        descendants[i].emplace_back(left);
                        descendants[i].emplace_back(right);
                    }
                    ++search;
                }
            }
        }
    }

    dataSize_t maxDepth = depthIndex.size() - 2;
    for(dataSize_t depth = maxDepth; depth >= 0; --depth){
        //Diagonal Update
        {
            #pragma omp parallel
            {
                #pragma omp for
                for(dataSize_t i = depthIndex[depth]; i < depthIndex[depth + 1]; ++i){
                    UpdateBlock(permDistance, size, i, i, i, supernodes);
                }
            }
        }

        //Panel Update
        #pragma omp parallel
        {
            #pragma omp for
            for(dataSize_t i = depthIndex[depth]; i < depthIndex[depth + 1]; ++i){
                std::vector<dataSize_t>& ances = ancestors[i];
                std::vector<dataSize_t>& desce = descendants[i];
                dataSize_t numDec = desce.size();
                dataSize_t numAns = ances.size();
                #pragma omp parallel
                {
                    #pragma omp for
                    for(dataSize_t j = 0; j < numDec + numAns; ++j){
                        if(j < numDec){    
                            UpdateBlock(permDistance, size, desce[j], i, i, supernodes);
                            UpdateBlock(permDistance, size, i, desce[j], i, supernodes);
                        }
                        else{
                            UpdateBlock(permDistance, size, ances[j - numDec], i, i, supernodes);
                            UpdateBlock(permDistance, size, i, ances[j - numDec], i, supernodes);
                        }
                    }
                }
            }
        }

        //MinPlus Outer Product
        #pragma omp parallel
        {
            #pragma omp for
            for(dataSize_t i = depthIndex[depth]; i < depthIndex[depth + 1]; ++i){
                std::vector<dataSize_t>& ances = ancestors[i];
                std::vector<dataSize_t>& desce = descendants[i];
                dataSize_t numDec = desce.size();
                dataSize_t numAns = ances.size();
                #pragma omp parallel
                {
                    #pragma omp for
                    for(dataSize_t j = 0; j < numAns + numDec ; ++j){
                        if(j < numDec){    
                            for(dataSize_t k = 0; k < numDec + numAns; ++k){
                                if(k < numDec){    
                                    UpdateBlock(permDistance, size, desce[j], desce[k], i, supernodes);
                                }
                                else{
                                    UpdateBlock(permDistance, size, desce[j], ances[k - numDec], i, supernodes);
                                }
                            }
                        }
                        else{
                            for(dataSize_t k = 0; k < numDec + numAns; ++k){
                                if(k < numDec){    
                                    UpdateBlock(permDistance, size, ances[j - numDec], desce[k], i, supernodes);
                                }
                                else{
                                    UpdateBlock(permDistance, size, ances[j - numDec], ances[k - numDec], i, supernodes, 
                                        &locks[ances[j - numDec] * numSupernodes + ances[k - numDec]]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    omp_set_nested(nested);
    for(dataSize_t i = 0; i < numSupernodes * numSupernodes; ++i){
        omp_destroy_lock(&locks[i]);
    }
    delete[] locks;

    //Inverse permutate again
    const dataSize_t* iperm = eTree.IPerm();
    for(dataSize_t i = 0; i < size; ++i){
        for(dataSize_t j = 0; j < size; ++j){
            distance[i * size + j] = permDistance[iperm[i] * size + iperm[j]];
        }
    }
    delete[] permDistance;
}