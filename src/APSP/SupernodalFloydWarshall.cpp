// FloydWarshall algorithm library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

// Work complexity : O(|V|^3) //
// Work complexity : O() //
// Time complexity with PRAM : O() //

// Aux
// Work : O(|E| + |V|^2 + |V|)
// Time : O(|E| + |V|   + |V| log |V|)

#include "Common/Type.hpp"
#include "Graph/CSRGraph.hpp"
#include "Algorithm/APSP.hpp"
#include "ETree/ETree.hpp"
#include <omp.h>

namespace{

    void UpdateBlock(weight_t* result, dataSize_t size, dataSize_t bI, dataSize_t bJ, dataSize_t bK, const std::vector<Supernode>& supernodes, omp_lock_t* lock){
        if(lock != nullptr){
            omp_set_lock(lock);
        }
        dataSize_t fromI = supernodes[bI].from;
        dataSize_t toI = supernodes[bI].to;
        dataSize_t fromJ = supernodes[bJ].from;
        dataSize_t toJ = supernodes[bJ].to;
        dataSize_t fromK = supernodes[bK].from;
        dataSize_t toK = supernodes[bK].to;
        toI = toI < size ? toI : size;
        toJ = toJ < size ? toJ : size;
        toK = toK < size ? toK : size;
        for(dataSize_t k = fromK; k < toK; ++ k){
            for(dataSize_t i = fromI; i < toI; ++i){
                if((result[i * size + k] == kWeightInf))
                    continue;
                for(dataSize_t j = fromJ; j < toJ; ++j){
                    if(result[k * size + j] == kWeightInf)
                        continue;
                    weight_t oldValue = result[i * size + j];
                    weight_t newValue = result[i * size + k] + result[k * size + j];
                    if(oldValue > newValue)
                        result[i * size + j] = newValue;
                }   
            }
        }
        if(lock != nullptr){
            omp_unset_lock(lock);
        }
    }
}

void APSP(const CSRGraph& input_graph, weight_t* distance){
    
    dataSize_t size = input_graph.Size();
    dataSize_t* rowPtr = input_graph.RowPtr();
    dataSize_t* colIdx = input_graph.ColIdx();
    weight_t* value = input_graph.Value();
    ETree eTree(input_graph);

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

    omp_lock_t* locks = new omp_lock_t[supernodes.size() * supernodes.size()];

    for(dataSize_t i = 0; i < supernodes.size() * supernodes.size(); ++i){
        omp_init_lock(&locks[i]);
    }

    std::vector<std::vector<dataSize_t>> ancestors; // Ancestors
    std::vector<std::vector<dataSize_t>> descendants; // Descendants
    ancestors.resize(supernodes.size());
    descendants.resize(supernodes.size());
    //Find ancestor and descendants
    #pragma omp parallel
    {
        #pragma omp for
        for(dataSize_t i = 0; i < supernodes.size(); ++i){
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
                    UpdateBlock(permDistance, size, i, i, i, supernodes, nullptr);
                }
            }
        }

        //Panel Update
        #pragma omp parallel
        {
            #pragma omp for
            for(dataSize_t i = depthIndex[depth]; i < depthIndex[depth + 1]; ++i){
                #pragma omp parallel
                {
                    #pragma omp for
                    for(dataSize_t j = 0; j < ancestors[i].size() + descendants[i].size(); ++j){
                        dataSize_t numDec = descendants[i].size();
                        if(j < numDec){    
                            UpdateBlock(permDistance, size, descendants[i][j], i, i, supernodes, nullptr);
                            UpdateBlock(permDistance, size, i, descendants[i][j], i, supernodes, nullptr);
                        }
                        else{
                            UpdateBlock(permDistance, size, ancestors[i][j - numDec], i, i, supernodes, nullptr);
                            UpdateBlock(permDistance, size, i, ancestors[i][j - numDec], i, supernodes, nullptr);
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
                dataSize_t numDec = descendants[i].size();
                #pragma omp parallel
                {
                    #pragma omp for
                    for(dataSize_t j = 0; j < ancestors[i].size() + descendants[i].size() ; ++j){
                        if(j < numDec){    
                        // #pragma omp parallel
                        // {
                        //     #pragma omp for
                            for(dataSize_t k = 0; k < ancestors[i].size() + descendants[i].size(); ++k){
                                if(k < numDec){    
                                    UpdateBlock(permDistance, size, descendants[i][j], descendants[i][k], i, supernodes, nullptr);
                                }
                                else{
                                    UpdateBlock(permDistance, size, descendants[i][j], ancestors[i][k - numDec], i, supernodes, nullptr);
                                }
                            }
                        // }
                        }
                        else{
                            for(dataSize_t k = 0; k < ancestors[i].size() + descendants[i].size(); ++k){
                                if(k < numDec){    
                                    UpdateBlock(permDistance, size, ancestors[i][j - numDec], descendants[i][k], i, supernodes, nullptr);
                                }
                                else{
                                    //TODO: Need to control the collision
                                    UpdateBlock(permDistance, size, ancestors[i][j - numDec], ancestors[i][k - numDec], i, supernodes, 
                                        &locks[ancestors[i][j - numDec] * supernodes.size() + ancestors[i][k - numDec]]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    omp_set_nested(nested);
    for(dataSize_t i = 0; i < supernodes.size() * supernodes.size(); ++i){
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