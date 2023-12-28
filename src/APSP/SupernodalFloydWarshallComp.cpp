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

    void UpdateBlock(weight_t* distance, dataSize_t size, dataSize_t bI, dataSize_t bJ, dataSize_t bK, const std::vector<Supernode>& supernodes,
                    dataSize_t* inEdge, dataSize_t* inEdgeNum,
                    dataSize_t* outEdge, dataSize_t* outEdgeNum,
                    dataSize_t numSupernodes){
        dataSize_t fromK = supernodes[bK].from;
        dataSize_t toK = supernodes[bK].to;
        toK = toK < size ? toK : size;
        for(dataSize_t k = fromK; k < toK; ++ k){
            dataSize_t inEdgeNumOrg = inEdgeNum[k * numSupernodes + bI];
            dataSize_t outEdgeNumOrg = outEdgeNum[k * numSupernodes + bJ];
            dataSize_t fromBase = k * size + supernodes[bI].from;
            dataSize_t toBase = k * size + supernodes[bJ].from;
            for(dataSize_t i = 0; i < inEdgeNumOrg; ++i){
                for(dataSize_t j = 0; j < outEdgeNumOrg; ++j){
                    dataSize_t from = inEdge[fromBase + i];
                    dataSize_t to = outEdge[toBase + j];
                    weight_t oldValue = distance[from * size + to];
                    weight_t newValue = distance[from * size + k] + distance[k * size + to];
                    if(oldValue > newValue){
                        if(oldValue == kWeightInf){
                            outEdge[from * size + supernodes[bJ].from + outEdgeNum[from * numSupernodes + bJ]] = to;
                            ++outEdgeNum[from * numSupernodes + bJ];
                            inEdge[to * size + supernodes[bI].from  + inEdgeNum[to * numSupernodes + bI]] = from;
                            ++inEdgeNum[to * numSupernodes + bI];
                        }
                        distance[from * size + to] = newValue;
                    }
                }   
            }
        }
    }

    void UpdateBlock(weight_t* distance, dataSize_t size, dataSize_t bI, dataSize_t bJ, dataSize_t bK, const std::vector<Supernode>& supernodes, omp_lock_t* lock,
                    dataSize_t* inEdge, dataSize_t* inEdgeNum,
                    dataSize_t* outEdge, dataSize_t* outEdgeNum,
                    dataSize_t numSupernodes){
        omp_set_lock(lock);
        UpdateBlock(distance, size, bI, bJ, bK, supernodes, inEdge, inEdgeNum, outEdge, outEdgeNum, numSupernodes);
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

    int nested = omp_get_nested();
    omp_set_nested(1);

    const std::vector<Supernode>& supernodes = eTree.Supernodes();
    const std::vector<dataSize_t>& depthIndex = eTree.DepthIndex();

    dataSize_t numSupernodes = supernodes.size();

    omp_lock_t* locks = new omp_lock_t[numSupernodes * numSupernodes];

    for(dataSize_t i = 0; i < numSupernodes * numSupernodes; ++i){
        omp_init_lock(&locks[i]);
    }

    dataSize_t* inEdge = new dataSize_t[size * size];
    dataSize_t* inEdgeNum = new dataSize_t[size * numSupernodes];
    dataSize_t* outEdge = new dataSize_t[size * size];
    dataSize_t* outEdgeNum = new dataSize_t[size * numSupernodes];
    
    for(dataSize_t i = 0; i < size; ++i){
        for(dataSize_t j = 0; j < numSupernodes; ++j){
            inEdgeNum[i * numSupernodes + j] = 0;
            outEdgeNum[i * numSupernodes + j] = 0;
        }
    }

    {
        //Invermap from vertex index to supernode index
        dataSize_t* nodeMap = new dataSize_t[size];
        for(dataSize_t i = 0; i < numSupernodes; ++i){
            for(dataSize_t j = supernodes[i].from; j < supernodes[i].to; ++j){
                nodeMap[j] = i;
            }
        }

        dataSize_t inNodeIdx = 0;
        dataSize_t outNodeIdx = 0;
        for(dataSize_t i = 0; i < size; ++i){
            inNodeIdx = nodeMap[i];
            for(dataSize_t j = rowPtrPerm[i]; j < rowPtrPerm[i + 1]; ++j){
                dataSize_t to = colIdxPerm[j];
                permDistance[i * size + to] = valuePerm[j];

                outNodeIdx = nodeMap[to];
                
                inEdge[to * size + supernodes[inNodeIdx].from + inEdgeNum[to * numSupernodes + inNodeIdx]]= i;
                ++inEdgeNum[to * numSupernodes + inNodeIdx];
                outEdge[i * size + supernodes[outNodeIdx].from + outEdgeNum[i * numSupernodes + outNodeIdx]]= to;
                ++outEdgeNum[i * numSupernodes + outNodeIdx];
            }
        }
        delete[] nodeMap;
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
        // printf("=================================================\n");
        // printf("Depth : %d\n", depth);
        // printf("Diagonal update : %d\n", depthIndex[depth + 1] - depthIndex[depth]);
        // for(dataSize_t i = depthIndex[depth]; i < depthIndex[depth + 1]; ++i){
        //     dataSize_t numDec = descendants[i].size();
        //     dataSize_t numAns = ancestors[i].size();
        //     printf("Pane/Minplus update : %d\n", numDec + numAns);
        // }
        //Diagonal Update
        {
            #pragma omp parallel
            {
                #pragma omp for
                for(dataSize_t i = depthIndex[depth]; i < depthIndex[depth + 1]; ++i){
                    UpdateBlock(permDistance, size, i, i, i, supernodes,
                        inEdge, inEdgeNum, outEdge, outEdgeNum, numSupernodes);
                }
            }
        }

        //Panel Update
        #pragma omp parallel
        {
            #pragma omp for
            for(dataSize_t i = depthIndex[depth]; i < depthIndex[depth + 1]; ++i){
                dataSize_t numDec = descendants[i].size();
                dataSize_t numAns = ancestors[i].size();
                #pragma omp parallel
                {
                    #pragma omp for
                    for(dataSize_t j = 0; j < numDec + numAns; ++j){
                        if(j < numDec){    
                            UpdateBlock(permDistance, size, descendants[i][j], i, i, supernodes,
                                inEdge, inEdgeNum, outEdge, outEdgeNum, numSupernodes);
                            UpdateBlock(permDistance, size, i, descendants[i][j], i, supernodes,
                                inEdge, inEdgeNum, outEdge, outEdgeNum, numSupernodes);
                        }
                        else{
                            UpdateBlock(permDistance, size, ancestors[i][j - numDec], i, i, supernodes,
                                inEdge, inEdgeNum, outEdge, outEdgeNum, numSupernodes);
                            UpdateBlock(permDistance, size, i, ancestors[i][j - numDec], i, supernodes,
                                inEdge, inEdgeNum, outEdge, outEdgeNum, numSupernodes);
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
                dataSize_t numAns = ancestors[i].size();
                #pragma omp parallel
                {
                    #pragma omp for
                    for(dataSize_t j = 0; j < numAns + numDec ; ++j){
                        if(j < numDec){    
                            for(dataSize_t k = 0; k < numDec + numAns; ++k){
                                if(k < numDec){    
                                    UpdateBlock(permDistance, size, descendants[i][j], descendants[i][k], i, supernodes,
                                        inEdge, inEdgeNum, outEdge, outEdgeNum, numSupernodes);
                                }
                                else{
                                    UpdateBlock(permDistance, size, descendants[i][j], ancestors[i][k - numDec], i, supernodes,
                                        inEdge, inEdgeNum, outEdge, outEdgeNum, numSupernodes);
                                }
                            }
                        }
                        else{
                            for(dataSize_t k = 0; k < numDec + numAns; ++k){
                                if(k < numDec){    
                                    UpdateBlock(permDistance, size, ancestors[i][j - numDec], descendants[i][k], i, supernodes,
                                        inEdge, inEdgeNum, outEdge, outEdgeNum, numSupernodes);
                                }
                                else{
                                    UpdateBlock(permDistance, size, ancestors[i][j - numDec], ancestors[i][k - numDec], i, supernodes, 
                                        &locks[ancestors[i][j - numDec] * numSupernodes + ancestors[i][k - numDec]],
                                        inEdge, inEdgeNum, outEdge, outEdgeNum, numSupernodes);
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

    delete[] inEdge;
    delete[] inEdgeNum;
    delete[] outEdge;
    delete[] outEdgeNum;
}