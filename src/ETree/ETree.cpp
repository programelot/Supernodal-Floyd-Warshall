// Elimination tree library for supernodal floyd warshall algorithm //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include "Common/Type.hpp"
#include "Graph/CSRGraph.hpp"
#include "ETree/Etree.hpp"
#include "Common/DebugAssert.hpp"
#include "Graph/SortCSRRow.hpp"
#include <metis.h>
#include <vector>

Supernode::Supernode(dataSize_t from, dataSize_t to, dataSize_t left, dataSize_t right, dataSize_t parent){
    this->from = from;
    this->to = to;
    this->left = left;
    this->right = right;
    this->parent = parent;
}

ETree::ETree(const CSRGraph& input_graph){

    //Minimum size of clear area to construct tree.
    constexpr dataSize_t clearAreaThreshold = 0;
    constexpr dataSize_t blocksizeThreshold = 16;//32;
    constexpr dataSize_t eTreeDepthMax = 50;

    int size = input_graph.Size();
    dataSize_t nnz = input_graph.NEdge();

    dataSize_t* rowPtr = input_graph.RowPtr();
    dataSize_t* colIdx = input_graph.ColIdx();
    weight_t* value = input_graph.Value();

    perm  = new dataSize_t[size];
    iperm = new dataSize_t[size];
    
    int resultND = METIS_NodeND(&size,rowPtr,colIdx,NULL,NULL,perm,iperm);
    DebugAssert(__FILE__, __LINE__, "Metis nested dissection failed", resultND == METIS_OK);

    //Permutate csr with iperm
    dataSize_t* rowPtrPerm = new dataSize_t[size + 1];
    dataSize_t* colIdxPerm = new dataSize_t[nnz];
    weight_t* valuePerm = new weight_t[nnz];    
    {
        dataSize_t base = 0;
        for(dataSize_t i = 0; i < size; ++i){
            dataSize_t rowOrg = perm[i];
            dataSize_t colSizeOrg = (rowPtr[rowOrg + 1] - rowPtr[rowOrg]);
            dataSize_t baseOrg = rowPtr[rowOrg];
            rowPtrPerm[i] = base;
            rowPtrPerm[i + 1] = base + colSizeOrg;
            for(dataSize_t j = 0; j < colSizeOrg; ++j){
                colIdxPerm[base + j] = iperm[colIdx[baseOrg + j]];
                valuePerm[base + j] = value[baseOrg + j];
            }
            SortCSRRow(&colIdxPerm[base], &valuePerm[base], colSizeOrg);
            base = rowPtrPerm[i + 1];
        }
    }
    permGraph = new CSRGraph(size, rowPtrPerm, colIdxPerm, valuePerm);

    //Count the size of clear area
    // O(|E|)
    dataSize_t leftCol = size;
    height = new dataSize_t[size];
    {
        for(dataSize_t i = 0; i < size; ++i){
            height[i] = 0;
        }
        for(dataSize_t row = 0; row < size; ++row){
            if(leftCol == 0) break;
            for(dataSize_t i = rowPtrPerm[row]; i < rowPtrPerm[row + 1]; ++i){
                if(height[colIdxPerm[i]]  == 0){
                    height[colIdxPerm[i]] = row;
                    --leftCol;
                }
            }
        }
    }
    //Find the maximum size of clear area until certain row
    //Can be used for all cases.
    // O(|V|^2)
    maxBy = new dataSize_t[size];
    #pragma omp parallel
    {
        #pragma omp for
        for(dataSize_t i = 0; i < size; ++i){
            //Can not be a block
            maxBy[i] = i; //Make size of block to zero
            for(dataSize_t j = i; j < size; ++j){
                //Can not be a block with row j.
                if(colIdxPerm[rowPtrPerm[j]] < i)
                    break;
                if(height[j] < i)
                    break;
                maxBy[i] = j + 1;
            }
        }
    }
    
    // O(|V| log |V|)
    // O(|V|)
    supernodes.emplace_back(0,size,size,size,size); //from, to, left, right , parent
    depthIndex.emplace_back(0);
    depthIndex.emplace_back(1);
    dataSize_t depth = 0;
    dataSize_t thisDepthSize;
    do{
        thisDepthSize = 0;
        for(dataSize_t i = depthIndex[depth]; i < depthIndex[depth + 1]; ++i){
            dataSize_t clearAreaSize = 0;
            dataSize_t cutPointA = 0;
            dataSize_t cutPointB = 0;
            //Cut will be [from, cutPointA) [cutPointA, cutPointB) [cutPointB, to)
            for(dataSize_t j = supernodes[i].from + 1; j < supernodes[i].to; ++j){
                dataSize_t cutCanB = maxBy[j]; //Cut candidate for point B
                if(cutCanB == j) continue; //Not a candidate to be cut
                if(cutCanB > supernodes[i].to){
                    cutCanB = supernodes[i].to;
                }
                dataSize_t basedJ = j - supernodes[i].from;
                dataSize_t clearAreaByJ = basedJ * (cutCanB - j);
                if(clearAreaSize < clearAreaByJ){
                    clearAreaSize = clearAreaByJ;
                    cutPointA = j;
                    cutPointB = cutCanB;
                }
            }
            if((clearAreaSize > clearAreaThreshold) &&
                ((cutPointA - supernodes[i].from) > blocksizeThreshold) &&
                ((cutPointB - cutPointA) > blocksizeThreshold) ){
                supernodes.emplace_back(supernodes[i].from, cutPointA, size, size, i);
                supernodes[i].left = thisDepthSize + depthIndex[depth + 1];
                ++thisDepthSize;
                supernodes.emplace_back(cutPointA, cutPointB, size, size, i);
                supernodes[i].right = thisDepthSize + depthIndex[depth + 1];
                supernodes[i].from = cutPointB;
                ++thisDepthSize;
            }
        }
        if(thisDepthSize == 0)
            break;
        depthIndex.emplace_back(depthIndex[depth + 1] + thisDepthSize);
        ++depth;
    }while(thisDepthSize > 0 && depth < eTreeDepthMax);

    delete[] valuePerm;
    delete[] colIdxPerm;
    delete[] rowPtrPerm;
}


ETree::~ETree(){
    delete[] perm;
    delete[] iperm;
    delete[] height;
    delete[] maxBy;
    delete permGraph;
}

const dataSize_t* const ETree::IPerm() const{
    return iperm;
}

const dataSize_t* const ETree::Perm() const{
    return perm;
}

const dataSize_t* const ETree::Height() const{
    return height;
}

const dataSize_t* const ETree::MaxBy() const{
    return maxBy;
}

const CSRGraph& const ETree::PermGraph() const{
    return *permGraph;
}

const std::vector<Supernode>& const ETree::Supernodes() const{
    return supernodes;
}

const std::vector<dataSize_t>& const ETree::DepthIndex() const{
    return depthIndex;
}
