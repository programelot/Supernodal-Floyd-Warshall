// Compressed sparse row graph library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#ifndef INCLUDE_CSRGRAPH_HPP
#define INCLUDE_CSRGRAPH_HPP

#include "Common/Type.hpp"

//Compressed sparse row Graphrix
//https://en.wikipedia.org/wiki/Sparse_Graphrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_forGraph)

class CSRGraph{
private:
    dataSize_t* rowPtr = nullptr;
    dataSize_t* colIdx = nullptr;
    weight_t* value = nullptr;

    dataSize_t size;
public:
    //Return size of the Graphrix.
    dataSize_t Size() const;
    dataSize_t NEdge() const;
    int ID;

    //Basic operations
    //nnz = number of non-zero
    CSRGraph(dataSize_t size, dataSize_t* const rowPtr, dataSize_t* const colIdx, weight_t* const value);
    CSRGraph(const CSRGraph& src);
    CSRGraph(CSRGraph&& src);
    ~CSRGraph();
    
    dataSize_t* RowPtr() const;
    dataSize_t* ColIdx() const;
    weight_t* Value() const;
};

#endif