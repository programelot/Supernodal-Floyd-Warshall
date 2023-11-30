// Compressed sparse row graph library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //
#include "Common/Type.hpp"

#ifndef INCLUDE_CSRGraph_HPP
#define INCLUDE_CSRGraph_HPP

//Compressed sparse row Graphrix
//https://en.wikipedia.org/wiki/Sparse_Graphrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_forGraph)

class CSRGraph{
private:
    size_t* rowPtr = nullptr;
    size_t* colIdx = nullptr;
    weight_t* value = nullptr;

    size_t size;
public:
    //Return size of the Graphrix.
    size_t Size() const;
    int ID;

    //Basic operations
    //nnz = number of non-zero
    CSRGraph(size_t size, size_t* const rowPtr, size_t* const colIdx, weight_t* const value);
    CSRGraph(const CSRGraph& src);
    CSRGraph(CSRGraph&& src);
    ~CSRGraph();
    
    size_t* RowPtr() const;
    size_t* ColIdx() const;
    weight_t* Value() const;
};

#endif