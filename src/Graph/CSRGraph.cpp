// Adjacent graph library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include "Graph/CSRGraph.hpp"
#include "Common/Type.hpp"

//Compressed sparse row Graphrix
//https://en.wikipedia.org/wiki/Sparse_Graphrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_forGraph)


size_t CSRGraph::Size() const{
    return size;
}

CSRGraph::CSRGraph(size_t size, size_t* const rowPtr, size_t* const colIdx, weight_t* const value){
    size_t nnz = rowPtr[size];
    this->size = size;
    this->rowPtr = new size_t[size + 1];
    this->colIdx = new size_t[nnz];
    this->value = new weight_t[nnz];
    for(int i = 0; i <= size; ++i){
        this->rowPtr[i] = rowPtr[i];
    }
    for(int i = 0; i < nnz; ++i){
        this->colIdx[i] = colIdx[i];
        this->value[i] = value[i];
    }
}

CSRGraph::CSRGraph(const CSRGraph& src){
    size_t size = src.Size();
    this->size = size;
    size_t nnz = src.RowPtr()[size];
    rowPtr = new size_t[size];
    colIdx = new size_t[nnz];
    value = new weight_t[nnz];
    size_t* srcRowPtr = src.RowPtr();
    size_t* srcColIdx = src.ColIdx();
    weight_t* srcVale = src.Value();
    for(int i = 0; i < size + 1; ++i){
        rowPtr[i] = srcRowPtr[i];
    }
    for(int i = 0; i < nnz; ++i){
        colIdx[i] = srcColIdx[i];
        value[i] = srcVale[i];
    }
}

CSRGraph::CSRGraph(CSRGraph&& src){
    this->size = src.Size();
    rowPtr = src.RowPtr();
    colIdx = src.ColIdx();
    value = src.Value();
    src.rowPtr = nullptr;
    src.colIdx = nullptr;
    src.value = nullptr;
}

CSRGraph::~CSRGraph(){
    if(rowPtr != nullptr) delete[] rowPtr;
    if(colIdx != nullptr) delete[] colIdx;
    if(value != nullptr) delete[] value;
    rowPtr = nullptr;
    colIdx = nullptr;
    value = nullptr;
}

size_t* CSRGraph::RowPtr() const{
    return rowPtr;
}

size_t* CSRGraph::ColIdx() const{
    return colIdx;

}

weight_t* CSRGraph::Value() const{
    return value;

}
