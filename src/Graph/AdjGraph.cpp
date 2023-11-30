// Adjacent graph library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include "Graph/AdjGraph.hpp"
#include "Common/Type.hpp"

size_t AdjGraph::Size() const{
    return size;
}

AdjGraph::AdjGraph(size_t size){
    size_t matSize = size * size;
    this->size = size;
    mat = new weight_t[matSize];
    for(int i = 0; i < matSize; ++i){
        mat[i] = kWeightInf;
    }
    //Same with
    //for(int i = 0; i < size; i){
    //mat[i * size + i] = 0;
    //}
    for(int i = 0; i < matSize; i += (size + 1)){
        mat[i] = 0;
    }
}

AdjGraph::AdjGraph(size_t size, weight_t* src){
    size_t matSize = size * size;
    this->size = size;
    mat = new weight_t[matSize];
    for(int i = 0; i < matSize; ++i){
        mat[i] = src[i];
    }
}

AdjGraph::AdjGraph(size_t size, weight_t** src){
    size_t matSize = size * size;
    this->size = size;
    mat = new weight_t[matSize];
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            mat[i * size + j] = src[i][j];
        }
    }
}

AdjGraph::AdjGraph(const AdjGraph& src, size_t extend){
    size_t srcSize = src.Size();
    size_t size = srcSize + extend;
    size_t matSize = size * size;
    this->size = size;
    mat = new weight_t[matSize];
    for(int i = 0; i < srcSize; ++i){
        for(int j = 0; j < srcSize; ++j){
            mat[i * size + j] = src[i * srcSize + j];
        }
        for(int j = srcSize; j < srcSize + extend; ++j){
            mat[i * size + j] = kWeightInf;
        }
    }
    //Same with 
    // for(int i = srcSize; i < srcSize + extend; ++i){
    //     for(int j = 0; j < srcSize + extend; ++j){
    //         mat[idxMap(i,j)] = kWeightInf;
    //     }
    //     mat[idxMap(i,i)] = 0;
    // }
    for(int i = srcSize * size; i < matSize; ++i){
        mat[i] = kWeightInf;
    }
    for(int i = srcSize * size + srcSize; i < matSize; i += size + 1){
        mat[i] = 0;
    }
}

AdjGraph::~AdjGraph(){
    delete[] mat;
}


weight_t& AdjGraph::operator[](size_t idx) const{
    return mat[idx];
}