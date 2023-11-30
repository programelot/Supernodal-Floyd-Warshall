// Graph converter library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //
#include "Common/Type.hpp"
#include "Graph/Converter.hpp"
#include "Graph/Graph.hpp"
#include "Graph/CSRGraph.hpp"
#include <vector>

Converter& Converter::Instance(){
    static Converter instance;
    return instance;
}

CSRGraph Converter::ToCSR(const Graph& g){
    size_t size = g.Size();
    Vertex* vertices = g.Vertices();
    size_t* rowPtr = new size_t[size + 1];
    rowPtr[0] = 0;
    //Compute rowPtr size
    for(int i = 0; i < size; ++i){
        rowPtr[i + 1] = rowPtr[i] + vertices[i].Edges().size();
    }
    size_t nnz = rowPtr[size];
    size_t* colIdx = new size_t[nnz];
    weight_t* value = new weight_t[nnz];
    for(int i = 0; i < size; ++i){
        size_t base = rowPtr[i];
        const std::vector<Edge>& edges = vertices[i].Edges();
        for(int j = 0; j < edges.size(); ++j){
            //Vertices exists in sequence
            //Compute by address
            colIdx[base + j] = (&edges[j].To() - &vertices[0]);
            value[base + j] = edges[j].Value();
        }
    }
    CSRGraph csr(size, rowPtr, colIdx, value);
    delete[] rowPtr;
    delete[] colIdx;
    delete[] value;
    return csr;
}
