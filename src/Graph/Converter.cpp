// Graph converter library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include "Common/Type.hpp"
#include "Graph/Converter.hpp"
#include "Graph/Graph.hpp"
#include "Graph/CSRGraph.hpp"
#include <vector>
#include "Common/DebugAssert.hpp"
#include "Common/SortCSRRow.hpp"

Converter& Converter::Instance(){
    static Converter instance;
    return instance;
}

CSRGraph Converter::ToCSR(const Graph& g){
    dataSize_t size = g.Size();
    Vertex* vertices = g.Vertices();
    dataSize_t* rowPtr = new dataSize_t[size + 1];
    rowPtr[0] = 0;
    //Compute rowPtr size
    for(int i = 0; i < size; ++i){
        rowPtr[i + 1] = rowPtr[i] + vertices[i].Edges().size();
    }
    dataSize_t nnz = rowPtr[size];
    dataSize_t* colIdx = new dataSize_t[nnz];
    weight_t* value = new weight_t[nnz];
    for(int i = 0; i < size; ++i){
        dataSize_t base = rowPtr[i];
        const std::vector<Edge>& edges = vertices[i].Edges();
        for(int j = 0; j < edges.size(); ++j){
            //Vertices exists in sequence
            //Compute by address
            colIdx[base + j] = (&edges[j].To() - &vertices[0]);
            value[base + j] = edges[j].Value();
        }
        if(edges.size() < 2) continue;
        //Sort colIdx, Value
        SortCSRRow(&colIdx[base], &value[base], edges.size());
        for(dataSize_t j = 0; j < edges.size() - 1; ++j){
            DebugAssert(__FILE__, __LINE__,"No multigraph supports", colIdx[base + j] != colIdx[base + j + 1]);
        }
    }
    CSRGraph csr(size, rowPtr, colIdx, value);
    delete[] rowPtr;
    delete[] colIdx;
    delete[] value;
    return csr;
}
