// Graph converter library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include "Common/Type.hpp"
#include "Graph/Converter.hpp"
#include "Graph/Graph.hpp"
#include "Graph/CSRGraph.hpp"
#include <vector>
#include "Common/DebugAssert.hpp"

namespace{
    template<typename T>
    void swap(T& a, T&b){
        T temp = a;
        a = b;
        b = temp;
    }

    //Quick sort + insert sort(In small cases)
    void sort(size_t* key, weight_t* data, size_t length){
        constexpr size_t kBubbleSortThreshHold = 16;
        if(length < kBubbleSortThreshHold){
            for(size_t i = 0; i < length; ++i){
                size_t minimum = i;
                for(size_t j = i + 1; j < length; ++j){
                    if(key[minimum] > key[j])
                        minimum = j;
                }
                if(i == minimum) continue;
                swap(key[i], key[minimum]);
                swap(data[i], data[minimum]);
            }
            return;
        }
        size_t pivot = key[0];
        size_t LIdx = 1; //LeftIndex
        size_t RIdx = length - 1; //RightIndex
        while(true){
            while(true){//Seek LeftReader that is bigger than pivot
                if(key[LIdx] > pivot) break;
                if(LIdx == RIdx) break;
                ++LIdx;
            }
            while(true){//Seek RightReader that is smaller than pivot
                if(key[RIdx] < pivot) break;
                if(LIdx == RIdx) break;
                --RIdx;
            }
            if(LIdx == RIdx) break; //CrossOverFinished
            //Swap LIndex and RIndex
            swap(key[LIdx], key[RIdx]);
            swap(data[LIdx], data[RIdx]);
        }
        if(pivot > key[LIdx]){
            swap(key[0], key[LIdx]);
            swap(data[0], data[LIdx]);
        }
        sort(key, data, LIdx);
        sort(key + LIdx, data + LIdx, length - LIdx);
    }
}

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
        if(edges.size() < 2) continue;
        //Sort colIdx, Value
        sort(&colIdx[base], &value[base], edges.size());
        for(size_t j = 0; j < edges.size() - 1; ++j){
            AssertDebug(__FILE__, __LINE__,"No multigraph supports", colIdx[base + j] != colIdx[base + j + 1]);
        }
    }
    CSRGraph csr(size, rowPtr, colIdx, value);
    delete[] rowPtr;
    delete[] colIdx;
    delete[] value;
    return csr;
}
