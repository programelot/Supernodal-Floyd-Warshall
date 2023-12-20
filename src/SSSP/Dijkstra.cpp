// Djikstra algorithm library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //
// Time complexity : O(|E|log|V|) //

#include "Graph/CSRGraph.hpp"
#include "Algorithm/SSSP.hpp"
#include "Common/Type.hpp"
#include "Heap/BinaryHeap.hpp"

void SSSP(dataSize_t src, const CSRGraph& input_graph, weight_t* distance){

    dataSize_t size = input_graph.Size();
    dataSize_t* rowPtr = input_graph.RowPtr();
    dataSize_t* colIdx = input_graph.ColIdx();
    weight_t* value = input_graph.Value();

    bool* visited = new bool[size];
    BinaryHeap heap(size);
    //Initialize
    for(dataSize_t i = 0; i < size; ++i){
        distance[i] = kWeightInf;
        visited[i] = false;
    }
    distance[src] = 0;
    visited[src] = true;
    for(dataSize_t i = rowPtr[src]; i < rowPtr[src + 1]; ++i){
        if(colIdx[i] == src) continue;
        heap.Insert(colIdx[i], value[i]);
    }
    while(!heap.IsEmpty()){
        HeapNode top = heap.Pop();
        dataSize_t selected = top.index;
        visited[selected] = true;
        distance[selected] = top.value;
        for(dataSize_t i = rowPtr[selected]; i < rowPtr[selected + 1]; ++i){
            if(visited[colIdx[i]]) continue;
            if(heap.Inserted(colIdx[i])){
                weight_t newValue = value[i] + distance[selected];
                weight_t oldData = heap.Get(colIdx[i]);
                if(newValue < oldData){
                    heap.Update(colIdx[i], newValue);
                }
            }
            else{
                heap.Insert(colIdx[i], value[i] + distance[selected]);
            }
        }
    }

    delete[] visited;  
}