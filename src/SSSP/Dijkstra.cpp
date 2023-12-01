// Djikstra algorithm library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include "Graph/CSRGraph.hpp"
#include "Algorithm/SSSP.hpp"
#include "Common/Type.hpp"
#include "Heap/BinaryHeap.hpp"

void SSSP(int src, const CSRGraph& input_graph, weight_t** distance){

    size_t size = input_graph.Size();
    size_t* rowPtr = input_graph.RowPtr();
    size_t* colIdx = input_graph.ColIdx();
    weight_t* value = input_graph.Value();


    *distance = new weight_t[size];
    weight_t* result = *distance;
    BinaryHeapTicket** ticket = new BinaryHeapTicket*[size];
    bool* visited = new bool[size];
    BinaryHeap heap;

    //Initialize
    for(int i = 0; i < size; ++i){
        result[i] = kWeightInf;
        visited[i] = false;
        ticket[i] = nullptr;
    }
    result[src] = 0;
    visited[src] = true;
    for(size_t i = rowPtr[src]; i < rowPtr[src + 1]; ++i){
        HeapNode data;
        data.index = colIdx[i];
        data.value = value[i];
        ticket[colIdx[i]] = heap.Insert(data);
    }
    while(!heap.isEmpty()){
        HeapNode top = heap.Pop();
        size_t selected = top.index;
        visited[selected] = true;
        result[selected] = top.value;
        for(size_t i = rowPtr[selected]; i < rowPtr[selected + 1]; ++i){
            if(visited[colIdx[i]]) continue;
            if(ticket[colIdx[i]] == nullptr){
                HeapNode data;
                data.index = colIdx[i];
                data.value = value[i] + result[selected];
                ticket[colIdx[i]] = heap.Insert(data);
            }
            else{
                HeapNode data;
                data.index = colIdx[i];
                data.value = value[i] + result[selected];
                HeapNode oldData = heap.Get(ticket[colIdx[i]]);
                if(data.value < oldData.value){
                    heap.Update(ticket[colIdx[i]], data);
                }
            }
        }
    }

    delete[] ticket;
    delete[] visited;  
}