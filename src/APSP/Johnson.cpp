// Johnson's algorithm library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include "Common/Type.hpp"
#include "Graph/CSRGraph.hpp"
#include "Algorithm/APSP.hpp"
#include "Heap/BinaryHeap.hpp"

namespace{
    weight_t min(weight_t a, weight_t b){
        return a > b? b : a;
    }
}

void APSP(const CSRGraph& input_graph, weight_t** distance){

    size_t size = input_graph.Size();
    *distance = new weight_t[size * size];
    weight_t* result = *distance;
    
    size_t* rowPtr = input_graph.RowPtr();
    size_t* colIdx = input_graph.ColIdx();
    weight_t* value = input_graph.Value();
    
    //Make argumented graph for Johnson's algorithm
    //Graph that has one more vertex than input graph.
    //This vertex attached to other vertex without cost(in other word, weight is zero)

    size_t nnz = rowPtr[size];
    size_t size_arg = size + 1;
    size_t nnz_arg = nnz + size * 2;
    size_t* rowPtr_Arg = new size_t[size_arg + 1];
    size_t* colIdx_Arg = new size_t[nnz_arg];
    weight_t* value_Arg = new weight_t[nnz_arg];
    
    for(size_t i = 0; i < size; ++i){
        rowPtr_Arg[i] = rowPtr[i] + i;
        for(size_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
            colIdx_Arg[j + i] = colIdx[j];
            value_Arg[j + i] = value[j];
        }
        colIdx_Arg[rowPtr[i + 1] + i] = size_arg - 1;
        value_Arg[rowPtr[i + 1] + i] = 0;
    }
    rowPtr_Arg[size] = rowPtr[size] + size;
    rowPtr_Arg[size_arg] = nnz_arg;
    size_t base = rowPtr_Arg[size];
    for(size_t j = 0; j < size; ++j){
        colIdx_Arg[base + j] = j + 1;
        value_Arg[base + j] = 0;
    }

    //Bellman ford algorithm from new vertex
    size_t src_arg = size;
    weight_t* result_Bellman = new weight_t[size_arg];
    //Initialize
    for(int i = 0; i < size_arg; ++i){
        result_Bellman[i] = kWeightInf;
    }
    result_Bellman[src_arg] = 0;
    for(size_t i = rowPtr_Arg[src_arg]; i < rowPtr_Arg[src_arg + 1]; ++i){
        result_Bellman[colIdx_Arg[i]] = value_Arg[i];
    }
    for(size_t repeat = 0; repeat < src_arg - 1; ++repeat){
        for(size_t i = 0; i < size_arg; ++i){
             //If i is not reachable, pass
            if(result_Bellman[i] == kWeightInf){ 
                continue;
            }
            for(size_t j = rowPtr_Arg[i]; j < rowPtr_Arg[i + 1]; ++j){
                result_Bellman[colIdx_Arg[j]] = 
                                    min(result_Bellman[colIdx_Arg[j]], 
                                    result_Bellman[i] + value_Arg[j]);
            }
        }
    }

    //Remove negative edge
    weight_t* value_positive = new weight_t[nnz];
    for(size_t i = 0; i < size; ++i){
        for(size_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
            value_positive[j] = result_Bellman[i] + value[j] - result_Bellman[colIdx[j]];
        }
    }

    //Dijkstra algorithm with positive edges

    BinaryHeapTicket** ticket = new BinaryHeapTicket*[size];
    bool* visited = new bool[size];
    BinaryHeap heap;
    
    for(size_t src = 0; src < size; ++src)
    {
        heap.clear();
        //Initialize
        for(int i = 0; i < size; ++i){
            result[src * size + i] = kWeightInf;
            visited[i] = false;
            ticket[i] = nullptr;
        }
        result[src * size + src] = 0;
        visited[src] = true;
        for(size_t i = rowPtr[src]; i < rowPtr[src + 1]; ++i){
            HeapNode data;
            data.index = colIdx[i];
            data.value = value_positive[i];
            ticket[colIdx[i]] = heap.Insert(data);
        }
        while(!heap.isEmpty()){
            HeapNode top = heap.Pop();
            size_t selected = top.index;
            visited[selected] = true;
            result[src * size + selected] = top.value;
            for(size_t i = rowPtr[selected]; i < rowPtr[selected + 1]; ++i){
                if(visited[colIdx[i]]) continue;
                if(ticket[colIdx[i]] == nullptr){
                    HeapNode data;
                    data.index = colIdx[i];
                    data.value = value_positive[i] + result[src * size  + selected];
                    ticket[colIdx[i]] = heap.Insert(data);
                }
                else{
                    HeapNode data;
                    data.index = colIdx[i];
                    data.value = value_positive[i] + result[src * size  + selected];
                    HeapNode oldData = heap.Get(ticket[colIdx[i]]);
                    if(data.value < oldData.value){
                        heap.Update(ticket[colIdx[i]], data);
                    }
                }
            }
        }
    }

    //Get correct distance from algorithm
    
    for(size_t i = 0; i < size; ++i){   
        for(size_t j = 0; j < size; ++j){
            result[i * size + j] += result_Bellman[j] - result_Bellman[i];         
        }
    }

    delete[] ticket;
    delete[] visited;  

    delete[] result_Bellman;
    delete[] rowPtr_Arg;
    delete[] colIdx_Arg;
    delete[] value_Arg;
    delete[] value_positive;
}