// Johnson's algorithm library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

// Parallelism : O(|V|) //
// Work complexity : O(|E||V|log|V|) //
// Time complexity with PRAM : O(|E||V|) //

#include "Common/Type.hpp"
#include "Graph/CSRGraph.hpp"
#include "Algorithm/APSP.hpp"
#include "Heap/BinaryHeap.hpp"

namespace{
    weight_t min(weight_t a, weight_t b){
        return a > b? b : a;
    }
}

void APSP(const CSRGraph& input_graph, weight_t* distance){

    dataSize_t size = input_graph.Size();
    
    dataSize_t* rowPtr = input_graph.RowPtr();
    dataSize_t* colIdx = input_graph.ColIdx();
    weight_t* value = input_graph.Value();
    
    //Make argumented graph for Johnson's algorithm
    //Graph that has one more vertex than input graph.
    //This vertex attached to other vertex without cost(in other word, weight is zero)

    dataSize_t nnz = input_graph.NEdge();
    dataSize_t size_arg = size + 1;
    dataSize_t nnz_arg = nnz + size * 2;
    dataSize_t* rowPtr_Arg = new dataSize_t[size_arg + 1];
    dataSize_t* colIdx_Arg = new dataSize_t[nnz_arg];
    weight_t* value_Arg = new weight_t[nnz_arg];
    
    for(dataSize_t i = 0; i < size; ++i){
        rowPtr_Arg[i] = rowPtr[i] + i;
        for(dataSize_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
            colIdx_Arg[j + i] = colIdx[j];
            value_Arg[j + i] = value[j];
        }
        colIdx_Arg[rowPtr[i + 1] + i] = size_arg - 1;
        value_Arg[rowPtr[i + 1] + i] = 0;
    }
    rowPtr_Arg[size] = nnz + size;
    rowPtr_Arg[size_arg] = nnz_arg;
    dataSize_t base = rowPtr_Arg[size];
    for(dataSize_t j = 0; j < size; ++j){
        colIdx_Arg[base + j] = j;
        value_Arg[base + j] = 0;
    }

    //Bellman ford algorithm from new vertex
    dataSize_t src_arg = size;
    weight_t* distance_Bellman = new weight_t[size_arg];
    //Initialize
    for(int i = 0; i < size_arg; ++i){
        distance_Bellman[i] = kWeightInf;
    }
    distance_Bellman[src_arg] = 0;
    bool quickExit = false;
    for(dataSize_t repeat = 0; repeat < src_arg - 1; ++repeat){
        if(quickExit) break;
        quickExit = true;
        for(dataSize_t i = 0; i < size_arg; ++i){
             //If i is not reachable, pass
            if(distance_Bellman[i] == kWeightInf){ 
                continue;
            }
            for(dataSize_t j = rowPtr_Arg[i]; j < rowPtr_Arg[i + 1]; ++j){
                weight_t oldVal = distance_Bellman[colIdx_Arg[j]];
                weight_t newVal = distance_Bellman[i] + value_Arg[j];
                if(newVal < oldVal){
                    distance_Bellman[colIdx_Arg[j]] = newVal;
                    quickExit = false;
                }
            }
        }
    }

    //Remove negative edge
    weight_t* value_positive = new weight_t[nnz];
    for(dataSize_t i = 0; i < size; ++i){
        for(dataSize_t j = rowPtr[i]; j < rowPtr[i + 1]; ++j){
            value_positive[j] = distance_Bellman[i] + value[j] - distance_Bellman[colIdx[j]];
        }
    }

    //Dijkstra algorithm with positive edges
    #pragma omp parallel
    {
        #pragma omp for
        for(dataSize_t src = 0; src < size; ++src)
        {
            bool* visited = new bool[size];
            BinaryHeap heap(size);

            //Initialize
            for(int i = 0; i < size; ++i){
                distance[src * size + i] = kWeightInf;
                visited[i] = false;
            }
            distance[src * size + src] = 0;
            visited[src] = true;
            for(dataSize_t i = rowPtr[src]; i < rowPtr[src + 1]; ++i){
                if(colIdx[i] == src) continue;
                heap.Insert(colIdx[i], value_positive[i]);
            }
            while(!heap.IsEmpty()){
                HeapNode top = heap.Pop();
                dataSize_t selected = top.index;
                visited[selected] = true;
                distance[src * size + selected] = top.value;
                for(dataSize_t i = rowPtr[selected]; i < rowPtr[selected + 1]; ++i){
                    if(visited[colIdx[i]]) continue;
                    if(heap.Inserted(colIdx[i])){
                        weight_t newValue = value_positive[i] + distance[src * size  + selected];
                        weight_t oldData = heap.Get(colIdx[i]);
                        if(newValue < oldData){
                            heap.Update(colIdx[i], newValue);
                        }
                    }
                    else{
                        heap.Insert(colIdx[i], value_positive[i] + distance[src * size  + selected]);
                    }
                }
            }
            delete[] visited;
        }
    }

    //Get correct distance from algorithm
    
    for(dataSize_t i = 0; i < size; ++i){   
        for(dataSize_t j = 0; j < size; ++j){
            distance[i * size + j] += distance_Bellman[j] - distance_Bellman[i];         
        }
    }

    delete[] distance_Bellman;
    delete[] rowPtr_Arg;
    delete[] colIdx_Arg;
    delete[] value_Arg;
    delete[] value_positive;
}