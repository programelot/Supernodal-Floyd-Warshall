// Binary heap using std::Vector library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //
#include "Common/Type.hpp"
#include "Heap/BinaryHeap.hpp"
#include "Common/DebugAssert.hpp"

void BinaryHeap::FixHeap(dataSize_t index){
    do{//UpHeap
        if(index > 0){//Has parent
            if(heap[index].value < heap[(index - 1)/2].value){
                swap(index, (index - 1)/2);
                index = (index - 1)/2;
                continue;
            }
        }
        break;
    }while(true);
    do{//DownHeap
        if(index * 2 + 1 < occupiedSize){
            if(index * 2 + 2 < occupiedSize){
                dataSize_t compareWith =
                    heap[2 * index + 1].value < heap[2 * index + 2].value?
                    2 * index + 1 : 2 * index + 2;
                if(heap[compareWith].value < heap[index].value){
                    swap(index, compareWith);
                    index = compareWith;
                    continue;
                }
            }
            else{
                if(heap[index].value > heap[2 * index + 1].value){
                    swap(index, 2 * index + 1);
                    index = 2 * index + 1;
                    continue;
                }
            }
        }
        break;
    }while(true);
}

void BinaryHeap::swap(dataSize_t a, dataSize_t b){
    if(a == b) return;
    {
        dataSize_t ai = heap[a].index;
        dataSize_t bi = heap[b].index;
        dataSize_t temp = indexMap[ai];
        indexMap[ai] = indexMap[bi];
        indexMap[bi] = temp;
    }
    {
        HeapNode temp = heap[a];
        heap[a] = heap[b];
        heap[b] = temp;
    }
}

dataSize_t BinaryHeap::Size(){
    return size;
}

BinaryHeap::BinaryHeap(dataSize_t size){
    heap = new HeapNode[size];
    indexMap = new dataSize_t[size];
    occupiedSize = 0;
    for(dataSize_t i = 0; i < size; ++i){
        indexMap[i] = size;
    }
    this->size = size;
}

BinaryHeap::BinaryHeap(BinaryHeap &&binaryHeap){
    heap = binaryHeap.heap;
    indexMap = binaryHeap.indexMap;
    occupiedSize = binaryHeap.occupiedSize;
    binaryHeap.heap = nullptr;
    binaryHeap.indexMap = nullptr;
}

BinaryHeap::~BinaryHeap(){
    Clear();
}

void BinaryHeap::Insert(int index, weight_t value){
    DebugAssert(__FILE__, __LINE__, "Can not insert index bigger than size of heap", index < size);
    heap[occupiedSize].index = index;
    heap[occupiedSize].value = value;
    indexMap[index] = occupiedSize;
    FixHeap(occupiedSize);
    ++occupiedSize;
}

HeapNode BinaryHeap::Pop(){
    HeapNode top;
    top.value = heap[0].value;
    top.index = heap[0].index;
    swap(0, occupiedSize - 1);
    indexMap[top.index] = size;
    --occupiedSize;
    FixHeap(0);
    return top;
}

void BinaryHeap::Update(int index, weight_t value){
    dataSize_t mappedIndex = indexMap[index];
    heap[mappedIndex].index = index;
    heap[mappedIndex].value = value;
    FixHeap(mappedIndex);
}

bool BinaryHeap::Inserted(int index){
    if(indexMap[index] == size) return false;
    return true;
}

weight_t BinaryHeap::Get(int index){
    dataSize_t mappedIndex = indexMap[index];
    return heap[mappedIndex].value;
}

bool BinaryHeap::IsEmpty(){
    return occupiedSize == 0;
}

void BinaryHeap::Clear(){
    if(heap != nullptr){
        delete[] heap;
        heap = nullptr;
    }
    if(indexMap != nullptr){
        delete[] indexMap;
        indexMap = nullptr;
    }
}

// #include <iostream>
// void BinaryHeap::print(){
//     for(int i = 0; i < heap.size(); ++i){
//         printf("(%d %f)", heap[i]->index, heap[i]->value);
//     }printf("\n");
// }