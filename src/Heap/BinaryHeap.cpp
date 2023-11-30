// Binary heap using std::Vector library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //
#include "Common/Type.hpp"
#include "Heap/BinaryHeap.hpp"
#include <cassert>

void BinaryHeap::FixHeap(size_t index){
    do{//UpHeap
        if(index > 0){//Has parent
            if(heap[index]->value < heap[(index - 1)/2]->value){
                swap(index, (index - 1)/2);
                index = (index - 1)/2;
                continue;
            }
        }
        break;
    }while(true);
    do{//DownHeap
        if(index * 2 + 1 < heap.size()){
            if(index * 2 + 2 < heap.size()){
                size_t compareWith =
                    heap[2 * index + 1]->value < heap[2 * index + 2]->value?
                    2 * index + 1 : 2 * index + 2;
                if(heap[compareWith]->value < heap[index]->value){
                    swap(index, compareWith);
                    index = compareWith;
                    continue;
                }
            }
            else{
                if(heap[index]->value > heap[2 * index + 1]->value){
                    swap(index, 2 * index + 1);
                    index = 2 * index + 1;
                    continue;
                }
            }
        }
        break;
    }while(true);
}

void BinaryHeap::swap(size_t a, size_t b){
    {
        HeapNode* temp = heap[a];
        heap[a] = heap[b];
        heap[b] = temp;
    }
    {
        BinaryHeapTicket* temp = heapTicket[a];
        heapTicket[a] = heapTicket[b];
        heapTicket[b] = temp;
        heapTicket[a]->index = a;
        heapTicket[b]->index = b;
    }
}

size_t BinaryHeap::Size(){
    return heap.size();
}

BinaryHeap::BinaryHeap(){}

BinaryHeap::BinaryHeap(BinaryHeap &&binaryHeap){
    heap = binaryHeap.heap;
    heapTicket = binaryHeap.heapTicket;
    binaryHeap.heap.clear();
    binaryHeap.heapTicket.clear();
}

BinaryHeap::~BinaryHeap(){
    for(int i = 0; i < Size(); ++i){
        delete heap[i];
        delete heapTicket[i];
    }
}

BinaryHeapTicket* BinaryHeap::Insert(HeapNode data){
    HeapNode* newValue = new HeapNode();
    newValue->index = data.index;
    newValue->value = data.value;
    heap.emplace_back(newValue);
    BinaryHeapTicket* binTicket = new BinaryHeapTicket();
    binTicket->index = Size() - 1;
    heapTicket.emplace_back(binTicket);
    FixHeap(Size() - 1);
    return binTicket;
}

HeapNode BinaryHeap::Pop(){
    assert("Can not get min from empty heap" && Size() > 0);
    HeapNode top;
    top.value = heap[0]->value;
    top.index = heap[0]->index;
    swap(0, Size() - 1);
    delete heap[Size() - 1];
    delete heapTicket[Size() - 1];
    heap.pop_back();
    heapTicket.pop_back();
    FixHeap(0);
    return top;
}

HeapNode BinaryHeap::GetMin(){
    assert("Can not get min from empty heap" && Size() > 0);
    HeapNode top;
    top.value = heap[0]->value;
    top.index = heap[0]->index;
    return top;
}

void BinaryHeap::Update(BinaryHeapTicket* ticket, HeapNode data){
    size_t index = ticket->index;
    heap[index]->index = data.index;
    heap[index]->value = data.value;
    FixHeap(index);
}

HeapNode BinaryHeap::Get(BinaryHeapTicket* ticket){
    size_t index = ticket->index;
    return *heap[index];
}

bool BinaryHeap::isEmpty(){
    return heap.empty();
}

// #include <iostream>
// void BinaryHeap::print(){
//     for(int i = 0; i < heap.size(); ++i){
//         printf("(%d %f)", heap[i]->index, heap[i]->value);
//     }printf("\n");
// }