// Binary heap using std::Vector library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#ifndef INCLUDE_BINARYHEAP_HPP
#define INCLUDE_BINARYHEAP_HPP

#include "Common/Type.hpp"
#include <vector>

typedef struct{
    //Two value to be stored
    dataSize_t index;
    weight_t value;
}HeapNode;

//Heap interface
class BinaryHeap{
private:
    //Heapfied list
    HeapNode* heap;
    dataSize_t* indexMap;
    dataSize_t size;
    dataSize_t occupiedSize;
    //Fix heap by upheap and downheap
    void FixHeap(dataSize_t index);
    //Swap data in heapifed list.
    //It will fix map at the same time.
    void swap(dataSize_t a, dataSize_t b);
public:
    BinaryHeap(dataSize_t size);
    BinaryHeap(BinaryHeap &&binaryHeap);
    ~BinaryHeap();

    void Insert(int index, weight_t value);
    HeapNode Pop();
    void Update(int index, weight_t value);
    bool Inserted(int index);
    weight_t Get(int index);
    dataSize_t Size();
    bool IsEmpty();
    void Clear();
    // void print();
};

#endif