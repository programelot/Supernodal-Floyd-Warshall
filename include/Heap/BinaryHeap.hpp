// Binary heap using std::Vector library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //
#include "Common/Type.hpp"

#ifndef INCLUDE_BINARYHEAP_HPP
#define INCLUDE_BINARYHEAP_HPP

#include <vector>

typedef struct{
    //Two value to be stored
    size_t index;
    weight_t value;
}HeapNode;

class BinaryHeapTicket{
private:
    friend class BinaryHeap;
    //Index that can be used to retrive position from HeapNode
    size_t index;
};

//Heap interface
class BinaryHeap{
private:
    //Heapfied list
    std::vector<HeapNode*> heap;
    std::vector<BinaryHeapTicket*> heapTicket;
    //Map table from insert order to actual index of heap
    //Fix heap by upheap and downheap
    void FixHeap(size_t index);
    //Swap data in heapifed list.
    //It will fix map at the same time.
    void swap(size_t a, size_t b);
public:
    BinaryHeap();
    BinaryHeap(BinaryHeap &&binaryHeap);
    ~BinaryHeap();

    BinaryHeapTicket* Insert(HeapNode data) ;
    HeapNode Pop();
    HeapNode GetMin();
    void Update(BinaryHeapTicket* ticket, HeapNode data);
    HeapNode Get(BinaryHeapTicket* ticket);
    size_t Size();
    bool isEmpty();
    // void print();
};

#endif