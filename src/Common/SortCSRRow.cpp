#include "Common/SortCSRRow.hpp"
#include "Common/Type.hpp"

namespace{
    template<typename T>
    void swap(T& a, T&b){
        T temp = a;
        a = b;
        b = temp;
    }
}

//Quick sort + insert sort(In small cases)
void SortCSRRow(dataSize_t* key, weight_t* data, dataSize_t length){
    constexpr dataSize_t kBubbleSortThreshHold = 16;
    if(length < kBubbleSortThreshHold){
        for(dataSize_t i = 0; i < length; ++i){
            dataSize_t minimum = i;
            for(dataSize_t j = i + 1; j < length; ++j){
                if(key[minimum] > key[j])
                    minimum = j;
            }
            if(i == minimum) continue;
            swap(key[i], key[minimum]);
            swap(data[i], data[minimum]);
        }
        return;
    }
    dataSize_t pivot = key[0];
    dataSize_t LIdx = 1; //LeftIndex
    dataSize_t RIdx = length - 1; //RightIndex
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
    SortCSRRow(key, data, LIdx);
    SortCSRRow(key + LIdx, data + LIdx, length - LIdx);
}