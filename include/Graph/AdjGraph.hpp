// Adjacent graph library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //
#include "Common/Type.hpp"

#ifndef INCLUDE_AdjGraph_HPP
#define INCLUDE_AdjGraph_HPP

//Adjacent Graphrix size of (N by N).
//N can be obtained by Size().
class AdjGraph{
private:
    weight_t* mat;
    size_t size;
public:
    //Return size of the Graphrix.
    size_t Size() const;

    //Basic operations
    AdjGraph(size_t size = 0);
    //Supports 1D/2D array for initialization
    AdjGraph(size_t size, weight_t* src);
    AdjGraph(size_t size, weight_t** src);
    AdjGraph(const AdjGraph& src, size_t extend = 0);
    ~AdjGraph();
    weight_t& operator[](size_t idx) const;
};

#endif