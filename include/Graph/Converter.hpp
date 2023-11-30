// Graph converter library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //
#include "Common/Type.hpp"
#include "Graph/Graph.hpp"
#include "Graph/CSRGraph.hpp"

#ifndef INCLUDE_Converter_HPP
#define INCLUDE_Converter_HPP
class Converter{
private:
    Converter(){}//Singleton
    static Converter instance;
public:
    static Converter& Instance();
    CSRGraph ToCSR(const Graph& g);
};
#endif