// Elimination tree library for supernodal floyd warshall algorithm //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#ifndef INCLUDE_SUPERNODE_HPP
#define INCLUDE_SUPERNODE_HPP

#include "Common/Type.hpp"
#include "Graph/CSRGraph.hpp"
#include <vector>

class Supernode{
public:
    dataSize_t from;
    dataSize_t to;
    dataSize_t left;
    dataSize_t right;
    dataSize_t parent;
    Supernode(dataSize_t from, dataSize_t to, dataSize_t left, dataSize_t right, dataSize_t parent);
};

class ETree{
private:
    dataSize_t* perm;
    dataSize_t* iperm;
    dataSize_t* height;
    dataSize_t* maxBy;
    CSRGraph* permGraph;
    std::vector<Supernode> supernodes;
    std::vector<dataSize_t> depthIndex; //index of each depth
public:
    const dataSize_t* const Perm() const;
    const dataSize_t* const IPerm() const;
    const dataSize_t* const Height() const;
    const dataSize_t* const MaxBy() const;
    const std::vector<Supernode>& const Supernodes() const;
    const std::vector<dataSize_t>& const DepthIndex() const;
    const CSRGraph& const PermGraph() const;
    ETree(const CSRGraph& input_graph);
    ~ETree();
};

#endif