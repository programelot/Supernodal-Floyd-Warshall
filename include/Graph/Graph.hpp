// Standard graph library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#ifndef INCLUDE_GRAPH_HPP
#define INCLUDE_GRAPH_HPP

#include "Common/Type.hpp"
#include <vector>

class Vertex;

class Edge{
private:
    Vertex* to;
    weight_t value;
public:
    const Vertex& To() const;
    weight_t Value() const;

    Edge(Vertex* to, weight_t value);
    ~Edge();
};

//A vertex information
//Can connect edge to other vertex
class Vertex{
private:
    std::vector<Edge> edges;
public:
    void Connect(Vertex& v, weight_t value);
    const std::vector<Edge>& Edges() const;
    //Remove unnecessary data
    void Shrink();
};

//Pointer based graph
class Graph{
private:
    Vertex* vertices;
    dataSize_t size;
public:
    Graph(dataSize_t size);
    Graph(Graph&& graph);
    ~Graph();
    void Connect(dataSize_t from, dataSize_t to, weight_t value);
    //Remove unnecessary data
    void Shrink();
    Vertex* Vertices() const;
    dataSize_t Size() const;
};

#endif