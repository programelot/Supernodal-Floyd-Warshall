// Standard graph library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //
#include "Common/Type.hpp"
#include <vector>

#ifndef INCLUDE_Graph_HPP
#define INCLUDE_Graph_HPP

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
    void connect(Vertex& const v, weight_t value);
    const std::vector<Edge>& Edges() const;
    //Remove unnecessary data
    void shrink();
};

//Pointer based graph
class Graph{
private:
    Vertex* vertices;
    size_t size;
public:
    Graph(size_t size);
    ~Graph();
    void connect(size_t from, size_t to, weight_t value);
    //Remove unnecessary data
    void shrink();
    Vertex* Vertices() const;
    size_t Size() const;
};

#endif