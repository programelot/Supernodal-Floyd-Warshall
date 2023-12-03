// Standard graph library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //
#include "Common/Type.hpp"
#include "Graph/Graph.hpp"
#include <vector>

Edge::Edge(Vertex* to, weight_t value){
    this->to = to;
    this->value = value;
}

Edge::~Edge(){}

const Vertex& Edge::To() const{
    return *to;
}

weight_t Edge::Value() const{
    return value;
}

void Vertex::Connect(Vertex& const v, weight_t value){
    //Ignore duplicated edges
    edges.emplace_back(&v, value);
}

const std::vector<Edge>& Vertex::Edges() const{
    return edges;
}

void Vertex::Shrink(){
    edges.shrink_to_fit();
}

Graph::Graph(size_t size){
    this->size = size;
    vertices = new Vertex[size];
}

Graph::Graph(Graph&& graph){
    size = graph.size;
    vertices = graph.vertices;
    graph.vertices = nullptr;
}

Graph::~Graph(){
    if(vertices != nullptr)
        delete[] vertices;
}
void Graph::Connect(size_t from, size_t to, weight_t value){
    vertices[from].Connect(vertices[to], value);
}
void Graph::Shrink(){
    for(int i = 0; i < size; ++i){
        vertices[i].Shrink();
    }
}

Vertex* Graph::Vertices() const{
    return vertices;
}

size_t Graph::Size() const{
    return size;
}