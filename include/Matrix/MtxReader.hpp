#include "Graph/Graph.hpp"

class MtxReader{
    public:
    static MtxReader& Instance();
    Graph Read(const char* path);
    private:
    MtxReader(){}
};