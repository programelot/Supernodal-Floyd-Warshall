// All Pair Shortest Path Algorithm testing Application //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include <string>

#include "Graph/Graph.hpp"
#include "Graph/Converter.hpp"
#include "Graph/CSRGraph.hpp"
#include "Matrix/MtxReader.hpp"
#include "Algorithm/APSP.hpp"

int main(int argc, char* argv[]){
    if(argc != 2){
        printf("Usage : NegativeCycleCheck MatrixFile\n");
        return 0;
    }

    Graph g = MtxReader::Instance().Read(argv[1]);
    CSRGraph csr = Converter::Instance().ToCSR(g);

    weight_t* result = nullptr;
    APSP(csr, &result);

    size_t N = g.Size();
    for(int i = 0; i < N ; ++i){
        if(result[i * N + i] < 0){
            printf("Negative cycle dectected\n");
            break;
        }
    }
    delete[] result;
    return 0;
}