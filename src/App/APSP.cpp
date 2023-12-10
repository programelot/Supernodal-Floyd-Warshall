// All Pair Shortest Path Algorithm testing Application //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include <fstream>
#include <time.h>
#include <string>

#include "Graph/Graph.hpp"
#include "Graph/Converter.hpp"
#include "Graph/CSRGraph.hpp"
#include "Matrix/MtxReader.hpp"
#include "Algorithm/APSP.hpp"
#include "Common/StringFunc.hpp"

int main(int argc, char* argv[]){
    if(argc != 3){
        printf("Usage : APSP MatrixFile Output\n");
        printf("Output will be ignored if the Output name is ignore\n");
        return 0;
    }
    clock_t start, finish;

    Graph g = MtxReader::Instance().Read(argv[1]);

    CSRGraph csr = Converter::Instance().ToCSR(g);
    weight_t* result = nullptr;
    start = clock();
    APSP(csr, &result);
    finish = clock();

    printf("%fsec\n", (double)(finish - start)/CLOCKS_PER_SEC);

    if(!StringFunc::Instance().StrCmp(argv[2], "ignore")){
        std::ofstream reportFile;
        reportFile.open(argv[2]);
        size_t N = g.Size();
        for(int i = 0; i < N ; ++i){
            for(int j = 0; j < N; ++j){
                std::string value = std::to_string(result[i * N + j]);
                reportFile.write(value.c_str(),value.size());
                reportFile.write(" ", 1);
            }reportFile.write("\n", 1);
        }
        reportFile.close();
    }
    delete[] result;
    return 0;
}