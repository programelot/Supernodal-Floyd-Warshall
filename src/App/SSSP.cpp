// Single Source Shortest Path Algorithm testing Application //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include <fstream>
#include <time.h>
#include <string>

#include "Graph/Graph.hpp"
#include "Graph/Converter.hpp"
#include "Graph/CSRGraph.hpp"
#include "Matrix/MtxReader.hpp"
#include "Algorithm/SSSP.hpp"
#include "Common/StringFunc.hpp"

int main(int argc, char* argv[]){
    if(argc != 3){
        printf("Usage : SSSP MatrixFile Output\n");
        printf("Output will be ignored if the Output name is ignore\n");
        return 0;
    }
    clock_t start, finish;

    Graph g = MtxReader::Instance().Read(argv[1]);

    CSRGraph csr = Converter::Instance().ToCSR(g);
    weight_t* result = new weight_t[csr.Size()];
    start = clock();
    SSSP(0, csr, result);
    finish = clock();
    

    printf("%fsec\n", (double)(finish - start)/CLOCKS_PER_SEC);

    if(!StringFunc::Instance().StrCmp(argv[2], "ignore")){
        std::ofstream reportFile;
        reportFile.open(argv[2]);
        dataSize_t N = g.Size();
        for(int i = 0; i < N ; ++i){
            std::string value = std::to_string(result[i]);
            reportFile.write(value.c_str(),value.size());
            reportFile.write(" ", 1);
        }
        reportFile.write("\n", 1);
        reportFile.close();
    }
    delete[] result;
    return 0;
}