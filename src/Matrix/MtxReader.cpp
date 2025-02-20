#include "Matrix/MtxReader.hpp"
#include "Common/StringFunc.hpp"
#include <fstream>
#include <string>
#include <vector>

MtxReader& MtxReader::Instance(){
    static MtxReader instance;
    return instance;
}

Graph MtxReader::Read(const char* path){
    std::string line;
    int n_vert = 0; // number of vertices
    int nnz = 0;    // number of non-zero

    std::ifstream file;
    file.open(path);


    // Read the header
    std::string header_object;
    std::string header_format;
    std::string header_field;
    std::string header_symmetry;
    {
        std::getline(file, line);
        std::vector<std::string> header = StringFunc::Instance().Parse(line.c_str());
        
        header_object   = header[1];
        header_format   = header[2];
        header_field    = header[3];
        header_symmetry = header[4];

        // Unsupported formats
        // Please read https://networkrepository.com/mtx-matrix-market-format.html
        // Following types are not supported
        if ( (strcmp(header_object.c_str(), "vector") == 0) ||
             (strcmp(header_field.c_str(), "complex") == 0) ||
             (strcmp(header_symmetry.c_str(), "hermitian") == 0) ){
            Graph g(0);
            return g;
        }
    }

    // Property
    while(std::getline(file, line)){
        const char* as_c_str = line.c_str();
        if (as_c_str[0] == '%') {
            continue; // This line is a comment
        }
        std::vector<std::string> line_parsed = StringFunc::Instance().Parse(as_c_str);
        int n_vert_check = 0;
        if (strcmp(header_format.c_str(), "coordinate") == 0 ){
            n_vert = std::stoi(line_parsed[0]);
            n_vert_check = std::stoi(line_parsed[1]);
            nnz = std::stoi(line_parsed[2]);
        }
        else{ //if (strcmp(header_format.c_str(), "array") == 0 ){
            n_vert = std::stoi(line_parsed[0]);
            n_vert_check = std::stoi(line_parsed[1]);
            nnz = INT_MAX;
        }
        if( n_vert_check > n_vert ) // If it is not square, force it to use bigger number of number of vertices
            n_vert = n_vert_check;
        break;    
    }

    bool is_pattern = strcmp(header_field.c_str(), "pattern") == 0;

    Graph g(n_vert);
    //Read out edges
    int edge_cnt = 0;
    while(std::getline(file, line) && (edge_cnt < nnz)){
        edge_cnt += 1;
        const char* as_c_str = line.c_str();
        std::vector<std::string> line_parsed = StringFunc::Instance().Parse(as_c_str);
        
        int i, j; //coordinates
        float value = 1;
        // Notice that mtx is one based matrix
        if (is_pattern) {
            i = std::stoi(line_parsed[0]) - 1;
            j = std::stoi(line_parsed[1]) - 1;
        }
        else{
            i = std::stoi(line_parsed[0]) - 1;
            j = std::stoi(line_parsed[1]) - 1;
            value = std::stof(line_parsed[2]);
        }
        g.Connect(i,j,value);
    }
    return g;
}