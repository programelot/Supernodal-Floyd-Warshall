#include <iostream>

#include "Graph/Converter.hpp"
#include "Graph/Graph.hpp"
#include "Graph/CSRGraph.hpp"
#include "Matrix/MtxReader.hpp"
#include "ETree/ETree.hpp"
#include "Image/BMP.hpp"

int min(int a, int b){
    return a < b? a : b;
}

enum class DrawType{
    original, permutated, eTree, analysis, error
};

int main(int argc, char* argv[]){
    unsigned char colorStrength = 255;
    unsigned char gray = 200;
    unsigned char darkGray = 125;

    if(argc != 4){
        printf("Usage : PrintEtree option MatrixFile Output\n");
        return 0;
    }
    DrawType drawType = DrawType::error;
    bool drawOriginal = false;
    bool drawDiagonal = true;
    bool drawFirstSeeker = false;
    bool drawSeperator = false;
    bool drawNoneZero = false;
    bool drawPossibleBox = false;
    bool drawTree = false;

    if(std::strcmp(argv[1], "org") == 0){
        drawType = DrawType::original;
    }
    else if(std::strcmp(argv[1], "perm") == 0){
        drawType = DrawType::permutated;
    }
    else if(std::strcmp(argv[1], "etree") == 0){
        drawType = DrawType::eTree;
    }
    else if(std::strcmp(argv[1], "analysis") == 0){
        drawType = DrawType::analysis;
    }

    if(drawType == DrawType::error){
        printf("Possible option need to be one of belows.\n");
        printf("org : original (not permutated) graph.\n");
        printf("perm : permutated graph.\n");
        printf("etree : permutated graph with only seperators and supernodes informations.\n");
        printf("analysis : permutated graph with all information related with elimination tree.\n");
        return 0;
    }
    switch(drawType){
        case DrawType::original:
            drawOriginal = true;
            drawDiagonal = true;
            drawFirstSeeker = false;
            drawSeperator = false;
            drawNoneZero = false;
            drawPossibleBox = false;
            drawTree = false;
            break;
        case DrawType::permutated:
            drawOriginal = false;
            drawDiagonal = true;
            drawFirstSeeker = false;
            drawSeperator = false;
            drawNoneZero = false;
            drawPossibleBox = false;
            drawTree = false;
            break;
        case DrawType::eTree:
            drawOriginal = false;
            drawDiagonal = true;
            drawFirstSeeker = false;
            drawSeperator = true;
            drawNoneZero = true;
            drawPossibleBox = false;
            drawTree = true;
            break;
        case DrawType::analysis:
            drawOriginal = false;
            drawDiagonal = true;
            drawFirstSeeker = true;
            drawSeperator = true;
            drawNoneZero = true;
            drawPossibleBox = true;
            drawTree = false;
            break;
    }
    Graph g = MtxReader::Instance().Read(argv[2]);
    CSRGraph csr = Converter::Instance().ToCSR(g);
    dataSize_t size = csr.Size();
    Bitmap bmp(argv[3], csr.Size(), csr.Size());
    printf("Drawing : ");
    switch(drawType){
        case DrawType::original:
        printf("Original graph\n");
        break;
        case DrawType::permutated:
        printf("Permutated graph\n");
        break;
        case DrawType::eTree:
        printf("Permutated graph in eTree\n");
        break;
        case DrawType::analysis:
        printf("Permutated graph with analysis informations\n");
        break;
    }
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            bmp[i * size + j].r = 255;
            bmp[i * size + j].g = 255;
            bmp[i * size + j].b = 255;
        }
    }
    
    ETree eTree(csr);
    const dataSize_t* perm = eTree.Perm();
    const dataSize_t* iperm = eTree.IPerm();
    const CSRGraph& csrPerm = eTree.PermGraph();

    const std::vector<Supernode>& supernodes = eTree.Supernodes();
    const std::vector<dataSize_t>& depthIndex = eTree.DepthIndex();
    const dataSize_t* maxBy = eTree.MaxBy();
    std::vector<dataSize_t> matSize;//Size of entire sub matrix
    matSize.resize(supernodes.size());

    //Draw empty space height map
    if(drawFirstSeeker){
        for(int i = 0; i < size; ++i){
            for(int j = 0; j < csrPerm.ColIdx()[csrPerm.RowPtr()[i]]; ++j){
                bmp[i * size + j].r = gray;
                bmp[i * size + j].g = gray;
                bmp[i * size + j].b = gray;
            }
            for(int j = 0; j < eTree.Height()[i]; ++j){
                bmp[j * size + i].r = gray;
                bmp[j * size + i].g = gray;
                bmp[j * size + i].b = gray;
            }
        }
    }

    //Draw empty space maxBy
    if(drawPossibleBox){
        for(int i = 0; i < size; ++i){
            for(int j = 0; j < min(maxBy[i] - i, i + 1); ++j){
                bmp[(i - j) * size + i].r = darkGray;
                bmp[(i - j) * size + i].g = darkGray;
                bmp[(i - j) * size + i].b = darkGray;

                bmp[i * size + (i - j)].r = darkGray;
                bmp[i * size + (i - j)].g = darkGray;
                bmp[i * size + (i - j)].b = darkGray;
            }
        }
    }

    //Draw diagonal line
    if(drawDiagonal){
        for(int i = 0; i < size; ++i){
            bmp[i * size + i].r = 0;
            bmp[i * size + i].g = colorStrength;
            bmp[i * size + i].b = 0;
        }
    }

    //Draw supernodes
    for(int i = depthIndex.size() - 2; i >= 0; --i){
        for(int j = depthIndex[i]; j < depthIndex[i + 1]; ++j){
            if(supernodes[j].left != size)// && upernodes[j].right != size
            {//Seperator
                int from = supernodes[j].from;
                int to = supernodes[j].to;
                matSize[j] = to - from + matSize[supernodes[j].left] + matSize[supernodes[j].right];
                if(drawSeperator){
                    for(int k = 0; k < matSize[j]; ++k){
                        bmp[(to - 1 - k) * size + (to - 1)].r = 0;
                        bmp[(to - 1 - k) * size + (to - 1)].g = 0;
                        bmp[(to - 1 - k) * size + (to - 1)].b = colorStrength;

                        bmp[(to - 1) * size + (to - 1 - k)].r = 0;
                        bmp[(to - 1) * size + (to - 1 - k)].g = 0;
                        bmp[(to - 1) * size + (to - 1 - k)].b = colorStrength;
                    }
                    
                    for(int k = from; k < to; ++k){
                        bmp[k * size + (to - matSize[j])].r = 0;
                        bmp[k * size + (to - matSize[j])].g = 0;
                        bmp[k * size + (to - matSize[j])].b = colorStrength;

                        bmp[(to - matSize[j]) * size + k].r = 0;
                        bmp[(to - matSize[j]) * size + k].g = 0;
                        bmp[(to - matSize[j]) * size + k].b = colorStrength;
                    }

                    for(int k = 0; k < matSize[j] - (to - from); ++k){
                        bmp[(from - k) * size + from].r = 0;
                        bmp[(from - k) * size + from].g = 0;
                        bmp[(from - k) * size + from].b = colorStrength;

                        bmp[from * size + (from - k)].r = 0;
                        bmp[from * size + (from - k)].g = 0;
                        bmp[from * size + (from - k)].b = colorStrength;
                    }
                }
            }
            else{//Not seperators
                int from = supernodes[j].from;
                int to = supernodes[j].to;
                matSize[j] = to - from;
                if(drawNoneZero){
                    for(int k = from; k < to; ++k){
                        bmp[k * size + from].r = 0;
                        bmp[k * size + from].g = 0;
                        bmp[k * size + from].b = colorStrength;

                        bmp[k * size + (to - 1)].r = 0;
                        bmp[k * size + (to - 1)].g = 0;
                        bmp[k * size + (to - 1)].b = colorStrength;

                        bmp[from * size + k].r = 0;
                        bmp[from * size + k].g = 0;
                        bmp[from * size + k].b = colorStrength;

                        bmp[(to - 1) * size + k].r = 0;
                        bmp[(to - 1) * size + k].g = 0;
                        bmp[(to - 1) * size + k].b = colorStrength;
                    }
                }
            }
            //Drwa Tress tructure
            if(drawTree){
                if(supernodes[j].parent == size){
                    continue;
                }
                int from = supernodes[j].to - matSize[j];
                int to = supernodes[j].to;
                for(int k = from; k < to; ++k){
                    bmp[k * size + from].r = 0;
                    bmp[k * size + from].g = colorStrength;
                    bmp[k * size + from].b = 0;
                    
                    bmp[k * size + (to - 1)].r = 0;
                    bmp[k * size + (to - 1)].g = colorStrength;
                    bmp[k * size + (to - 1)].b = 0;

                    bmp[from * size + k].r = 0;
                    bmp[from * size + k].g = colorStrength;
                    bmp[from * size + k].b = 0;
                    
                    bmp[(to - 1) * size + k].r = 0;
                    bmp[(to - 1) * size + k].g = colorStrength;
                    bmp[(to - 1) * size + k].b = 0;
                }
            }

        }
    }

    //Draw points
    if(drawOriginal){
        for(int i = 0; i < size; ++i){
            for(int j = csr.RowPtr()[i]; j < csr.RowPtr()[i + 1]; ++j){
                bmp[i * size + csr.ColIdx()[j]].r = colorStrength;
                bmp[i * size + csr.ColIdx()[j]].g = 0;
                bmp[i * size + csr.ColIdx()[j]].b = 0;
            }
        }
    }
    else{
        for(int i = 0; i < size; ++i){
            for(int j = csrPerm.RowPtr()[i]; j < csrPerm.RowPtr()[i + 1]; ++j){
                bmp[i * size + csrPerm.ColIdx()[j]].r = colorStrength;
                bmp[i * size + csrPerm.ColIdx()[j]].g = 0;
                bmp[i * size + csrPerm.ColIdx()[j]].b = 0;
            }
        }
    }
}