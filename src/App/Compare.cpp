#include <fstream>
#include "Common/DebugAssert.hpp"
#include <string>
#include <iostream>
#include "Common/StringFunc.hpp"

namespace{

    float max(float a, float b){
        return a > b? a : b;
    }
}

int main(int argc, const char* argv[]){
    if(argc != 3){
        printf("Usage: Compare.exe file1 file2");
    }
    std::ifstream file1, file2;
    printf("Open file %s\n", argv[1]);
    file1.open(argv[1]);
    AssertDebug(__FILE__, __LINE__,"File can not be opend", file1.is_open());

    printf("Open file %s\n", argv[2]);
    file2.open(argv[2]);
    AssertDebug(__FILE__, __LINE__,"File can not be opend", file2.is_open());

    float error = 0;
    std::string value1, value2;
    while(!file1.eof()){
        getline(file1, value1, ' ');
        if(file2.eof()){
            printf("Different size of result\n");
            break;
        }
        getline(file2, value2, ' ');
        if(value1.empty() || value2.empty()){
            printf("empty line deceted\n");
            break;
        }
        // printf("%d : ", value1.size());
        // for(int i = 0; i < value1.size(); i++){
        //     printf("%d ", value1[i]);
        // }printf(" // ");
        // printf("%d : ", value2.size());
        // for(int i = 0; i < value2.size(); i++){
        //     printf("%d ", value2[i]);
        // }printf("\n");
        if(StringFunc::Instance().StrCmp(value1, value2))
            continue;
        float f1 = std::stof(value1);
        float f2 = std::stof(value2);
        float diff = abs(f1 - f2);
        float max = abs(f1) > abs(f2) ? abs(f1): abs(f2);
        error += diff/max;
    }
    printf("Error detected : %f %%\n",error * 100);
    
    file1.close();
    file2.close();
}