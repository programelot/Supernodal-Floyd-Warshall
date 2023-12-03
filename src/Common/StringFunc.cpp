// String parser library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#include <vector>
#include <string>
#include "Common/StringFunc.hpp"

StringFunc& StringFunc::Instance(){
    static StringFunc instance;
    return instance;
}

bool StringFunc::isWhitespace(char c){
    switch(c){
        default: 
            return false;
        case '\n':
        case ' ': 
            return true;
    }
}

std::vector<std::string> StringFunc::Parse(const std::string &str){
    std::vector<std::string> parsed;
    size_t cut = 0;
    while(isWhitespace(str[cut])){
        ++cut;
    }
    for(size_t i = cut + 1; i < str.size(); ++i){
        if(isWhitespace(str[i])){
            if(cut == i){//Sequenced white space will be ignored.
                ++cut;
                continue; 
            } 
            parsed.emplace_back(str.substr(cut, i - cut));
            cut = i + 1;
        }
    }
    if(cut != str.size()){
        parsed.emplace_back(str.substr(cut, str.size() - cut));
    }
    return parsed;
}

bool StringFunc::StrCmp(const std::string& lhs, const std::string& rhs){
    if(lhs.size() != rhs.size()) return false;
    for(size_t i = 0; i < lhs.size(); ++i)
        if(lhs[i] != rhs[i]) return false;
    return true;
}