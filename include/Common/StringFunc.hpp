// String parser library //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#ifndef INCLUDE_STRINGPARSER_HPP
#define INCLUDE_STRINGPARSER_HPP

#include <vector>
#include <string>

class StringFunc{
private:
    StringFunc(){}//Singleton
    bool isWhitespace(char c);
public:
    static StringFunc& Instance();
    std::vector<std::string> Parse(const std::string &str);
    bool StrCmp(const std::string& lhs, const std::string& rhs);
};
#endif