#include "Common/DebugAssert.hpp"
#include <cassert>
#include <cstdio>
#include <cstdlib>

void AssertDebug(const char* file, int line, const char* message, bool test){

#ifdef DEBUG
    if(test) return;
    fprintf(stderr, "Assertion failed at %s, line %d \n %s\n", file, line, message);
    std::abort();
#endif

}