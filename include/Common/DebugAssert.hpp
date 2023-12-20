// Assert library only actives it when it is on debug mode. //
// Author : Hyunmo Sung //
// Contact: programelot@gmail.com //

#ifndef INCLUDE_DEBUGASSERT_HPP
#define INCLUDE_DEBUGASSERT_HPP

void DebugAssert(const char* file, int line, const char* message, bool test);

#endif