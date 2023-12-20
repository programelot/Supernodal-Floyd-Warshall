#ifndef INCLUDE_TYPE_HPP
#define INCLUDE_TYPE_HPP

#include <limits>
#include <stdint.h>

using weight_t = float;
using dataSize_t = int32_t;
constexpr weight_t kWeightInf = std::numeric_limits<float>::infinity();

#endif