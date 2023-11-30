cmake_minimum_required(VERSION 3.10)

project(Supernodal_floyd_warshall)

#Define compilers
IF(NOT DEFINED CMAKE_C_COMPILER)
    SET(CMAKE_C_COMPILER clang)
ENDIF(NOT DEFINED CMAKE_C_COMPILER)

IF(NOT DEFINED CMAKE_CXX_COMPILER)
    SET(CMAKE_CXX_COMPILER clang++)
ENDIF(NOT DEFINED CMAKE_CXX_COMPILER)

SET(SRC_ROOT ${CMAKE_SOURCE_DIR}/src)

#find_package(CUDA)

LINK_LIBRARIES(${CMAKE_SOURCE_DIR}/modules/metis/build/Debug/metis.lib)
LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/modules/metis/build/Debug)
include_directories(${CMAKE_SOURCE_DIR}/modules/metis/build/src/include)

include_directories(${CMAKE_SOURCE_DIR}/include)

#Libraries
## Graph
add_library(AdjGraph ${SRC_ROOT}/Graph/AdjGraph.cpp)
add_library(CSRGraph ${SRC_ROOT}/Graph/CSRGraph.cpp)
add_library(Graph ${SRC_ROOT}/Graph/Graph.cpp)
add_library(Converter ${SRC_ROOT}/Graph/Converter.cpp)

## Heap
add_library(BinaryHeap ${SRC_ROOT}/Heap/BinaryHeap.cpp)

## SSSP
add_library(Bellman-Ford ${SRC_ROOT}/SSSP/Bellman-Ford.cpp)
add_library(Dijkstra ${SRC_ROOT}/SSSP/Dijkstra.cpp)


#Applications
add_executable(SSSP_Bellman ${SRC_ROOT}/App/SSSP.cpp)
add_executable(SSSP_Dijkstra ${SRC_ROOT}/App/SSSP.cpp)
target_link_libraries(SSSP_Bellman CSRGraph Graph Converter Bellman-Ford)
target_link_libraries(SSSP_Dijkstra CSRGraph Graph Converter BinaryHeap Dijkstra)

add_compile_definitions(-g)