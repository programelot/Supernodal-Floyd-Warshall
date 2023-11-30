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
add_library(CSRGraph ${SRC_ROOT}/Graph/CSRGraph.cpp)
add_library(Graph ${SRC_ROOT}/Graph/Graph.cpp)
add_library(Converter ${SRC_ROOT}/Graph/Converter.cpp)

## Heap
add_library(BinaryHeap ${SRC_ROOT}/Heap/BinaryHeap.cpp)

## SSSP
add_library(Bellman-Ford-SSSP ${SRC_ROOT}/SSSP/Bellman-Ford.cpp)
target_link_libraries(Bellman-Ford-SSSP CSRGraph)

add_library(Dijkstra-SSSP ${SRC_ROOT}/SSSP/Dijkstra.cpp)
target_link_libraries(Dijkstra-SSSP CSRGraph BinaryHeap)

## APSP
add_library(Bellman-Ford-APSP ${SRC_ROOT}/APSP/Delegation.cpp)
target_link_libraries(Bellman-Ford-APSP CSRGraph Bellman-Ford-SSSP)

add_library(Dijkstra-APSP ${SRC_ROOT}/APSP/Delegation.cpp)
target_link_libraries(Dijkstra-APSP CSRGraph BinaryHeap Dijkstra-SSSP)

#Applications
add_executable(SSSP_Bellman ${SRC_ROOT}/App/SSSP.cpp)
target_link_libraries(SSSP_Bellman CSRGraph Graph Converter Bellman-Ford-SSSP)

add_executable(SSSP_Dijkstra ${SRC_ROOT}/App/SSSP.cpp)
target_link_libraries(SSSP_Dijkstra CSRGraph Graph Converter  Dijkstra-SSSP)

add_executable(APSP_Bellman ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_Bellman CSRGraph Graph Converter Bellman-Ford-APSP)

add_executable(APSP_Dijkstra ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_Dijkstra CSRGraph Graph Converter Dijkstra-APSP)

add_compile_definitions(-g)