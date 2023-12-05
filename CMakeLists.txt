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

## Common
add_library(StringFunc ${SRC_ROOT}/Common/StringFunc.cpp)

## Heap
add_library(BinaryHeap ${SRC_ROOT}/Heap/BinaryHeap.cpp)

## Graph
add_library(CSRGraph ${SRC_ROOT}/Graph/CSRGraph.cpp)
add_library(Graph ${SRC_ROOT}/Graph/Graph.cpp)
add_library(Converter ${SRC_ROOT}/Graph/Converter.cpp)

## Matrix
add_library(MtxReader ${SRC_ROOT}/Matrix/MtxReader.cpp)
target_link_libraries(MtxReader StringFunc)

## SSSP
add_library(Bellman-Ford-SSSP ${SRC_ROOT}/SSSP/Bellman-Ford.cpp)
target_link_libraries(Bellman-Ford-SSSP CSRGraph)

add_library(Dijkstra-SSSP ${SRC_ROOT}/SSSP/Dijkstra.cpp)
target_link_libraries(Dijkstra-SSSP CSRGraph BinaryHeap)

## APSP
add_library(Bellman-Ford-APSP ${SRC_ROOT}/APSP/Delegation.cpp)
target_link_libraries(Bellman-Ford-APSP CSRGraph Bellman-Ford-SSSP)

add_library(Dijkstra-APSP ${SRC_ROOT}/APSP/Delegation.cpp)
target_link_libraries(Dijkstra-APSP CSRGraph Dijkstra-SSSP)

add_library(FloydWarshall ${SRC_ROOT}/APSP/FloydWarshall.cpp)
target_link_libraries(FloydWarshall CSRGraph)

add_library(Johnson ${SRC_ROOT}/APSP/Johnson.cpp)
target_link_libraries(Johnson CSRGraph BinaryHeap)

#Applications
add_executable(SSSP_Bellman ${SRC_ROOT}/App/SSSP.cpp)
target_link_libraries(SSSP_Bellman CSRGraph Graph Converter MtxReader Bellman-Ford-SSSP)

add_executable(SSSP_Dijkstra ${SRC_ROOT}/App/SSSP.cpp)
target_link_libraries(SSSP_Dijkstra CSRGraph Graph Converter MtxReader  Dijkstra-SSSP)

add_executable(APSP_Bellman ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_Bellman CSRGraph Graph Converter MtxReader Bellman-Ford-APSP)

add_executable(APSP_Dijkstra ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_Dijkstra CSRGraph Graph Converter MtxReader Dijkstra-APSP)

add_executable(APSP_FloydWarshall ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_FloydWarshall CSRGraph Graph Converter MtxReader FloydWarshall)

add_executable(APSP_Johnson ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_Johnson CSRGraph Graph Converter MtxReader Johnson)

add_executable(Compare ${SRC_ROOT}/App/Compare.cpp)
target_link_libraries(Compare StringFunc)

add_executable(NegativeCycleCheck ${SRC_ROOT}/App/NegativeCycleCheck.cpp)
target_link_libraries(NegativeCycleCheck CSRGraph Graph Converter MtxReader FloydWarshall)

add_compile_definitions(-g)