cmake_minimum_required(VERSION 3.10)

project(Supernodal_floyd_warshall)

#Define compilers
IF(NOT DEFINED CMAKE_C_COMPILER)
    SET(CMAKE_C_COMPILER clang)
ENDIF(NOT DEFINED CMAKE_C_COMPILER)

IF(NOT DEFINED CMAKE_CXX_COMPILER)
    SET(CMAKE_CXX_COMPILER clang++)
ENDIF(NOT DEFINED CMAKE_CXX_COMPILER)

IF(NOT DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE "Debug")
ENDIF(NOT DEFINED CMAKE_BUILD_TYPE)

IF(NOT DEFINED USE_THREAD)
    SET(USE_THREAD "N")
ENDIF(NOT DEFINED USE_THREAD)

SET(SRC_ROOT ${CMAKE_SOURCE_DIR}/src)

#find_package(CUDA)

LINK_LIBRARIES(${CMAKE_SOURCE_DIR}/modules/metis/build/Debug/metis.lib)
LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/modules/metis/build/Debug)
include_directories(${CMAKE_SOURCE_DIR}/modules/metis/build/src/include)

include_directories(${CMAKE_SOURCE_DIR}/include)

#Libraries

## Common
add_library(StringFunc ${SRC_ROOT}/Common/StringFunc.cpp)
add_library(DebugAssert ${SRC_ROOT}/Common/DebugAssert.cpp)

## Image
add_library(BMP ${SRC_ROOT}/Image/BMP.cpp)

## Heap
add_library(BinaryHeap ${SRC_ROOT}/Heap/BinaryHeap.cpp)
target_link_libraries(BinaryHeap DebugAssert)

## Graph
add_library(SortCSRRow ${SRC_ROOT}/Graph/SortCSRRow.cpp)
add_library(CSRGraph ${SRC_ROOT}/Graph/CSRGraph.cpp)
add_library(Graph ${SRC_ROOT}/Graph/Graph.cpp)
add_library(Converter ${SRC_ROOT}/Graph/Converter.cpp)
target_link_libraries(Converter DebugAssert SortCSRRow)

## Supernod
add_library(ETree ${SRC_ROOT}/ETree/ETree.cpp)

## Matrix
add_library(MtxReader ${SRC_ROOT}/Matrix/MtxReader.cpp)
target_link_libraries(MtxReader StringFunc DebugAssert)

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

add_library(FloydWarshallCompressive ${SRC_ROOT}/APSP/FloydWarshallCompressive.cpp)
target_link_libraries(FloydWarshallCompressive CSRGraph)

add_library(FloydWarshallDualBuffer ${SRC_ROOT}/APSP/FloydWarshallDualBuffer.cpp)
target_link_libraries(FloydWarshallDualBuffer CSRGraph)

add_library(FloydWarshallDualCompressive ${SRC_ROOT}/APSP/FloydWarshallDualCompressive.cpp)
target_link_libraries(FloydWarshallDualCompressive CSRGraph)

add_library(FloydWarshallDualCompInOut ${SRC_ROOT}/APSP/FloydWarshallDualCompInOut.cpp)
target_link_libraries(FloydWarshallDualCompInOut CSRGraph)

add_library(BlockedFloydWarshall ${SRC_ROOT}/APSP/BlockedFloydWarshall.cpp)
target_link_libraries(BlockedFloydWarshall CSRGraph)

add_library(BlockedFloydWarshallCompressive ${SRC_ROOT}/APSP/BlockedFloydWarshallCompressive.cpp)
target_link_libraries(BlockedFloydWarshallCompressive CSRGraph)

add_library(Johnson ${SRC_ROOT}/APSP/Johnson.cpp)
target_link_libraries(Johnson CSRGraph BinaryHeap)

add_library(SupernodalFloydWarshall ${SRC_ROOT}/APSP/SupernodalFloydWarshall.cpp)
target_link_libraries(SupernodalFloydWarshall ETree)

# Applications
## SSSP
add_executable(SSSP_Bellman ${SRC_ROOT}/App/SSSP.cpp)
target_link_libraries(SSSP_Bellman CSRGraph Graph Converter MtxReader Bellman-Ford-SSSP StringFunc)

add_executable(SSSP_Dijkstra ${SRC_ROOT}/App/SSSP.cpp)
target_link_libraries(SSSP_Dijkstra CSRGraph Graph Converter MtxReader  Dijkstra-SSSP StringFunc)

##APSP
add_executable(APSP_Bellman ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_Bellman CSRGraph Graph Converter MtxReader Bellman-Ford-APSP StringFunc)

add_executable(APSP_Dijkstra ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_Dijkstra CSRGraph Graph Converter MtxReader Dijkstra-APSP StringFunc)

add_executable(APSP_Johnson ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_Johnson CSRGraph Graph Converter MtxReader Johnson StringFunc)

add_executable(APSP_FloydWarshall ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_FloydWarshall CSRGraph Graph Converter MtxReader FloydWarshall StringFunc)

add_executable(APSP_FloydWarshallCompressive ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_FloydWarshallCompressive CSRGraph Graph Converter MtxReader FloydWarshallCompressive StringFunc)

add_executable(APSP_FloydWarshallDualBuffer ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_FloydWarshallDualBuffer CSRGraph Graph Converter MtxReader FloydWarshallDualBuffer StringFunc)

add_executable(APSP_FloydWarshallDualCompressive ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_FloydWarshallDualCompressive CSRGraph Graph Converter MtxReader FloydWarshallDualCompressive StringFunc)

add_executable(APSP_FloydWarshallDualCompInOut ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_FloydWarshallDualCompInOut CSRGraph Graph Converter MtxReader FloydWarshallDualCompInOut StringFunc)

add_executable(APSP_BlockedFloydWarshall ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_BlockedFloydWarshall CSRGraph Graph Converter MtxReader BlockedFloydWarshall StringFunc)

add_executable(APSP_BlockedFloydWarshallCompressive ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_BlockedFloydWarshallCompressive CSRGraph Graph Converter MtxReader BlockedFloydWarshallCompressive StringFunc)

add_executable(APSP_SupernodalFloydWarshall ${SRC_ROOT}/App/APSP.cpp)
target_link_libraries(APSP_SupernodalFloydWarshall CSRGraph Graph Converter MtxReader SupernodalFloydWarshall StringFunc)

## Tools
add_executable(Compare ${SRC_ROOT}/App/Compare.cpp)
target_link_libraries(Compare StringFunc DebugAssert)

add_executable(NegativeCycleCheck ${SRC_ROOT}/App/NegativeCycleCheck.cpp)
target_link_libraries(NegativeCycleCheck CSRGraph Graph Converter MtxReader Bellman-Ford-APSP)

add_executable(PrintETree ${SRC_ROOT}/App/PrintETree.cpp)
target_link_libraries(PrintETree CSRGraph Graph Converter ETree MtxReader BMP)

IF (${USE_THREAD} STREQUAL "Y")
    add_definitions(-DUSE_THREAD)
    add_definitions(-openmp)

    target_link_libraries(NegativeCycleCheck CSRGraph Graph Converter MtxReader FloydWarshall)
ELSE()
ENDIF (${CMAKE_BUILD_TYPE} STREQUAL "Y")

add_compile_definitions(-O3)
IF (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    add_compile_definitions(-g)
    add_definitions(-DDEBUG)
ENDIF (${CMAKE_BUILD_TYPE} STREQUAL "Debug")