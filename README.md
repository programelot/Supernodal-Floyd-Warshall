# Supernodal-Floyd-Warshall

# Introduction
The Supernodal Floyd-Warshall algorithm is an algorithm that can solve all-pair-shortest-path in a sparse graph.\
I made this in 2020 the first semester of my graduate university project.\
The code is so un-organized and complex so I reconstructed it from the ground again at the end of 2023.\
This project includes many other graph-related algorithms like Dijkstra's algorithm.\
For a traveler who travels the optimization field, I left this implementation in here.

# Prerequirement

I rebuilt this on Windows desktop.\
Therefore, I installed LLVM in the local desktop.\
Please use following commands to install it.\
This project uses openMP either.\
So, it must be enabled.

```bash
git clone https://github.com/llvm/llvm-project.git
cd .\llvm-project\
git checkout release/17.x
cmake -S llvm -B build_release -DCMAKE_BUILD_TYPE=Release -DLLVM_ENABLE_PROJECTS="clang;clang-tools-extra;compiler-rt;mlir;polly;lldb;openmp" -DLLVM_TARGETS_TO_BUILD=X86 -Thost=x64
cmake --build .\build_release --target INSTALL --config Release
```

# Implementation
Detailed implementation will be uploaded soon.

# Build and test

```bash
git clone https://github.com/programelot/Supernodal-Floyd-Warshall.git
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_THREAD=Y ..
cmake --build .
./Debug/APSP_SupernodalFloydWarshall.exe path_to_delaunay_n13.mtx path_to_output.txt    
```


# Source
* APP : Directory includes all applications
    * APSP.cpp : All-source-shortest-path problem test
    * Compare.cpp : Compare two outputs of APSP
    * NegativeCycleCheck.cpp : Check the negative cycle of the graph
    * SSSP.cpp : Single-source-shortest-path problem test
* APSP : APSP related libraries
    * BlockedFloydWarshall.cpp : Blocked Floyd-Warshall algorithm
    * BlockedFloydWarshallCompressive.cpp : Blocked Floyd-Warshall algorithm with compressive edge counting
    * Delegation.cpp : Delegation algorithm (use SSSP to solve APSP)
    * FloydWarshall.cpp : Normal Floyd-Warshall algorithm
    * FloydWarshallCompressive.cpp : Normal Floyd-Warshall algorithm with compressive edge counting
    * FloydWarshallDualBuffer.cpp : Normal Floyd-Warshall algorithm with duble buffering
    * FloydWarshallDualCompInOut.cpp : Normal Floyd-Warshall algorithm with duble buffering and compressive edge counting for both In&Out edges
    * FloydWarshallDualCompressive.cpp :  Normal Floyd-Warshall algorithm with duble buffering and compressive edge counting for only outgoing edges
    * Johnson.cpp : Johnson's algorithm
    * SupernodalFloydWarshall.cpp : Supernodal FloydWarshall algorithm
    * SupernodalFloydWarshallComp.cpp : Supernodal FloydWarshall algorithm with compressive edge counting
* Common : Common librarys for entire projects
    * DebugAssert.cpp : Aseert for debug mode.
    * StringFunc.cpp : String related functions.
* ETree
    * ETree.cpp : Elimination tree library
* Graph : Graph related libraries
    * Converter.cpp : Converter between regular graph and CSR graph
    * CSRGraph.cpp : CSR graph format
    * Graph.cpp : Regular graph format
    * SortCSRRow.cpp : Simple sorting algorithm for a row of CSR graph.
* Heap
    * BinaryHeap.cpp : Binary heap library
* Image
    * BMP.cpp : BMP drawing library
* Matrix
    * MtxReader.cpp : MtxFileReader
* SSSP
    * Bellman-Ford.cpp : Bellman-Ford algorithm
    * Dijkstra.cpp : Dijkstra algorithm

# Result (sec)

* Input : delaunay_n13
* CPU : I7-12700F 2.10 GHz
* RAM : 32GB
* OS : Windows 10

## Delegation
### Bellman-ford algorithm
1. Single thread : 30.212
2. Multiple threads : 2.373

### Dijkstra algorithm
1. Single thread : 13.846
2. Multiple threads : 1.126

## Johnson's algorithm
1. Single thread: 14.137
2. Multiple threads: 1.275

## Floyd-Warshall algorithm
1. Single thread: 363.867

## Floyd-Warshall algorithm (With Compressive edge counting)
1. Single thread: 414.157

## Floyd-Warshall algorithm (With Dual buffering)
1. Single thread: 544.475
2. Multiple threads: 128.34

## Floyd-Warshall algorithm (With Dual buffering and compressive outgoing edge counting)
1. Single thread: 483.309
2. Multiple threads: 95.764

## Floyd-Warshall algorithm (With Dual buffering and compressive In&Out edge counting)
1. Single thread: 552.198
2. Multiple threads: 80.538

## Blocked Floyd-Warshall algorithm
1. Single thread: 409.435
2. Multiple threads: 54.791

## Blocked Floyd-Warshall algorithm (With Compressive edge counting)
1. Single thread: 367.1
2. Multiple threads: 48.478

## Supernodal Floyd-Warshall
1. Single thread: 31.285
2. Multiple threads: 4.646

## Supernodal Floyd-Warshall (Without compressive edge counting)
1. Single thread: 42.091
2. Multiple threads: 5.948

# Notes
## Graph
There are a lot of ways to implement graphs.\
I used an adjacent matrix and CSR in this project.\
Input graph needs to be formed in CSR format.\
Output of the algorithm is dense matrix so it will be an adjacent matrix.\
Only exception will be a list that will be used for other sub algorithms.

## Future works
Algorithm achieved big speed-up incomparison with classic Floyd-Warshall algorithm.\
However, it was slower than Johnson's algorithm.\
I only did a comparison it with the classic Floyd-Warshall algorithm before.\
It made more sight by comparing it with Johnson's algorithm.\
There could be a better implementation since I made some part of it by my own.\
I made some operations like "Finding seperators in the nested dissection" by my own.\
I will optimize it more when time allows and method founds.

# Original source
The original source I made in 2020 can be found in the archive branch.\
It also includes some documents that I used for the lecture.\
https://github.com/programelot/Supernodal-Floyd-Warshall/tree/archived

# Public paper
**I am not an author of the Paper.** \
https://dl.acm.org/doi/abs/10.1145/3332466.3374533
