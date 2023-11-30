# Supernodal-Floyd-Warshall

# Introduction
The Supernodal Floyd-Warshall algorithm is an algorithm that can solve all-pair-shortest-path in a sparse graph.\
I made this in 2020 the first semester of my graduate university project.\
The code is so un-organized and complex so I reconstructed it from the ground again at the end of 2023.\
This project includes many other graph-related algorithms like Dijkstra's algorithm.\
For a traveler who travels the optimization field, I left this implementation in here.

# Prerequirement

I rebuilt this on Windows desktop.
Therefore, I installed LLVM in the local desktop.
Please use following commands to install it.

```bash
git clone https://github.com/llvm/llvm-project.git
cd .\llvm-project\
git checkout release/17.x
cmake -S llvm -B build_release -DCMAKE_BUILD_TYPE=Release -DLLVM_ENABLE_PROJECTS="clang;clang-tools-extra;compiler-rt;mlir;polly;lldb" -DLLVM_TARGETS_TO_BUILD=X86 -Thost=x64
cmake --build .\build_release --target INSTALL --config Release
```


# Notes
## Graph
There are a lot of ways to implement graphs.\
I used Adjacent matrix and CSR in this project.\
Input graph need to be formed in CSR format.\
Output of the algorithm is dense matrix so it will be adjacent matrix.\
Only exception will be a list that will be used for other sub algorithms.

## Improvement
There could be a better implementation since I am not metis-friendly.\
I made some parts of the operation in my own way of estimation. (Finding nested Dissection)\
I will fix it when I can do it later.

# Original source
The original source I made in 2020 can be found in the archive branch.\
It also includes some documents that I used for the lecture.\
https://github.com/programelot/Supernodal-Floyd-Warshall/tree/archived

# Public paper
**I am not an author of the Paper.** \
https://dl.acm.org/doi/abs/10.1145/3332466.3374533
