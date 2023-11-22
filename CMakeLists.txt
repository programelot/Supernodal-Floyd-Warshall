cmake_minimum_required(VERSION 3.10)

set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")
find_package(CUDA)

LINK_LIBRARIES(papi)
LINK_DIRECTORIES(/opt/papi/5.7.1.0/lib)
include_directories(/opt/papi/5.7.1.0/include)

LINK_LIBRARIES(/opt/metis/lib/libmetis.a) #LINK_LIBRARIES(/scratch/sunghyunmo/metis/lib/libmetis.a)
LINK_DIRECTORIES(/opt/metis/include/lib) #LINK_DIRECTORIES(/scratch/sunghyunmo/metis/include/lib)
include_directories(/opt/metis/include) #include_directories(/scratch/sunghyunmo/metis/include)

SET(CMAKE_CXX_FLAGS -pthread)

add_executable(main main.cpp)
add_compile_options(main -lpthread -lpathDoubling -O3 -std=c+17)
target_link_libraries(main pathDoubling)
cuda_add_library(pathDoubling pathDoubling.cu)
