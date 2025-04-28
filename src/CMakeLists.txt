cmake_minimum_required(VERSION 3.10)
project(LouvainAlgorithm VERSION 1.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Add optimization flags - use appropriate flags for MinGW
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# Set output directories to put executables directly in build folder
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR})

# Main Louvain algorithm executable
add_executable(test_louvain
    test_louvain.cpp
    louvain_seq.cpp
    louvain_parallel.cpp
    louvain_parallel_bl.cpp
    louvain_parallel_static.cpp
    louvain_parallel_static_bl.cpp
    louvain_parallel_vfc.cpp
    louvain_parallel_vfc_bl.cpp
    graph.cpp
    core_type.cpp
)