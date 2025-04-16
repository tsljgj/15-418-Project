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

# Add optimization flags
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# Main Louvain algorithm executable
add_executable(test_louvain
    test_louvain.cpp
    louvain_seq.cpp
    louvain_parallel.cpp
)
target_link_libraries(test_louvain PRIVATE OpenMP::OpenMP_CXX)