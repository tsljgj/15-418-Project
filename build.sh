#!/bin/bash

# Create build directory if it doesn't exist
mkdir -p build

# Navigate to build directory
cd build

# Configure CMake (point to src directory)
# No need to specify -lpthread as it will be handled properly per platform
cmake ../src 

# Build the project
cmake --build .

echo "Build complete. Executable is in the build directory."