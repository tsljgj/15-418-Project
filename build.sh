#!/bin/bash

# Create build directory if it doesn't exist
mkdir -p build

# Navigate to build directory
cd build

# Configure CMake (point to src directory)
cmake ../src

# Build the project
cmake --build .

echo "Build complete. Executable is in the build directory."