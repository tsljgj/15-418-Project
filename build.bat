@echo off
REM Create build directory if it doesn't exist
if not exist build mkdir build

REM Navigate to build directory
cd build

REM Configure CMake (point to src directory)
cmake ..\src

REM Build the project
cmake --build . --config Release

echo Build complete. Executable is in the build\Release directory.