@echo off
REM Create build directory if it doesn't exist
if not exist build mkdir build

REM Navigate to build directory
cd build

REM Configure CMake with MinGW generator (point to src directory)
cmake -G "MinGW Makefiles" ..\src

REM Build the project
cmake --build . --config Release

REM Inform about executable location

cd ..
echo Build complete. 
echo Executable is located at: %cd%\build\test_louvain.exe