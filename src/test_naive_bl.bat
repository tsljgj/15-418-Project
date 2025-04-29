@echo off
setlocal enabledelayedexpansion

:: Define P:E ratios to test
set RATIOS=1:2 1:6 2:4 3:2 2:12 4:8 6:4

:: Find all input files
echo Finding input files in the inputs directory...
set INPUT_FILES=

for %%F in (inputs\*.txt) do (
    set INPUT_FILES=!INPUT_FILES! %%F
)

echo Found the following input files: %INPUT_FILES%
echo.

:: Run the test for each input file with all ratios
for %%F in (%INPUT_FILES%) do (
    echo ===========================================
    echo Processing file: %%F
    echo ===========================================
    
    for %%R in (%RATIOS%) do (
        for /f "tokens=1,2 delims=:" %%A in ("%%R") do (
            set P_CORES=%%A
            set E_CORES=%%B
            
            echo Testing with !P_CORES! P-cores and !E_CORES! E-cores...
            build\test_louvain.exe %%F -B -pc !P_CORES! -ec !E_CORES!
            echo.
        )
    )
    
    echo.
)

echo All tests completed.