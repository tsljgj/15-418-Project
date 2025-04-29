@echo off
setlocal enabledelayedexpansion

:: Define algorithms and P:E ratios
set ALGORITHMS=naive,static,static_bl,vfc
set RATIOS=1:2,1:6,2:4,3:2,2:12,4:8,6:4

:: Find all input files
echo Finding input files in the inputs directory...
set INPUT_FILES=

for %%F in (inputs\*.txt) do (
    set INPUT_FILES=!INPUT_FILES! %%F
)

echo Found the following input files: %INPUT_FILES%

:: Run the checker for each input file with the specified parameters
for %%F in (%INPUT_FILES%) do (
    echo.
    echo ===========================================
    echo Processing file: %%F
    echo ===========================================
    python src/checker.py %%F --algorithm %ALGORITHMS% --p-e-ratio %RATIOS%
    echo.
)

echo All tests completed.