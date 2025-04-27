#!/bin/bash

# Loop through all .txt files in the inputs/ directory
for file in inputs/*.txt
do
    echo "Running checker on $file"
    python3 src/checker.py "$file" --algorithm sequential,naive,vfc --threads 2,4,8 --runs 1
done
