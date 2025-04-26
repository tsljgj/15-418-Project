# Louvain Community Detection

This repository contains sequential and parallel implementations of the Louvain community detection algorithm for large graphs.

## Overview

The Louvain algorithm is a hierarchical community detection method that optimizes modularity. This implementation includes:

- Sequential version (with P-core and E-core options)
- Naive parallel version using graph partitioning
- Optimized parallel version with Vertex Following and Coloring (VFC)

## Building the Project

### Windows (with MinGW)

Run the build script to compile the project:

```bash
build.bat
```

This will create the executable at `build/test_louvain.exe`.

## Running the Algorithm

You can run the algorithm in sequential or parallel mode:

```bash
# Sequential version (system decides which cores to use)
./build/test_louvain <graph_file> -S

# Sequential version on P-cores (performance cores)
./build/test_louvain <graph_file> -S -p

# Sequential version on E-cores (efficiency cores)
./build/test_louvain <graph_file> -S -e

# Naive parallel version with 4 threads
./build/test_louvain <graph_file> -P -n 4

# VFC parallel version with 8 threads
./build/test_louvain <graph_file> -V -n 8
```

### Command-line options:
- `-S`: Run sequential algorithm (default)
- `-P`: Run naive parallel algorithm with graph partitioning
- `-V`: Run optimized parallel algorithm with Vertex Following and Coloring
- `-n <num_threads>`: Number of threads to use (default: 1)
- `-p`: Run on P-cores (performance cores)
- `-e`: Run on E-cores (efficiency cores)

## Benchmarking with checker.py

The repository includes a benchmarking script `checker.py` that automatically runs the sequential algorithm on both P-cores and E-cores as baselines, and then tests your specified algorithms with different thread counts.

### Basic Usage

```bash
python src/checker.py inputs/community_graph_5e5.txt --algorithm naive
```

This will:
1. Run sequential algorithm on both P-cores and E-cores as baselines
2. Run the naive parallel algorithm with default thread counts
3. Report performance metrics including speedups relative to both baselines

### All Available Options

```bash
python src/checker.py <input_file> [--algorithm ALGORITHM] [--threads THREADS] [--runs RUNS] [--executable EXECUTABLE_PATH]
```

### Parameters

- `<input_file>`: Required. Path to the graph file. Default test file is `inputs/community_graph_5e5.txt`
- `--algorithm`: Comma-separated list of algorithms to test. Options: `naive,vfc`. Default: `naive,vfc`
- `--threads`: Comma-separated list of thread counts to test. Default: `1,2,4,8`
- `--runs`: Number of runs for each configuration for more stable results. Default: `1`
- `--executable`: Path to the test_louvain executable. Default: auto-detected in build directory

### Examples

```bash
# Run the naive parallel algorithm with default settings
python src/checker.py inputs/community_graph_5e5.txt --algorithm naive

# Test VFC algorithm with different thread counts
python src/checker.py inputs/community_graph_5e5.txt --algorithm vfc --threads 2,4,8,16

# Test both algorithms with custom thread counts and multiple runs
python src/checker.py inputs/community_graph_5e5.txt --algorithm naive,vfc --threads 1,2,4,8,16 --runs 3

# Specify a custom executable location
python src/checker.py inputs/community_graph_5e5.txt --algorithm naive --executable ./custom/path/test_louvain.exe
```

## File Formats

The input graph file should be in the following format:
```
<num_nodes> <num_edges>
<node1> <node2> [weight]
<node3> <node4> [weight]
...
```

- First line: number of nodes and edges
- Each subsequent line: a pair of nodes representing an edge, with an optional weight
- If weight is omitted, it defaults to 1.0
- The graph is treated as undirected