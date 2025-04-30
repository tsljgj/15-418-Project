# 15-418 Project: Parallel Louvain Community Detection

This project implements parallel versions of the Louvain community detection algorithm for graph analysis. Louvain is a widely used algorithm for discovering communities in large networks by optimizing modularity. The project was completed for **CMU 15-418: Parallel Computer Architecture and Programming, Spring 2025**.

**Note**: Please check out the `windows` branch of this repository if you want to run on a heterogeneous system (specifically Windows with an i9-14900K with 8 P cores and 16 E cores).

## Project Documentation

- [Final Report](https://drive.google.com/file/d/1VzAY6lRszolpZRC6qHhj5d_nXmYvDBSe/view?usp=sharing)
- [Project Slides](https://drive.google.com/file/d/1h51lXF6JB_K5gC3VmE7N-P_WmjGbdS3B/view?usp=sharing)

## System Requirements

- Linux operating system (homogeneous system)
- C++ compiler with C++11 support
- CMake 3.10 or higher
- OpenMP support
- Python 3 with colorama (for running benchmarks)

## Building the Project

Build the project using the provided script:

```bash
chmod +x build.sh
./build.sh
```

The executable will be created in the `build` directory.

## Algorithm Implementations

This repository contains three implementations of the Louvain community detection algorithm:

1. **Sequential (-S)**: The baseline sequential implementation
2. **Naive Parallel (-P)**: A parallelized version using simple graph partitioning
3. **VFC Parallel (-V)**: An optimized parallel version using Vertex Following and Coloring techniques

### Key Differences Between Implementations

- **Sequential**: Standard implementation that processes nodes one at a time
- **Naive Parallel (-P)**: 
  - Partitions the graph into roughly equal parts
  - Each thread processes its assigned partition
  - Uses atomic operations to handle shared data structures

- **VFC Parallel (-V)**:
  - Uses graph coloring to create conflict-free "isolate sets" of vertices
  - Implements vertex following to preprocess the graph (merging single-degree vertices with neighbors)
  - Minimizes thread contention for better scalability

## Running the Algorithm

### Default Run (Sequential Version)

To run the algorithm with default settings (sequential mode):

```bash
./build/test_louvain <graph_file>
```

### With Command-line Options

```bash
# Sequential version (explicitly specified)
./build/test_louvain <graph_file> -S

# Naive parallel version with 4 threads
./build/test_louvain <graph_file> -P -n 4

# VFC parallel version with 8 threads
./build/test_louvain <graph_file> -V -n 8
```

### Command-line Options:
- `-S`: Run sequential algorithm (default if no option specified)
- `-P`: Run naive parallel algorithm with graph partitioning
- `-V`: Run optimized parallel algorithm with Vertex Following and Coloring
- `-n <num_threads>`: Number of threads to use (default: 1)

## Running Benchmarks

The repository includes a benchmarking script to compare the performance of different implementations:

```bash
# Basic benchmark with default settings
python src/checker.py inputs/community_graph.txt

# Custom benchmark with specific algorithms, thread counts, and runs
python src/checker.py inputs/community_graph.txt --algorithm sequential,naive,vfc --threads 1,2,4,8 --runs 3
```

### Benchmark Options:
- `--algorithm`: Comma-separated list of algorithms to test (sequential,naive,vfc)
- `--threads`: Comma-separated list of thread counts to test
- `--runs`: Number of runs for each configuration
- `--executable`: Path to the test_louvain executable

## Input Graph Format

Input graphs should be in the following format:
- First line: `<number_of_nodes> <number_of_edges>`
- Following lines: `<source_node> <target_node> [weight]`
  
If the weight is omitted, it defaults to 1.0. All graphs are treated as undirected.

## Project Structure

```
project_root/
├── src/                        # Algorithm source files
│   ├── louvain_seq.h/cpp       # Sequential implementation
│   ├── louvain_parallel.h/cpp  # Naive parallel implementation
│   ├── louvain_parallel_vfc.h/cpp # VFC parallel implementation
│   ├── graph.h/cpp             # Graph data structures
│   ├── hierarchy.h             # Hierarchy data structure
│   ├── test_louvain.cpp        # Main executable
│   └── checker.py              # Benchmarking script
├── build/                      # Build output (created by build script)
├── inputs/                     # Input graph files
└── build.sh                    # Build script
```