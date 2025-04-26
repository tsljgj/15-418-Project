# 15-418 Project: Parallel Louvain Community Detection

Welcome to our 15-418 Project repository! This project implements parallel versions of the Louvain community detection algorithm for graph analysis.

## Project Proposal

For a detailed overview and further information, please review our [Project Proposal](https://docs.google.com/document/d/1WRLmLnNcaVsU4yoNn0Z73KKaIhtWlAypexBPJcd5nok/edit?usp=sharing).

## Milestone Report

For a detailed overview and further information, please review our [Milestone Report](https://docs.google.com/document/d/1c3C9YgUiPjCyHc-xrMc9cLpcsJeDmjJKt6cI9FQebZM/edit?usp=sharing).

## Algorithm Implementations

This repository contains three implementations of the Louvain community detection algorithm:

1. **Sequential**: The baseline sequential implementation
2. **Naive Parallel**: A parallelized version using simple graph partitioning
3. **VFC Parallel**: An optimized parallel version using Vertex Following and Coloring techniques

### Key Differences Between Parallel Implementations

- **Naive Parallel (-P)**: 
  - Partitions the graph into roughly equal parts
  - Each thread processes its assigned partition
  - Uses atomic operations to handle shared data structures
  - Simple but can suffer from load imbalance and high synchronization overhead

- **VFC Parallel (-V)**:
  - Uses graph coloring to create conflict-free "isolate sets" of vertices
  - Implements vertex following to preprocess the graph (merging single-degree vertices with neighbors)
  - Minimizes thread contention by processing independent vertex sets in parallel
  - More sophisticated synchronization strategy for better scalability

## Building the Project

### Prerequisites

- CMake 3.10 or higher
- C++ compiler with C++11 support
- OpenMP support
- Python 3 with colorama (for running benchmarks)

### Building

1. Build the project using the provided script:
   ```bash
   chmod +x build.sh
   ./build.sh
   ```

2. The executable will be created in the `build` directory.

## Running the Algorithm

You can run the algorithm in sequential or parallel mode:

```bash
# Sequential version
./build/test_louvain <graph_file> -S

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

## Benchmarking

The repository includes a benchmarking script to compare the performance of different algorithm implementations:

```bash
# Basic benchmark with default settings
python src/checker.py inputs/community_graph.txt

# Custom benchmark specifying algorithms, thread counts, and runs
python src/checker.py inputs/community_graph.txt --algorithm sequential,naive,vfc --threads 1,2,4,8 --runs 3
```

### Benchmark options:
- `--algorithm`: Comma-separated list of algorithms to test (sequential,naive,vfc)
- `--threads`: Comma-separated list of thread counts to test
- `--runs`: Number of runs for each configuration
- `--executable`: Path to the test_louvain executable

## Graph Generator

The repository includes tools for generating test graphs in the `graph_generator` directory:

```bash
# Generate a random graph
./graph_generator/graph_generator random_graph.txt random 1000 5000

# Generate a community graph (good for testing community detection)
./graph_generator/graph_generator community_graph.txt community 1000 10 0.3 0.02

# Generate a preferential attachment graph
./graph_generator/graph_generator preferential_graph.txt preferential 1000 10 3

# Generate a small-world graph
./graph_generator/graph_generator smallworld_graph.txt smallworld 1000 6 0.1
```

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
│   ├── test_louvain.cpp        # Main executable
│   └── checker.py              # Benchmarking script
├── build/                      # Build output (created by build script)
├── graph_generator/            # Graph generator tools
├── inputs/                     # Input graph files
└── build.sh                    # Build script
```
