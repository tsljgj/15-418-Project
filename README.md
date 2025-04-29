# Louvain Community Detection on Heterogeneous Systems

This repository implements sequential and parallel versions of the Louvain community detection algorithm for large graphs, optimized for heterogeneous processors (i9-14900K with 8 P-cores and 16 E-cores).

## System Requirements

- **Operating System**: Windows
- **Compiler**: MinGW (GNU Compiler Collection for Windows)
- **Build System**: CMake 3.10 or higher
- **Processor**: Intel hybrid architecture (tested on i9-14900K)
- **Python Dependencies** (for benchmarking scripts):
  - Python 3.x
  - colorama package (`pip install colorama`)
  - pandas (optional)

## File Structure

```
louvain-community-detection/
├── README.md              # Project overview and usage instructions
├── build.bat              # Build script for Windows with MinGW
├── build/                 # Created after compilation
│   └── test_louvain.exe   # Main executable
├── src/                   # Source code
│   ├── CMakeLists.txt     # CMake configuration
│   ├── checker.py         # Benchmarking script
│   ├── core_type.cpp      # Heterogeneous core management
│   ├── core_type.h
│   ├── graph.cpp          # Graph data structure
│   ├── graph.h
│   ├── hierarchy.h        # Community hierarchy structure
│   ├── louvain_seq.cpp    # Sequential algorithm
│   ├── louvain_seq.h
│   ├── louvain_parallel.cpp          # Naive parallel implementation
│   ├── louvain_parallel.h
│   ├── louvain_parallel_bl.cpp       # Big.LITTLE-aware implementation
│   ├── louvain_parallel_bl.h
│   ├── louvain_parallel_static.cpp   # Static scheduling implementation
│   ├── louvain_parallel_static.h
│   ├── louvain_parallel_static_bl.cpp # Static Big.LITTLE implementation
│   ├── louvain_parallel_static_bl.h
│   ├── louvain_parallel_vfc.cpp      # VFC implementation
│   ├── louvain_parallel_vfc.h
│   ├── test_all.bat       # Batch testing script
│   ├── test_louvain.cpp   # Main executable source
│   ├── test_naive_bl.bat  # Big.LITTLE testing script
│   └── test_static.py     # Advanced heterogeneous testing
└── inputs/                # Input graph datasets
```

## Building the Project

This project uses MinGW for compilation on Windows:

```bash
build.bat
```

This creates the executable at `build/test_louvain.exe`.

## Available Algorithms

- **Sequential**: Baseline implementation with options for P-cores or E-cores
- **Naive Parallel**: Simple parallel implementation using graph partitioning
- **VFC Parallel**: Optimized parallel algorithm with Vertex Following and Coloring
- **Big.LITTLE Parallel**: Heterogeneous-aware implementation for mixed core types
- **Static Scheduling**: Parallel version with workload-balancing
- **Static Big.LITTLE**: Combined static scheduling with heterogeneous awareness

## Running the Algorithm

### Quick Start (Default Run)

```bash
# Sequential version (system decides which cores to use)
build\test_louvain.exe inputs\community_1e4_ultrasparse.txt
```

### Algorithm Options

```bash
# Sequential version on P-cores
build\test_louvain.exe inputs\community_1e4_ultrasparse.txt -S -p

# Sequential version on E-cores
build\test_louvain.exe inputs\community_1e4_ultrasparse.txt -S -e

# Naive parallel version with 4 threads (system decides cores)
build\test_louvain.exe inputs\community_1e4_ultrasparse.txt -P -a 4

# VFC parallel version with 8 threads (system decides cores)
build\test_louvain.exe inputs\community_1e4_ultrasparse.txt -V -a 8

# Parallel using specific core counts (2 P-cores, 4 E-cores)
build\test_louvain.exe inputs\community_1e4_ultrasparse.txt -P -pc 2 -ec 4

# Big.LITTLE aware version with specific core counts
build\test_louvain.exe inputs\community_1e4_ultrasparse.txt -B -pc 2 -ec 6
```

## Benchmarking with checker.py

The repository includes a benchmarking script that automatically runs tests with different configurations:

```bash
# Requirements: Python 3.x with colorama package
pip install colorama

# Basic usage (tests naive algorithm with default settings)
python src/checker.py inputs\community_1e4_ultrasparse.txt --algorithm naive

# Test multiple algorithms with various thread counts
python src/checker.py inputs\community_1e4_ultrasparse.txt --algorithm naive,vfc --threads 1,2,4,8,16
```

### Command-line Options for checker.py

- `<input_file>`: Graph file path
- `--algorithm`: Comma-separated list (options: naive,vfc,naive_bl,static,static_bl)
- `--threads`: Thread counts to test
- `--p-e-ratio`: P:E core ratios to test (e.g., 4:12,8:8,2:16)
- `--runs`: Number of runs per configuration (default: 1)

## Input Graph Format

```
<num_nodes> <num_edges>
<node1> <node2> [weight]
<node3> <node4> [weight]
...
```

Weight defaults to 1.0 if omitted, and the graph is undirected.