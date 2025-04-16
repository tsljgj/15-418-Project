# 15-418 Project

Welcome to our 15-418 Project repository!

## Project Proposal

For a detailed overview and further information, please review our [Project Proposal](https://docs.google.com/document/d/1WRLmLnNcaVsU4yoNn0Z73KKaIhtWlAypexBPJcd5nok/edit?usp=sharing).

## Milestone Report

For a detailed overview and further information, please review our [Milestone Report](https://docs.google.com/document/d/1c3C9YgUiPjCyHc-xrMc9cLpcsJeDmjJKt6cI9FQebZM/edit?usp=sharing).


## Building the Project

### Prerequisites

- CMake 3.10 or higher
- C++ compiler with C++11 support
- OpenMP support

### Building

1. Build the project using the provided script:
   ```bash
   chmod +x build.sh
   ./build.sh
   ```

2. The executable will be created in the `build` directory.

## Running the Algorithm

The algorithm can be run in either sequential or parallel mode:

```bash
# Sequential version
./build/test_louvain <graph_file> -S

# Parallel version with 4 threads
./build/test_louvain <graph_file> -P -n 4
```

Command-line options:
- `-S`: Run sequential algorithm (default)
- `-P`: Run parallel algorithm with graph partitioning
- `-n <num_threads>`: Number of threads/partitions to use (default: 1)

## Graph Generator

The repository includes tools for generating test graphs in the `graph_generator` directory. You can use these to create different types of graphs for testing:

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
├── build/                      # Build output (created by build script)
├── graph_generator/            # Graph generator tools
├── inputs/                     # Input graph files (optional)
└── build.sh                    # Build script
```
