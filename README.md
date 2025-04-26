# Running the Algorithm

You can run the algorithm in sequential or parallel mode:

```bash
# Sequential version on default P-cores
./build/test_louvain <graph_file> -S

# Sequential version on E-cores
./build/test_louvain <graph_file> -S -e

# Sequential version on P-cores (explicitly specified)
./build/test_louvain <graph_file> -S -p

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
- `-p`: Run on P-cores (performance cores, default for sequential algorithm)
- `-e`: Run on E-cores (efficiency cores, optional for sequential algorithm)