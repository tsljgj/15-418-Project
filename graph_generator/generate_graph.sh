#!/bin/bash

# Simple script to generate one of each type of graph
# All graphs are saved in the current directory with default parameters

echo "Generating test graphs..."
echo "-------------------------------------"

# Generate random graph (500 nodes, 2500 edges)
# echo "Generating random graph..."
# ./graph_generator random_graph.txt random 500 2500 1.0 10.0

# Generate community graph
echo "Generating community graph..."
# ./graph_generator small-world_5e5_dense_intra.txt community 50000 50 0.9 0.01 1.0 15.0
# ./graph_generator small-world_2k_v_dense.txt community 2000 10 0.99 0.1 1.0 15.0
./graph_generator random_5e5_dense.txt community 50000 5 0.95 0.2 1.0 15.0
./graph_generator small-world_1e5_med_com.txt community 10000 200 0.6 0.05 1.0 15.0
# ./graph_generator preferential_2k_dense_overlap.txt community 2000 50 0.8 0.3 1.0 15.0
./graph_generator preferential_1e5_med_sparse.txt community 10000 150 0.5 0.02 1.0 15.0
# ./graph_generator preferential_1k_hi_intra_inter.txt community 1000 20 0.7 0.6 1.0 15.0
./graph_generator Hierarchical_1e5_few_dense.txt community 10000 30 0.9 0.005 1.0 15.0
# ./graph_generator Hierarchical_2k_many_sparse.txt community 2000 200 0.2 0.05 1.0 15.0
# ./graph_generator Hierarchical_5e5_med_overlap.txt community 500000 120 0.55 0.1 1.0 15.0

# Generate preferential attachment graph (500 nodes)
# echo "Generating preferential attachment graph..."
# ./graph_generator preferential_graph.txt preferential 500 10 5 1.0 10.0

# Generate small-world graph (500 nodes)
# echo "Generating small-world graph..."
# ./graph_generator smallworld_graph.txt smallworld 500 8 0.1 1.0 10.0

echo "-------------------------------------"
# echo "Generated 4 graph types in the current directory:"
# echo "  - random_graph.txt"
echo "  - community_graph_5e5.txt"
# echo "  - preferential_graph.txt"
# echo "  - smallworld_graph.txt"
