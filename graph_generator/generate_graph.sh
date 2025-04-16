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
./graph_generator community_graph_5e5.txt community 50000 50 0.3 0.02 1.0 10.0

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
