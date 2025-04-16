#ifndef LOUVAIN_PARALLEL_H
#define LOUVAIN_PARALLEL_H

#include "louvain_seq.h"
#include <omp.h>
#include <vector>
#include <unordered_map>

// Structure to hold graph partitions for parallel processing
struct GraphPartition {
    std::vector<int> nodeIds;                 // Original node IDs in this partition
    std::vector<std::vector<Edge>> localAdj;  // Local adjacency list
    std::vector<double> localDegree;          // Local degree vector
    std::unordered_map<int, int> globalToLocal; // Maps global node ID to local index
    std::unordered_map<int, int> localToGlobal; // Maps local index to global node ID
};

// Function to partition a graph into numPartitions parts
std::vector<GraphPartition> partitionGraph(const Graph &g, int numPartitions);

// Aggregates graph g given a partition.
// Creates a new graph where each node corresponds to a community from "partition".
Graph aggregateGraph(const Graph &g, const std::vector<int> &partition, std::unordered_map<int,int> &commMap);

// Runs local move phase on a partition and returns updated community assignments
std::vector<int> parallelLocalMovePhase(const Graph &g, const GraphPartition &partition, 
                                       const std::vector<int> &initialCommunities);

// Parallel implementation of the hierarchical Louvain algorithm
void louvainParallel(const Graph &g, Hierarchy &hierarchy, int numThreads);

#endif // LOUVAIN_PARALLEL_H