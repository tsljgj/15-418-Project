#ifndef LOUVAIN_PARALLEL_VFC_H
#define LOUVAIN_PARALLEL_VFC_H

#include "louvain_seq.h"
#include "graph.h"
#include "hierarchy.h"
#include "core_type.h"

#include <omp.h>
#include <vector>
#include <unordered_map>
#include <mutex>

// Create isolate sets for conflict-free parallel processing
// Isolate sets contain vertices that don't share neighbors
std::vector<std::vector<int>> createIsolateSets(const Graph &g);

// Apply the Vertex Following heuristic (merge single-degree vertices with their neighbors)
// Returns the compressed graph and updates the mapping from original to compressed vertices
Graph applyVertexFollowing(const Graph &g, std::vector<int> &mapping);

// Perform a local move phase using isolate sets for conflict-free parallelism
bool localMovePhaseWithIsolateSets(
    Graph &g,
    std::vector<int> &communities, 
    const std::vector<std::vector<int>> &isolateSets,
    double threshold
);

// Create an aggregated graph based on communities with parallel processing
Graph aggregateGraphFast(
    const Graph &g,
    const std::vector<int> &communities, 
    std::unordered_map<int, int> &commMap
);

// Run the Louvain algorithm with Vertex Following and Coloring heuristics,
// optionally pinning threads to a specified number of P-cores and E-cores
void louvainParallelVFC(
    const Graph &g,
    Hierarchy &hierarchy,
    int numThreads,
    int pCoreCount = 0,
    int eCoreCount = 0
);

#endif // LOUVAIN_PARALLEL_VFC_H
