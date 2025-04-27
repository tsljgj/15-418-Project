#ifndef LOUVAIN_PARALLEL_H
#define LOUVAIN_PARALLEL_H

#include "graph.h"
#include <vector>
#include <unordered_map>
#include "hierarchy.h"

// Partition descriptor (though we won't copy adjacency in this simple version)
struct GraphPartition {
    std::vector<int> nodeIds;
};

// The parallel Louvain entry point:
void louvainParallel(const Graph &g, Hierarchy &H, int numThreads, int pCoreCount = 0, int eCoreCount = 0);

#endif // LOUVAIN_PARALLEL_H