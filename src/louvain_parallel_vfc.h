#ifndef LOUVAIN_PARALLEL_VFC_H
#define LOUVAIN_PARALLEL_VFC_H

#include "graph.h"
#include "hierarchy.h"

// Run parallel Louvain with Vertex-Following + conflict-free isolate sets.
// Populates hierarchy.partitions just like the other variants.
void louvainParallelVFC(const Graph &g0, Hierarchy &hierarchy, int numThreads);

#endif // LOUVAIN_PARALLEL_VFC_H
