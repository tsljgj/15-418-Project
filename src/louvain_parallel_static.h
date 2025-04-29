#ifndef LOUVAIN_PARALLEL_STATIC_H
#define LOUVAIN_PARALLEL_STATIC_H

#include "graph.h"
#include "hierarchy.h"

void louvainParallelStatic(const Graph &g, Hierarchy &H, int numThreads, int pCoreCount = 0, int eCoreCount = 0);

#endif // LOUVAIN_PARALLEL_STATIC_H