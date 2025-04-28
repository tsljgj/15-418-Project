// louvain_parallel_bl.h
#ifndef LOUVAIN_PARALLEL_STATIC_BL_H
#define LOUVAIN_PARALLEL_STATIC_BL_H

#include "graph.h"
#include "hierarchy.h"

// The parallel Louvain baseline entry point:
void louvainParallelStaticBL(const Graph &g, Hierarchy &H, int numThreads);

#endif // LOUVAIN_PARALLEL_STATIC_BL_H