// louvain_parallel_static_bl.h 
#ifndef LOUVAIN_PARALLEL_STATIC_BL_H
#define LOUVAIN_PARALLEL_STATIC_BL_H

#include "graph.h"
#include "hierarchy.h"

// The parallel Louvain baseline entry point:
void louvainParallelStaticBL(const Graph &g, Hierarchy &H, int numThreads, int pCoreCount = 0, int eCoreCount = 0);

#endif // LOUVAIN_PARALLEL_STATIC_BL_H