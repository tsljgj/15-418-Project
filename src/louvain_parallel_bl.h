// louvain_parallel_bl.h
#ifndef LOUVAIN_PARALLEL_BL_H
#define LOUVAIN_PARALLEL_BL_H

#include "graph.h"
#include "hierarchy.h"

void louvainParallelBL(const Graph &g, Hierarchy &H, int numThreads, 
                      int pCoreCount = 0, int eCoreCount = 0);

#endif // LOUVAIN_PARALLEL_BL_H