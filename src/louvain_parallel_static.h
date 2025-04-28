#ifndef LOUVAIN_PARALLEL_STATIC_H
#define LOUVAIN_PARALLEL_STATIC_H

#include "graph.h"
#include "hierarchy.h"

void louvainParallelStatic(const Graph &g, Hierarchy &H, int numThreads);

#endif // LOUVAIN_PARALLEL_STATIC_H