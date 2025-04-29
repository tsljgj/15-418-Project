#ifndef LOUVAIN_PARALLEL_VFC_BL_H
#define LOUVAIN_PARALLEL_VFC_BL_H

#include "graph.h"
#include "hierarchy.h"

void louvainParallelVFCBL(const Graph &g, Hierarchy &H, int numThreads, int PCore, int Ecore);

#endif // LOUVAIN_PARALLEL_VFC_BL_H