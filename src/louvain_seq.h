// louvain_seq.h (updated)
#ifndef LOUVAIN_SEQ_H
#define LOUVAIN_SEQ_H

#include "graph.h"      // brings in Edge, Graph, readGraph, computeModularity
#include "hierarchy.h"  // brings in Hierarchy
#include "core_type.h"  // brings in CoreType

// Runs the sequential Louvain algorithm (local move phase plus aggregation)
// and builds a full hierarchical decomposition stored in "hierarchy".
// The full hierarchy is stored as a sequence of partitions.
// Added coreType parameter to specify which type of core to run on.
void louvainHierarchical(const Graph &g, Hierarchy &hierarchy, CoreType coreType = ANY_CORE);

#endif // LOUVAIN_SEQ_H