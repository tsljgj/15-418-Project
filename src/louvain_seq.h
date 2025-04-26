#ifndef LOUVAIN_SEQ_H
#define LOUVAIN_SEQ_H

#include <vector>
#include <string>
#include "graph.h"      // brings in Edge, Graph, readGraph, computeModularity
#include "hierarchy.h"  // brings in Hierarchy

// Enum for specifying CPU core type
enum CoreType {
    ANY_CORE,  // No specific core preference
    P_CORE,    // Performance core
    E_CORE     // Efficiency core
};

// Helper function to set thread affinity to specific core type
void setCoreAffinity(CoreType coreType);

// Runs the sequential Louvain algorithm (local move phase plus aggregation)
// and builds a full hierarchical decomposition stored in "hierarchy".
// The full hierarchy is stored as a sequence of partitions.
// Added coreType parameter to specify which type of core to run on.
void louvainHierarchical(const Graph &g, Hierarchy &hierarchy, CoreType coreType = P_CORE);

#endif // LOUVAIN_SEQ_H
