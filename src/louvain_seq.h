#ifndef LOUVAIN_SEQ_H
#define LOUVAIN_SEQ_H

#include <vector>
#include <string>

// Enum for specifying CPU core type
enum CoreType {
    ANY_CORE,  // No specific core preference
    P_CORE,    // Performance core
    E_CORE     // Efficiency core
};

// Structure for a (weighted) edge.
struct Edge {
    int neighbor;    // Index of neighbor node.
    double weight;   // Weight of the edge.
};

// Graph structure: contains number of nodes (n), total weight m (where m is the sum of weights of all undirected edges),
// an adjacency list (each edge stored as an Edge), and the weighted degree for each node.
struct Graph {
    int n;                       
    double m;                    
    std::vector<std::vector<Edge>> adj;  
    std::vector<double> degree;            
};

// Structure for storing the hierarchical decomposition.
// The vector "partitions" will hold one partition per level; level 0 is the partition
// for the current graph, level 1 for the first aggregated graph, etc.
struct Hierarchy {
    std::vector<std::vector<int>> partitions;
};

// Set thread affinity to a specific CPU ID
void setThreadAffinityToCpu(int cpuId);

// Get the appropriate CPU ID range for a core type
std::vector<int> getCpuIdsForCoreType(CoreType coreType);

// Set thread affinity based on core type
void setThreadAffinityCoreType(CoreType coreType);

// Helper function to set thread affinity to specific core type (legacy function)
void setCoreAffinity(CoreType coreType);

// Reads a graph from a file. File format:
//   First line: <number_of_nodes> <number_of_edges>
//   Each subsequent line: <node> <neighbor> [weight]
// If the weight is missing, it defaults to 1.0. The graph is undirected.
Graph readGraph(const std::string &filename);

// Computes the modularity Q of a given partition on graph g (weighted version).
// Partition: community[i] gives the community id for node i.
double computeModularity(const Graph &g, const std::vector<int> &community);

// Runs the sequential Louvain algorithm (local move phase plus aggregation)
// and builds a full hierarchical decomposition stored in "hierarchy".
// The full hierarchy is stored as a sequence of partitions.
// Added coreType parameter to specify which type of core to run on.
void louvainHierarchical(const Graph &g, Hierarchy &hierarchy, CoreType coreType = ANY_CORE);

#endif // LOUVAIN_SEQ_H