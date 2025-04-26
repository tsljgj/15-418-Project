#include "louvain_seq.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <thread>
#include <windows.h>

// Thread affinity settings for Windows
void setCoreAffinity(CoreType coreType) {
    // Get number of available CPUs
    int numCPUs = std::thread::hardware_concurrency();
    
    // Print the current thread ID
    DWORD threadId = GetCurrentThreadId();
    std::cout << "Current Thread ID: " << threadId << std::endl;
    
    // Create affinity mask
    DWORD_PTR affinityMask = 0;
    
    if (coreType == ANY_CORE) {
        // Use all available cores - let the system decide
        for (int i = 0; i < numCPUs; i++) {
            affinityMask |= (static_cast<DWORD_PTR>(1) << i);
        }
        std::cout << "No specific core affinity set - using any available core (system decides)\n";
        std::cout << "Available cores: ";
        for (int i = 0; i < numCPUs; i++) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
        
        // For the "any core" option, we can either:
        // 1. Set the affinity mask to include all cores (which we're doing)
        // 2. Not set an affinity mask at all (which would have the same effect)
        
        // Option 2 would be:
        // return;  // Don't set any affinity, let system decide
        
    } else if (numCPUs >= 8) {  // Assuming i9-14900K with P-cores and E-cores
        if (coreType == P_CORE) {
            // P-cores are typically the first half on hybrid architectures
            for (int i = 0; i < numCPUs / 2; i++) {
                affinityMask |= (static_cast<DWORD_PTR>(1) << i);
            }
            std::cout << "Setting affinity to P-cores for Thread ID: " << threadId << std::endl;
            
            // Print which cores we're targeting
            std::cout << "Target cores: ";
            for (int i = 0; i < numCPUs / 2; i++) {
                std::cout << i << " ";
            }
            std::cout << std::endl;
        } else if (coreType == E_CORE) {
            // E-cores are typically the second half
            for (int i = numCPUs / 2; i < numCPUs; i++) {
                affinityMask |= (static_cast<DWORD_PTR>(1) << i);
            }
            std::cout << "Setting affinity to E-cores for Thread ID: " << threadId << std::endl;
            
            // Print which cores we're targeting
            std::cout << "Target cores: ";
            for (int i = numCPUs / 2; i < numCPUs; i++) {
                std::cout << i << " ";
            }
            std::cout << std::endl;
        }
    } else {
        // If fewer cores or cannot determine, use all available
        for (int i = 0; i < numCPUs; i++) {
            affinityMask |= (static_cast<DWORD_PTR>(1) << i);
        }
        std::cout << "Warning: Cannot determine P/E cores. Using all available cores for Thread ID: " << threadId << std::endl;
    }
    
    // Set the CPU affinity for the current thread
    HANDLE currentThread = GetCurrentThread();
    BOOL result = SetThreadAffinityMask(currentThread, affinityMask);
    if (result == 0) {
        std::cerr << "Error setting thread affinity: " << GetLastError() << " for Thread ID: " << threadId << std::endl;
    } else {
        std::cout << "Affinity set successfully for Thread ID: " << threadId << std::endl;
        
        // Get and print the actual affinity mask that was set
        DWORD_PTR processAffinityMask, systemAffinityMask;
        if (GetProcessAffinityMask(GetCurrentProcess(), &processAffinityMask, &systemAffinityMask)) {
            std::cout << "Process affinity mask: 0x" << std::hex << processAffinityMask << std::dec << std::endl;
        }
    }
}

// -------------------------------------------------------------------------
// Graph Reading function (supports weighted and unweighted graphs).
Graph readGraph(const std::string &filename) {
    Graph g;
    std::ifstream in(filename);
    if (!in.is_open()) {
        std::cerr << "Error opening input file: " << filename << "\n";
        exit(EXIT_FAILURE);
    }
    std::string line;
    // First line: number of nodes and (approximate) number of edges.
    if (std::getline(in, line)) {
        std::istringstream iss(line);
        iss >> g.n;
        // The edge count from the header is not strictly used.
    } else {
        std::cerr << "Error reading first line of file.\n";
        exit(EXIT_FAILURE);
    }
    // Initialize adjacency list and degree vector.
    g.adj.resize(g.n);
    g.degree.assign(g.n, 0.0);

    int u, v;
    double w;
    int edgeCount = 0;
    while (std::getline(in, line)) {
        if (line.empty())
            continue;
        std::istringstream iss(line);
        // Read u, v; if weight is not provided, default to 1.0.
        if (!(iss >> u >> v))
            continue;
        if (!(iss >> w))
            w = 1.0;
        if (u < 0 || u >= g.n || v < 0 || v >= g.n)
            continue; // skip invalid nodes.
        // Undirected edge: add both directions.
        g.adj[u].push_back({v, w});
        g.adj[v].push_back({u, w});
        g.degree[u] += w;
        g.degree[v] += w;
        edgeCount++;
    }
    in.close();

    double totalDegree = 0.0;
    for (int i = 0; i < g.n; ++i)
        totalDegree += g.degree[i];
    g.m = totalDegree / 2.0; // Each edge counted twice.
    std::cout << "Graph loaded: " << g.n << " nodes, " << edgeCount 
              << " edge lines read, total weighted m = " << g.m << "\n";
    return g;
}

// -------------------------------------------------------------------------
// Computes modularity Q using the weighted formula.
// Q = (1/(2*m)) * sum_c [ L_c - (d_c^2)/(2*m) ]
double computeModularity(const Graph &g, const std::vector<int> &community) {
    int n = g.n;
    double m_total = 2.0 * g.m; // Total weight (each edge appears twice)
    std::unordered_map<int, double> comm_tot;  // Total degree per community.
    std::unordered_map<int, double> comm_in;   // Sum of internal edge weights.

    for (int i = 0; i < n; ++i) {
        int c = community[i];
        comm_tot[c] += g.degree[i];
        for (const auto &edge : g.adj[i]) {
            if (edge.neighbor > i && community[edge.neighbor] == c)
                comm_in[c] += edge.weight;
        }
    }
    double Q = 0.0;
    for (const auto &entry : comm_tot) {
        int c = entry.first;
        double tot = entry.second;
        double in = comm_in[c];
        Q += (in / g.m) - std::pow(tot / m_total, 2);
    }
    return Q;
}

// -------------------------------------------------------------------------
// Helper: Performs one pass (sweep) of local node moves.
// For each node, computes the gain in modularity if moving it to a neighboring community
// and reassigns the node if the best gain is positive.
// It updates the "community" vector and a helper vector "comm_tot" which holds the total weighted degree for each community.
// Returns true if any node was moved.
static bool localMovePass(const Graph &g, std::vector<int> &community, std::vector<double> &comm_tot) {
    bool moved = false;
    int n = g.n;
    double m_total = 2.0 * g.m;
    for (int i = 0; i < n; ++i) {
        int current_comm = community[i];
        double k_i = g.degree[i];
        // Compute weight sum from node i to each neighbor community.
        std::unordered_map<int, double> neighCommWeight;
        for (const auto &edge : g.adj[i]) {
            int nb_comm = community[edge.neighbor];
            neighCommWeight[nb_comm] += edge.weight;
        }
        // Temporarily remove node i from its community.
        comm_tot[current_comm] -= k_i;
        community[i] = -1;
        double best_gain = 0.0;
        int best_comm = current_comm;
        // Evaluate gain for each neighboring community.
        for (const auto &entry : neighCommWeight) {
            int c = entry.first;
            double k_i_in = entry.second;
            double gain = k_i_in - (comm_tot[c] * k_i) / m_total;
            if (gain > best_gain) {
                best_gain = gain;
                best_comm = c;
            }
        }
        if (best_comm == -1)
            best_comm = current_comm;
        community[i] = best_comm;
        comm_tot[best_comm] += k_i;
        if (best_comm != current_comm)
            moved = true;
    }
    return moved;
}

// -------------------------------------------------------------------------
// Local move phase: repeatedly perform local node moves until convergence.
// Returns the partition for the current graph (each node's community).
static std::vector<int> localMovePhase(const Graph &g) {
    int n = g.n;
    std::vector<int> partition(n);
    for (int i = 0; i < n; ++i)
        partition[i] = i; // initially, each node in its own community.
    std::vector<double> comm_tot(n, 0.0);
    for (int i = 0; i < n; ++i)
        comm_tot[partition[i]] = g.degree[i];
    
    double currentQ = computeModularity(g, partition);
    int iter = 0;
    bool improvement = true;
    while (improvement) {
        iter++;
        bool moved = localMovePass(g, partition, comm_tot);
        double newQ = computeModularity(g, partition);
        std::cout << "Local move pass " << iter << ": modularity = " << newQ << "\n";
        if (std::fabs(newQ - currentQ) < 1e-6)
            break;
        currentQ = newQ;
        improvement = moved;
    }
    return partition;
}

// -------------------------------------------------------------------------
// Aggregates graph g given a partition.
// Creates a new graph where each node corresponds to a community from "partition".
// The mapping "commMap" (old community label -> new node index) is filled by this function.
Graph aggregateGraph(const Graph &g, const std::vector<int> &partition, std::unordered_map<int,int> &commMap) {
    Graph g_new;
    commMap.clear();
    // Create mapping: assign consecutive indices for unique communities.
    for (int i = 0; i < g.n; ++i) {
        int comm = partition[i];
        if (commMap.find(comm) == commMap.end()) {
            commMap[comm] = commMap.size();
        }
    }
    g_new.n = commMap.size();
    g_new.adj.resize(g_new.n);
    g_new.degree.assign(g_new.n, 0.0);

    // Build aggregated edge weights: for each edge in g, add its weight between the corresponding communities.
    std::vector<std::unordered_map<int, double>> aggEdges(g_new.n);
    for (int i = 0; i < g.n; ++i) {
        int ci = commMap[ partition[i] ];
        for (const auto &edge : g.adj[i]) {
            int cj = commMap[ partition[edge.neighbor] ];
            aggEdges[ci][cj] += edge.weight;
        }
    }
    g_new.m = 0.0;
    for (int i = 0; i < g_new.n; ++i) {
        for (const auto &entry : aggEdges[i]) {
            int j = entry.first;
            double w = entry.second;
            g_new.adj[i].push_back({j, w});
            g_new.degree[i] += w;
            if (i <= j)
                g_new.m += w;  // Count edge once.
        }
    }
    return g_new;
}

// -------------------------------------------------------------------------
// Recursive function that computes the full hierarchical decomposition.
// For each level, the partition is recorded in the Hierarchy object.
// The base case is reached when no aggregation is possible (i.e. each node remains separate).
static void runLouvainHierarchy(const Graph &g, Hierarchy &H) {
    // Level partition for current graph.
    std::vector<int> partition = localMovePhase(g);
    double Q = computeModularity(g, partition);
    std::cout << "Level modularity: " << Q << ", ";
    // Count unique communities.
    std::unordered_set<int> uniqueComm(partition.begin(), partition.end());
    std::cout << "communities = " << uniqueComm.size() << "\n";
    H.partitions.push_back(partition);
    
    // If every node is alone (no aggregation possible), stop.
    if (uniqueComm.size() == (size_t)g.n)
        return;
    
    // Aggregate graph.
    std::unordered_map<int,int> commMap;
    Graph g_new = aggregateGraph(g, partition, commMap);
    std::cout << "Aggregated graph: " << g.n << " -> " << g_new.n << " nodes.\n";
    
    // Recursive call on the aggregated graph.
    runLouvainHierarchy(g_new, H);
}

// -------------------------------------------------------------------------
// Top-level interface for the hierarchical sequential Louvain algorithm.
// This function computes the full hierarchical decomposition and stores it in hierarchy.
void louvainHierarchical(const Graph &g, Hierarchy &hierarchy, CoreType coreType) {
    // Set core affinity before starting the algorithm
    setCoreAffinity(coreType);
    
    hierarchy.partitions.clear();
    runLouvainHierarchy(g, hierarchy);
}