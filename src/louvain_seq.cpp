// louvain_seq.cpp (updated)
#include "louvain_seq.h"
#include "graph.h"
#include "hierarchy.h"
#include "core_type.h"

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

// Top-level interface for the hierarchical sequential Louvain algorithm.
// This function computes the full hierarchical decomposition and stores it in hierarchy.
void louvainHierarchical(const Graph &g, Hierarchy &hierarchy, CoreType coreType) {
    // Set core affinity before starting the algorithm
    setThreadAffinityCoreType(coreType);
    
    hierarchy.partitions.clear();
    runLouvainHierarchy(g, hierarchy);
}