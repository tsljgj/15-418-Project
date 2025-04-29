#include "louvain_parallel_static.h"
#include "louvain_seq.h"  // for fallback when numThreads == 1
#include "graph.h"
#include "hierarchy.h"
#include "core_type.h"

#include <omp.h>
#include <vector>
#include <numeric>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <cmath>

// Helper structure to store node index and its degree for workload balancing
struct NodeDegree {
    int index;
    double degree;
    
    bool operator<(const NodeDegree& other) const {
        return degree > other.degree; // Sort in descending order
    }
};

// Static parallel sweep with degree-based workload distribution
static bool parallelLocalSweepStatic(const Graph &g, std::vector<int> &C) {
    int n = g.n;
    double two_m = 2.0 * g.m;
    int numThreads = omp_get_max_threads();
    
    // 1) Build degree copy and community totals
    std::vector<double> k(n), comm_tot(n, 0.0);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        k[i] = g.degree[i];
    }
    for (int i = 0; i < n; i++) {
        comm_tot[C[i]] += k[i];
    }
    
    // 2) Create nodes sorted by degree for workload balancing
    std::vector<NodeDegree> sortedNodes(n);
    for (int i = 0; i < n; i++) {
        sortedNodes[i] = {i, k[i]};
    }
    std::sort(sortedNodes.begin(), sortedNodes.end());
    
    // 3) Distribute nodes to threads to balance workload
    std::vector<std::vector<int>> threadNodes(numThreads);
    std::vector<double> threadWorkload(numThreads, 0.0);
    
    // Use greedy algorithm to assign nodes to threads
    for (const auto& node : sortedNodes) {
        // Find thread with minimum current workload
        int minThread = 0;
        double minLoad = threadWorkload[0];
        for (int t = 1; t < numThreads; t++) {
            if (threadWorkload[t] < minLoad) {
                minLoad = threadWorkload[t];
                minThread = t;
            }
        }
        
        // Assign node to this thread
        threadNodes[minThread].push_back(node.index);
        threadWorkload[minThread] += node.degree;
    }
    
    // Print thread workload distribution
    // #pragma omp single
    // {
    //     std::cout << "Thread workload distribution: ";
    //     for (int t = 0; t < numThreads; t++) {
    //         std::cout << "Thread " << t << ": " << threadWorkload[t] << " (" 
    //                   << threadNodes[t].size() << " nodes), ";
    //     }
    //     std::cout << std::endl;
    // }
    
    // 4) Compute best moves in parallel with static assignment
    std::vector<int> bestC(n);
    bool moved = false;
    
    #pragma omp parallel reduction(||:moved)
    {
        int tid = omp_get_thread_num();
        const std::vector<int>& myNodes = threadNodes[tid];
        
        for (int i : myNodes) {
            int curr = C[i];
            double ki = k[i];
            
            // Sum weights to each neighbor community
            std::unordered_map<int,double> w2c;
            for (auto &e : g.adj[i]) {
                w2c[C[e.neighbor]] += e.weight;
            }
            
            // Find best community
            double bestGain = 0.0;
            int bestComm = curr;
            for (auto &p : w2c) {
                double gain = p.second - (comm_tot[p.first] * ki) / two_m;
                if (gain > bestGain) {
                    bestGain = gain;
                    bestComm = p.first;
                }
            }
            bestC[i] = bestComm;
            if (bestComm != curr) moved = true;
        }
    }
    
    // 5) Apply all moves (sequentially)
    for (int i = 0; i < n; i++) {
        C[i] = bestC[i];
    }
    return moved;
}

// Aggregate graph by community vector (same as in louvain_parallel.cpp)
static Graph aggregateGraph(const Graph &g,
                            const std::vector<int> &C,
                            std::vector<int> &old2new) {
    int n = g.n;
    std::unordered_map<int,int> remap;
    old2new.resize(n);
    int cid = 0;
    // assign new IDs
    for (int i = 0; i < n; i++) {
        auto it = remap.find(C[i]);
        if (it == remap.end()) {
            remap[C[i]] = cid++;
        }
        old2new[i] = remap[C[i]];
    }
    // build new graph
    Graph g2;
    g2.n = cid;
    g2.adj.assign(cid, {});
    g2.degree.assign(cid, 0.0);

    std::vector<std::unordered_map<int,double>> W(cid);
    for (int i = 0; i < n; i++) {
        int ci = old2new[i];
        for (auto &e : g.adj[i]) {
            int cj = old2new[e.neighbor];
            W[ci][cj] += e.weight;
        }
    }
    double tot = 0.0;
    for (int i = 0; i < cid; i++) {
        for (auto &p : W[i]) {
            int j = p.first;
            double w = p.second;
            if (w > 0) {
                g2.adj[i].push_back({j, w});
                g2.degree[i] += w;
                if (i <= j) tot += w;
            }
        }
    }
    g2.m = tot;
    return g2;
}

// Top‐level parallel Louvain with static scheduling
void louvainParallelStatic(const Graph &g0, Hierarchy &H, int numThreads, int pCoreCount, int eCoreCount) {
    // Fast fallback if only 1 thread requested
    if (numThreads <= 1) {
        louvainHierarchical(g0, H);
        return;
    }
    
    // Create vector of core assignments
    std::vector<int> coreAssignments;
    bool useSpecificCores = (pCoreCount > 0 || eCoreCount > 0);
    
    // If specific core counts are provided, set up thread affinity
    if (useSpecificCores) {
        coreAssignments = assignParallelCores(pCoreCount, eCoreCount);
        
        if (coreAssignments.empty()) {
            // Error already reported in assignParallelCores
            return;
        }
        
        // Print the core assignments
        // std::cout << "Core assignments: ";
        // for (size_t i = 0; i < coreAssignments.size(); i++) {
        //     std::cout << coreAssignments[i];
        //     if (i < coreAssignments.size() - 1) {
        //         std::cout << ", ";
        //     }
        // }
        // std::cout << std::endl;
        
        // The total number of threads will be the sum of the core counts
        numThreads = pCoreCount + eCoreCount;
    } else {
        std::cout << "Using " << numThreads << " threads with system-decided core affinity\n";
    }
    
    omp_set_num_threads(numThreads);
    
    // Set thread affinity ONCE at the beginning
    if (useSpecificCores) {
        #pragma omp parallel
        {
            int threadId = omp_get_thread_num();
            if (threadId < static_cast<int>(coreAssignments.size())) {
                setThreadAffinityToCpu(coreAssignments[threadId]);
                
                // // Print confirmation from each thread
                // #pragma omp critical
                // {
                //     std::cout << "Thread " << threadId << " pinned to CPU " 
                //               << coreAssignments[threadId] << std::endl;
                // }
            }
        }
    }

    // Main algorithm
    Graph g = g0;
    std::vector<int> C(g.n);
    std::iota(C.begin(), C.end(), 0);
    H.partitions.clear();

    const int MAX_LEVELS = 10;
    const int MAX_SWEEPS = 100;
    const double EPS = 1e-6;

    for (int level = 0; level < MAX_LEVELS; level++) {
        double curQ = computeModularity(g, C);
        // std::cout << "Level " << level << " starting modularity: " << curQ << std::endl;

        // Local‐move sweeps with static scheduling
        for (int sweep = 0; sweep < MAX_SWEEPS; sweep++) {
            bool moved = parallelLocalSweepStatic(g, C);
            double newQ = computeModularity(g, C);
            // std::cout << "  Sweep " << sweep << ": modularity = " << newQ << std::endl;
            
            if (!moved || std::fabs(newQ - curQ) < EPS) {
                break;
            }
            curQ = newQ;
        }

        // Record partition at this level
        H.partitions.push_back(C);
        // std::cout << "Level " << level << " final modularity: " << curQ << std::endl;

        // Aggregate for next level
        std::vector<int> old2new;
        Graph g2 = aggregateGraph(g, C, old2new);
        if (g2.n == g.n || g2.n <= 1) {
            // std::cout << "No further aggregation possible. Stopping." << std::endl;
            break;
        }
        
        std::cout << "Aggregated graph: " << g.n << " -> " << g2.n << " nodes" << std::endl;
        g = std::move(g2);
        C.assign(g.n, 0);
        std::iota(C.begin(), C.end(), 0);
    }
}